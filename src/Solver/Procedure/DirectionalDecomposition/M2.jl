include("M2Memory.jl")

struct M2{RealType,Grid,TKernel,ItSolver,TTime,StoppingCriterion} <:
       AbstractSolverMethod{RealType}
    grid::Grid
    Kernel::TKernel

    linear_solve_params::ItSolver
    time_collection::TTime
    stopping_criteria::StoppingCriterion
    assembly_time::RealType
end

stopping_criteria(method::M2) = method.stopping_criteria
time_collection(method::M2) = method.time_collection
linear_solve_params(method::M2) = method.linear_solve_params
assembly_time(method::M2) = method.assembly_time

function opLn(PDE, grid, opA, opD, σ)
    points = collect_points(grid)
    Vᵢ = trapping_potential(PDE, 1) #[Todo] . [Fix this]
    opA * spdiagm(Vᵢ(points)) - σ * opD
end

function M2(PDE::SPDE, conf::SolverConfig,
            grid::Grid;
            linear_solver_params = IterativeLinearSolver(backend_type(conf)),
            stopping_criteria = NormBased(backend_type(conf))) where {Grid<:AbstractPDEGrid,
                                                                      SPDE<:SchrodingerPDE}
    backend = backend_type(conf)
    FloatType = BackendReal(backend)
    ComplexType = Complex{FloatType}
    IntType = BackendInt(backend)

    start_time = time()

    FloatCPUVector = Vector{FloatType}
    IntCPUVector = Vector{IntType}
    ComplexCPUVector = Vector{Complex{FloatType}}
    ComplexCPUArray = Array{Complex{FloatType},2}
    #FloatGPUVector = CuArray{FloatType}
    #IntGPUVector = CuArray{IntType}
    ComplexGPUVector = CuVector{Complex{FloatType}}
    ComplexGPUArray = CuArray{Complex{FloatType},2}

    Stats = initialize_stats(FloatCPUVector, IntCPUVector, PDE, conf, grid)

    AI, AJ, AV = get_A_format_COO(FloatType, PeriodicAbstractGrid(grid),
                                  space_order(conf))

    DI, DJ, DV = get_D_format_COO(FloatType, grid,
                                  space_order(conf))

    opA = sparse(AI, AJ, ComplexCPUVector(AV))

    opD = sparse(DI, DJ, ComplexCPUVector(DV))

    time_comp = get_time_composition(time_order(conf))

    σset = Set(get_σ(PDE))
    TimeMultipliers = grid.τ * coefficients(time_comp)

    time_substeps = Tuple(grid.τ * collect(time_comp))

    dictionary_keys = Array{Tuple{FloatType,FloatType},1}(undef, 0)

    sizehint!(dictionary_keys, length(σset) * length(TimeMultipliers))

    if iscpu(backend)
        Memory = M2Memory(ComplexCPUVector, ComplexCPUArray, PDE, conf, grid,
                          opA, opD)

        dictionary_values = Array{Kernel{SparseMatrixCSR{1,ComplexType,IntType},
                                         UniformScaling{Bool},
                                         UniformScaling{Bool}},1}(undef,
                                                                  0)

        sizehint!(dictionary_values, length(σset) * length(TimeMultipliers))

        for (σ, βτ) in product(σset,
                               TimeMultipliers)
            Ln = opLn(PDE, grid, opA, opD, σ)

            opB = opA + im * (βτ / 4) * Ln

            Ker = Kernel(opB, I,
                         I)

            push!(dictionary_keys, (σ, βτ))
            push!(dictionary_values, Ker)
        end

    elseif isgpu(backend)
        Memory = M2Memory(ComplexGPUVector, ComplexGPUArray, PDE, conf, grid,
                          CuSparseMatrixCSR(opA), CuSparseMatrixCSR(opD);
                          energy_solver_params = linear_solver_params)

        dictionary_values = Array{Kernel{CuSparseMatrixCSR{ComplexType,Int32},
                                         UniformScaling{Bool},
                                         UniformScaling{Bool}},1}(undef,
                                                                  0)
        sizehint!(dictionary_values, length(σset) * length(TimeMultipliers))

        for (σ, βτ) in product(σset,
                               TimeMultipliers)
            Ln = opLn(PDE, grid, opA, opD, σ)

            opB = opA + im * (βτ / 4) * Ln

            Ker = Kernel(CuSparseMatrixCSR(opB), I,
                         I)

            push!(dictionary_keys, (σ, βτ))
            push!(dictionary_values, Ker)
        end
    end

    KernelDict = Dictionary(dictionary_keys, dictionary_values)

    cstate = current_state!(Memory)
    evaluate_ψ!(PDE, grid, cstate)

    power_startup = system_power(Memory, grid)
    energy_startup = system_energy(Memory, PDE, grid)

    startup_stats!(Stats, power_startup, energy_startup)

    assembly_time = time() - start_time

    Meth = M2(grid, KernelDict, linear_solver_params, time_substeps, stopping_criteria,
              FloatType(assembly_time))
    #Return must be Method, Memory, Stats
    Meth, Memory, Stats
end

@inline function update_component!(method::M2, memory, stats, PDE, τ, σ, component_index)
    iτhalf = (τ / 4) * im
    Γ = junction_coefficient(PDE)

    current_state = current_state!(memory)
    grid_measure = sqrt(measure(method.grid))
    solved = false

    ψ = view(current_state, :, component_index)
    zₗ = memory.component_temp

    #dst src
    copy!(zₗ, ψ)

    current_state_abs2 = memory.current_state_abs2
    temporary_abs2 = memory.temp_state_abs2

    stage1 = memory.stage1 #Is assumed that the norm of stage1 is the norm of the difference
    stage2 = memory.stage2

    b0_temp = memory.b0_temp
    b_temp = memory.b_temp

    SolverMem = solver_memory!(memory)
    #End of Memory temporary arrays

    #Method operators
    solver_params = linear_solve_params(method)
    stopping_criteria_m2 = stopping_criteria(method)

    Kernel = method.Kernel[(σ, τ)]

    opB = Kernel.opB
    opA = memory.opA
    #End of Method operators

    #N optimized in PDE for easy calculations
    N = get_optimized(PDE)
    #End of N optimized

    #mul!(b0_temp, opC, ψ)
    @. current_state_abs2 = abs2(current_state)
    junction!(PDE, memory, b0_temp, component_index) # evaluation of the junction Jⁿ
    b0_temp .*= Γ # Γ * Jⁿ

    exit_iterations = 0
    for l in 1:get_max_iterations(stopping_criteria_m2)
        @. temporary_abs2 = abs2(zₗ)
        stage2 .= zₗ + ψ # zₗ + ψ
        stage2 .*= 0.5
        stage1 .= N(current_state_abs2, temporary_abs2, component_index) # evaluation of N
        @. b_temp = stage1 * stage2 # N ⊙ (zₗ+ψ)

        b_temp .+= b0_temp # N ⊙ (zₗ+ψ) + Γ * Jⁿ
        b_temp .*= -iτhalf # -iτ/4 * (N ⊙ (zₗ+ψ) + Γ * Jⁿ)

        b_temp .+= ψ # ψ - iτ/4 * (N ⊙ (zₗ+ψ) + Γ * Jⁿ)

        mul!(stage1, opA, b_temp) # A * (ψ - iτ/4 * (N ⊙ (zₗ+ψ) + Γ * Jⁿ))

        gmres!(SolverMem, opB, stage1; atol = get_atol(solver_params),
               rtol = get_rtol(solver_params),
               itmax = get_max_iterations(solver_params))
        update_solver_info!(stats, SolverMem.stats.timer, SolverMem.stats.niter)

        #NormBased
        copy!(stage2, zₗ) #dst src
        copy!(zₗ, SolverMem.x)#dst src
        zₗ .= 2 * SolverMem.x - ψ
        @. stage1 = stage2 - zₗ

        znorm = grid_measure * norm(stage1)
        solved = znorm <=
                 get_atol(stopping_criteria_m2) + get_rtol(stopping_criteria_m2) * znorm
        if solved
            exit_iterations = l
            break
        end
    end
    copy!(ψ, zₗ)

    if !solved
        @warn "Convergence not reached in $(get_max_iterations(stopping_criteria_m2)) iterations..."
    end
    exit_iterations
end

function step!(method::M2, memory, stats, PDE, conf::SolverConfig)
    start_timer = time()
    grid = method.grid

    σ_forward = get_σ(PDE)
    σ_backward = reverse(σ_forward)

    for τ in time_collection(method)
        #Forward
        for (component_index, σ) in enumerate(σ_forward)
            steps = update_component!(method, memory, stats, PDE, τ, σ, component_index)
            update_component_update_steps!(stats, steps)
        end

        #Backward
        for (component_index, σ) in zip(length(σ_backward):-1:1, σ_backward)
            steps = update_component!(method, memory, stats, PDE, τ, σ, component_index)
            update_component_update_steps!(stats, steps)
        end
    end
    work_timer = time() - start_timer

    update_stats!(stats, memory, grid, PDE,
                  work_timer)
    work_timer
end

export M2