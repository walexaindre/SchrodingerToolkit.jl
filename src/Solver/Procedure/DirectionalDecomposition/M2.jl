include("M2Memory.jl")

struct M2{RealType,Grid,TKernel,ItSolver,TTime,StoppingCriterion} <:
       AbstractSolverMethod{RealType}
    grid::Grid
    Kernel::TKernel

    linear_solve_params::ItSolver
    time_collection::TTime
    stopping_criteria::StoppingCriterion
    origin::RealType
end

stopping_criteria(method::M2) = method.stopping_criteria
time_collection(method::M2) = method.time_collection
linear_solve_params(method::M2) = method.linear_solve_params

struct M2Kernel{KernelType}
    LKernel::KernelType
end

function M2Kernel(PDE, curr_state, opA, opD)
    ncomp = ncomponents(PDE)
    points = typeof(curr_state)(collect_points(grid))

    for comp_i in 1:ncomp
        Vᵢ = trapping_potential(PDE, comp_i)
        σ = get_σ(PDE, comp_i)
        opA * Vᵢ(points) - σ * opD
    end
end

function SparseArrays.spdiagm(v::CuArray{Tv}) where {Tv}
    nzVal = v
    N = Int32(length(v))

    colPtr = CuArray(one(Int32):(N + one(Int32)))
    rowVal = CuArray(one(Int32):N)
    dims = (N, N)
    CuSparseMatrixCSC(colPtr, rowVal, nzVal, dims)
end

Base.display(spdiagm([1, 2, 3]))

function M2(PDE::SPDE, conf::SolverConfig,
            grid::Grid;
            linear_solver_params = IterativeLinearSolver(backend_type(conf)),
            stopping_criteria = NormBased(backend_type(conf))) where {Grid<:AbstractPDEGrid,
                                                                      SPDE<:SchrodingerPDE}
    backend = backend_type(conf)
    FloatType = BackendReal(backend)
    ComplexType = Complex{FloatType}
    IntType = BackendInt(backend)

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

    opA = sparse(AI, AJ, ComplexCPUVector(AV))

    opD = sparse(DI, DJ, ComplexCPUVector(DV))

    time_comp = get_time_composition(time_order(conf))

    σset = Set(get_σ(PDE))
    TimeMultipliers = grid.τ * coefficients(time_comp)

    time_substeps = Tuple(grid.τ * collect(time_comp))

    dictionary_keys = Array{Tuple{FloatType,FloatType},1}(undef, 0)

    sizehint!(dictionary_keys, length(σset) * length(TimeMultipliers))
end

@inline function update_component!(method::M2, memory, stats, PDE, τ, σ, component_index)
    iτhalf = (τ / 2) * im
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
    stopping_criteria_m1 = stopping_criteria(method)

    Kernel = method.Kernel[(σ, τ)]

    opC = Kernel.opC
    opB = Kernel.opB
    opA = memory.opA
    #End of Method operators

    #N optimized in PDE for easy calculations
    N = get_optimized(PDE)
    #End of N optimized

    mul!(b0_temp, opC, ψ)
    for _ in 1:get_max_iterations(stopping_criteria_m1)
        @. current_state_abs2 = abs2(current_state)
        @. temporary_abs2 = abs2(zₗ)

        @. stage1 = N(current_state_abs2, temporary_abs2, component_index) # evaluation of N
        @. b_temp = stage1 * zₗ # N ⊙ zₗ

        junction!(PDE, memory, stage2, component_index) # evaluation of the junction Jⁿ
        stage2 .= Γ # Γ * Jⁿ

        b_temp .+= stage2 # N ⊙ zₗ + Γ * Jⁿ
        b_temp .*= -iτhalf # -iτ/2 * (N ⊙ zₗ + Γ * Jⁿ)
        b_temp .+= ψ # ψ - iτ/2 * (N ⊙ zₗ + Γ * Jⁿ)

        mul!(stage1, opA, b_temp) # A * (ψ - iτ/2 * (N ⊙ zₗ + Γ * Jⁿ))

        gmres!(SolverMem, opB, b_temp; atol = get_atol(solver_params),
               rtol = get_rtol(solver_params),
               itmax = get_max_iterations(solver_params))
        copy!(zₗ, SolverMem.x)
        update_solver_info!(stats, SolverMem.stats.timer, SolverMem.stats.niter)

        #NormBased
        copy!(stage2, zₗ)
        copy!(zₗ, SolverMem.x)
        @. stage1 = stage2 - zₗ

        znorm = grid_measure * norm(stage1)
        solved = znorm <=
                 get_atol(stopping_criteria_m1) + get_rtol(stopping_criteria_m1) * znorm
        if solved
            break
        end
    end
    copy!(ψ, zₗ)
    nothing
end

function step!(method::M2, memory, stats, PDE, conf::SolverConfig)
    start_timer = time()
    grid = method.grid

    σ_forward = get_σ(PDE)
    σ_backward = reverse(σ_forward)

    for τ in time_collection(method)
        #Forward
        for (component_index, σ) in enumerate(σ_forward)
            update_component!(method, memory, stats, PDE, τ, σ, component_index)
        end
        #Backward

        for (component_index, σ) in zip(length(σ_backward):-1:1, σ_backward)
            update_component!(method, memory, stats, PDE, τ, σ, component_index)
        end
    end

    power_per_component = system_power(memory, grid)
    energy = system_energy(memory, PDE, grid)

    work_timer = time() - start_timer

    update_stats!(stats, work_timer, power_per_component, energy)
    work_timer
end

export M2