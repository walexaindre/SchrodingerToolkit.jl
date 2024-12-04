include("M4Memory.jl")

struct M4{RealType,Grid,TKernel,TTime,StoppingCriterion} <:
       AbstractSolverMethod{RealType}
    grid::Grid
    Kernel::TKernel
    time_collection::TTime
    stopping_criteria::StoppingCriterion
    assembly_time::RealType
end

time_collection(method::M4) = method.time_collection
assembly_time(method::M4) = method.assembly_time
stopping_criteria(method::M4) = method.stopping_criteria

function M4(PDE::SPDE, conf::SolverConfig,
            grid::Grid;
            stopping_criteria = NormBased(backend_type(conf))) where {Grid<:AbstractPDEGrid,
                                                                      SPDE<:SchrodingerPDE}
    start_time = time()

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
        Memory = M4Memory(ComplexCPUVector, ComplexCPUArray, PDE, conf, grid,
                          opA, lu(opA), opD)

        dictionary_values = Vector{Kernel{SparseArrays.UMFPACK.UmfpackLU{ComplexType,
                                                                         IntType},
                                          UniformScaling{Bool},
                                          SparseMatrixCSC{ComplexType,IntType}}}(undef,
                                                                               0)

        sizehint!(dictionary_values, length(σset) * length(TimeMultipliers))

        for (σ, βτ) in product(σset,
                               TimeMultipliers)
            opB = (4im * opA + βτ * σ * opD)
            dropzeros!(opB)
            opC = (4im * opA - βτ * σ * opD)
            dropzeros!(opC)

            Ker = Kernel(lu(opB), I, opC)

            push!(dictionary_keys, (σ, βτ))
            push!(dictionary_values, Ker)
        end

    elseif isgpu(backend)
        gpuopA = CuSparseMatrixCSR(opA)
        Memory = M4Memory(ComplexGPUVector, ComplexGPUArray, PDE, conf, grid,
                          gpuopA, lu(gpuopA), CuSparseMatrixCSR(opD))

        dictionary_values = Vector{Kernel{CudssSolver{ComplexType},
                                          UniformScaling{Bool},
                                          CuSparseMatrixCSR{ComplexType,Int32}}}(undef,
                                                                                 0)
        sizehint!(dictionary_values, length(σset) * length(TimeMultipliers))

        for (σ, βτ) in product(σset,
                               TimeMultipliers)
            opB = (4im * opA + βτ * σ * opD)
            dropzeros!(opB)
            opC = (4im * opA - βτ * σ * opD)
            dropzeros!(opC)

            Ker = Kernel(lu(CuSparseMatrixCSR(opB)), I,
                         CuSparseMatrixCSR(opC))

            push!(dictionary_keys, (σ, βτ))
            push!(dictionary_values, Ker)
        end
    end

    KernelDict = Dictionary(dictionary_keys, dictionary_values)

    cstate = current_state!(Memory)
    evaluate_ψ!(PDE, grid, cstate)

    power_startup = system_mass(Memory, grid)
    energy_startup = system_energy(Memory, PDE, grid)

    startup_stats!(Stats, power_startup, energy_startup)

    assembly_time = time() - start_time

    Meth = M4(grid, KernelDict, time_substeps, stopping_criteria,
              FloatType(assembly_time))
    #Return must be Method, Memory, Stats
    Meth, Memory, Stats
end

@inline function update_component!(method::M4, memory, stats, PDE, τ, σ, component_index)
    current_state = current_state!(memory)
    grid_measure = sqrt(measure(method.grid))
    solved = false
    step_count = 0

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

    #End of Memory temporary arrays

    #Method operators
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
    for l in 1:get_max_iterations(stopping_criteria_m1)
        @. current_state_abs2 = abs2(current_state)
        @. temporary_abs2 = abs2(zₗ)
        @. stage1 = zₗ + ψ
        stage2 .= N(current_state_abs2, temporary_abs2, component_index)
        @. b_temp = stage1 * stage2
        mul!(stage1, opA, b_temp)
        @. b_temp = τ * stage1 + b0_temp
        sol_etime = CUDA.@elapsed ldiv!(stage1, opB, b_temp)
        update_solver_info!(stats, sol_etime, 0)

        #NormBased
        copy!(stage2, zₗ)
        copy!(zₗ, stage1)
        @. stage1 = stage2 - zₗ

        znorm = grid_measure * norm(stage1)
        solved = znorm <=
                 get_atol(stopping_criteria_m1) + get_rtol(stopping_criteria_m1) * znorm
        if solved
            step_count = l
            break
        end
    end
    copy!(ψ, zₗ)

    if !solved
        @warn "Convergence not reached in $(get_max_iterations(stopping_criteria_m1)) iterations..."
    end
    step_count
end

function step!(method::M4, memory, stats, PDE, conf::SolverConfig)
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

export M4, assembly_time