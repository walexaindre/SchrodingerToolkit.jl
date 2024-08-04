include("M1Memory.jl")

struct M1{RealType,Grid,TKernel,ItSolver,TTime} <:
       AbstractSolverMethod{RealType}
    grid::Grid
    Kernel::TKernel
    linear_solve_params::ItSolver
    time_collection::TTime
    origin::RealType
end

time_collection(method::M1) = method.time_collection

function M1(PDE::SPDE, conf::SolverConfig,
            grid::Grid;
            linear_solver_params = IterativeLinearSolver(backend_type(conf))) where {Grid<:AbstractPDEGrid,
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
        Memory = M1Memory(ComplexCPUVector, ComplexCPUArray, PDE, conf, grid,
                          opA, opD)

        dictionary_values = Array{Kernel{SparseMatrixCSR{1,ComplexType,IntType},
                                         UniformScaling{Bool},
                                         SparseMatrixCSR{1,ComplexType,IntType}},1}(undef,
                                                                                    0)

        sizehint!(dictionary_values, length(σset) * length(TimeMultipliers))

        for (σ, βτ) in product(σset,
                               TimeMultipliers)
            opB = (4im * opA + βτ * σ * opD)
            dropzeros!(opB)
            opC = (4im * opA - βτ * σ * opD)
            dropzeros!(opC)

            Ker = Kernel(opB, I, opC)

            push!(dictionary_keys, (σ, βτ))
            push!(dictionary_values, Ker)
        end

    elseif isgpu(backend)
        Memory = M1Memory(ComplexGPUVector, ComplexGPUArray, PDE, conf, grid,
                          CuSparseMatrixCSR(opA), CuSparseMatrixCSR(opD);
                          energy_solver_params = linear_solver_params)

        dictionary_values = Array{Kernel{CuSparseMatrixCSR{ComplexType,Int32},
                                         UniformScaling{Bool},
                                         CuSparseMatrixCSR{ComplexType,Int32}},1}(undef,
                                                                                  0)
        sizehint!(dictionary_values, length(σset) * length(TimeMultipliers))

        for (σ, βτ) in product(σset,
                               TimeMultipliers)
            opB = (4im * opA + βτ * σ * opD)
            dropzeros!(opB)
            opC = (4im * opA - βτ * σ * opD)
            dropzeros!(opC)

            Ker = Kernel(CuSparseMatrixCSR(opB), I,
                         CuSparseMatrixCSR(opC))

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

    Meth = M1(grid, KernelDict, linear_solver_params, time_substeps, zero(FloatType))
    #Return must be Method, Memory, Stats
    Meth, Memory, Stats
end




function step!(method::M1, memory, stats, PDE, conf::SolverConfig)
 

    start_timer = time()

    σ_forward = get_σ(PDE)
    σ_backward = reverse(σ_forward)

    for τ in time_collection(method)
        #Forward
        for (component_index, σ) in enumerate(σ_forward)
            update_component!(PDE, Method, Memory, Stats, Style, component_index, τ,
                              σ)
        end
        #Backward

        for (component_index, σ) in zip(length(σ_backward):-1:1, σ_backward)
            update_component!(PDE, Method, Memory, Stats, Style, component_index, τ,
                              σ)
        end
    end

    power_per_component = system_power(memory, grid)
    energy = system_energy(memory, PDE, grid)

    work_timer = time() - start_timer
    
    update_stats!(stats, work_timer, power_per_component, energy)
    work_timer
end

export M1