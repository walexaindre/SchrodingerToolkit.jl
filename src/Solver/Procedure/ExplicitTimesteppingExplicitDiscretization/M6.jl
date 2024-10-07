include("M6Memory.jl")

struct M6{RealType,Grid,TKernel,ItSolver,TTime,StoppingCriterion} <:
       AbstractSolverMethod{RealType}
    grid::Grid
    Kernel::TKernel
    linear_solve_params::ItSolver
    time_collection::TTime
    stopping_criteria::StoppingCriterion
    assembly_time::RealType
end

stopping_criteria(method::M6) = method.stopping_criteria
time_collection(method::M6) = method.time_collection
linear_solve_params(method::M6) = method.linear_solve_params
assembly_time(method::M6) = method.assembly_time

"""
    Fully explicit method with Parallel stencil
"""
function M6(PDE::SPDE, conf::SolverConfig,
            grid::Grid;
            stopping_criteria = NormBased(backend_type(conf))) where {Grid<:AbstractPDEGrid,
                                                                      SPDE<:SchrodingerPDE}
    start_time = time()

    if !is_explicit(conf)
        throw(ArgumentError("M6 is an explicit method, please use configure your parameters accordingly"))
    end

    backend = backend_type(conf)

    if isgpu(backend)
        @init_parallel_stencil()
    end

    FloatType = BackendReal(backend)
    ComplexType = Complex{FloatType}
    IntType = BackendInt(backend)

    FloatCPUVector = Vector{FloatType}
    IntCPUVector = Vector{IntType}
    ComplexCPUVector = Vector{Complex{FloatType}}
    ComplexCPUArray = Array{Complex{FloatType},2}

    ComplexGPUVector = CuVector{Complex{FloatType}}
    ComplexGPUArray = CuArray{Complex{FloatType},2}

    Stats = initialize_stats(FloatCPUVector, IntCPUVector, PDE, conf, grid)

    nelem = size(grid)



    assembly_time = time() - start_time
end

@inline function update_component!(method::M6, memory, stats, PDE, τ, σ, component_index)
    stopping_criteria_m6 = stopping_criteria(method)

    exit_iterations = 0
    for l in 1:get_max_iterations(stopping_criteria_m6)
        if solved
            exit_iterations = l
            break
        end
    end
    exit_iterations
end