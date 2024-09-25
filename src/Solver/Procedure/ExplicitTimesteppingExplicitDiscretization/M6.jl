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

"""
    Fully explicit method    
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













    assembly_time = time() - start_time
end