IterativeLinearSolver(::Type{ComputeBackend}) where {RealType,IntType,
ComputeBackend<:AbstractBackend{IntType,RealType}} = IterativeLinearSolver(700 *
                                                                            eps(RealType),
                                                                            700 *
                                                                            eps(RealType),
                                                                            IntType(10000))

NormBased(::Type{ComputeBackend}) where {RealType,IntType,
ComputeBackend<:AbstractBackend{IntType,RealType}} = NormBased(700 * eps(RealType),
                                                                700 * eps(RealType),
                                                                IntType(10000))

@inline get_atol(param::NormBased{IntType,RealType}) where {IntType,RealType} = param.atol
@inline get_rtol(param::NormBased{IntType,RealType}) where {IntType,RealType} = param.rtol
@inline get_max_iterations(param::NormBased{IntType,RealType}) where {IntType,RealType} = param.max_steps

@inline get_atol(param::IterativeLinearSolver{IntType,RealType}) where {IntType,RealType} = param.atol
@inline get_rtol(param::IterativeLinearSolver{IntType,RealType}) where {IntType,RealType} = param.rtol
@inline get_max_iterations(param::IterativeLinearSolver{IntType,RealType}) where {IntType,RealType} = param.max_iterations

Base.show(io::IO, param::IterativeLinearSolver) = print(io,
                                                        "\nIterativeLinearSolver\nrtol: $(param.rtol)\natol: $(param.atol)\nmaxiter: $(param.max_iterations)")
Base.show(io::IO, param::NormBased) = print(io,
                                            "\nStoppingCriteria\nNormBased\natol: $(param.atol)\nrtol: $(param.rtol)\nmaxiter: $(param.max_steps)")
Base.show(io::IO, param::FixedSteps) = print(io,
                                             "\nStoppingCriteria\nFixedSteps\nnsteps: $(param.nsteps)")

export IterativeLinearSolver, NormBased, FixedSteps, DirectLinearSolver, AbstractMemory