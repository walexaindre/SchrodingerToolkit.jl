function DefaultSolver(::Type{ComputeBackend},
                       dim::IntType, spaceorder = ntuple(Returns(:ord2), dim),
                       timeorder = :tord2_1_1) where {IntType,FloatType,ComplexType,
                                                      VectorType,
                                                      VectorComplexType,MatrixType,
                                                      MatrixComplexType,
                                                      ComputeBackend<:AbstractBackend{IntType,
                                                                                      FloatType,
                                                                                      ComplexType,
                                                                                      VectorType,
                                                                                      VectorComplexType,
                                                                                      MatrixType,
                                                                                      MatrixComplexType}}
    SolverConfig(timeorder,
                 spaceorder,
                 ComputeBackend,
                 zero(IntType),
                 false,
                 "/data",
                 NormBased(ComputeBackend),
                 (true, IntType(5)),
                 (false, zero(IntType), zero(IntType)))
end

"""

    With this function you can advance one step in time.

"""
function step! end