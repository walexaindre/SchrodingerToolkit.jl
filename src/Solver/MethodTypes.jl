abstract type AbstractStoppingCriterion end

struct FixedSteps{IntType}<:AbstractStoppingCriterion
    nsteps::IntType
end

struct NormBased{IntType,FloatType}<:AbstractStoppingCriterion
    atol::FloatType
    rtol::FloatType
    max_steps::IntType
end

abstract type AbstractLinearSolver end

struct IterativeLinearSolver{IntType,RealType}<:AbstractLinearSolver
    rtol::RealType
    atol::RealType
    max_iterations::IntType
end

struct DirectLinearSolver<:AbstractLinearSolver
end

abstract type AbstractSolverMethod{RealType} end

struct Kernel{LinOp1,LinOp2,LinOp3}
    opB::LinOp1
    preB::LinOp2
    opC::LinOp3
end


abstract type MethodProperty end

struct IsDirect<:MethodProperty end
struct IsIterative<:MethodProperty end