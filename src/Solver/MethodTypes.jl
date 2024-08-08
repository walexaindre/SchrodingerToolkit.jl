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

"""
    AbstractSolverMethod{RealType}

    This is the abstract type for all solver methods. It is parametrized by the type of the real numbers used in the solver.

    # How to define a new solver method
    To define a new solver method, you need to define a new struct that is a subtype of this abstract type. The struct must have the following fields:
    - `Grid`: The grid on which the solver will operate.

    # Expected initialization
    The solver method must be initialized with the following arguments:
    - `config`: The basic configuration of the solver.
    - `grid`: The grid on which the solver will operate.

    # Required methods
    The solver method must implement the following methods:
    - `step!`: Perform a single step of the solver.

    # Methods related to statistics
    [TODO]
"""
abstract type AbstractSolverMethod{RealType} end

struct Kernel{LinOp1,LinOp2,LinOp3}
    opB::LinOp1
    preB::LinOp2
    opC::LinOp3
end

abstract type AbstractMemory end


abstract type MethodProperty end

struct IsDirect<:MethodProperty end
struct IsIterative<:MethodProperty end