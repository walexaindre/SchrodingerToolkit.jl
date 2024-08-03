struct PaulMethod{RealType,Grid,TKernel,ItSolver,ItSolver2,TTime} <: AbstractSolverMethod{RealType}
    Mesh::Grid
    Kernel::TKernel
    linear_solve_params::ItSolver
    preconditioner_drop_tol::RealType
    preconditioner_solve_params::ItSolver2
    time_collection::TTime
end