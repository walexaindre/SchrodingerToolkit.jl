#Types => Here we include all the types that are used in the package

#Core
include("Core/Offset/OffsetTypes.jl")
include("Core/GridModel/GridModelTypes.jl")
include("Core/GridPDE/GridPDETypes.jl")
include("Core/GridDiscretization/GridDiscretizationTypes.jl")
include("Core/TimeComposition/TimeCompositionTypes.jl")
include("Core/ModelPDE/ModelPDETypes.jl")
include("Core/RuntimeStatistics/RuntimeStatisticsTypes.jl")
include("Core/Diagrams/DiagramsRecipes.jl")


#Solver
include("Solver/MethodTypes.jl")
include("Solver/SolverTypes.jl")
include("Solver/Backend/ComputeBackend.jl")


#Methods and declarations => Here we include all the methods and declarations that are used in the package

#Core
include("Core/VecOps/VecOps.jl")
include("Core/Offset/Offset.jl")
include("Core/GridModel/GridModel.jl")
include("Core/GridPDE/GridPDE.jl")
include("Core/ModelPDE/ModelPDE.jl")
include("Core/TimeComposition/TimeComposition.jl")
include("Core/RuntimeStatistics/RuntimeStatistics.jl")
include("Core/Preconditioner/Preconditioner.jl")
include("Core/Diagrams/Diagrams.jl")
include("Core/TupleSort/TupleSort.jl")
include("Core/GridDiscretization/GridDiscretization.jl")



#Solver
include("Solver/Method.jl")
include("Solver/Solver.jl")


#PDE Models
include("Solver/Procedure/index.jl")