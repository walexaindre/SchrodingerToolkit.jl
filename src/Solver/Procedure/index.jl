#Crank Nicolson
include("CrankNicolson/M1.jl")
include("CrankNicolson/M4.jl")

#Directional Splitting
include("DirectionalDecomposition/M2.jl")
include("DirectionalDecomposition/M3.jl")

#Explicit Timestepping
include("ExplicitTimestepping/M5.jl")

#Explict Timestepping with Explicit discretization
include("ExplicitTimesteppingExplicitDiscretization/M6.jl")