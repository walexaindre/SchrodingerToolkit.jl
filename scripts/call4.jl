include("examples/3D/PDE_ex3D_1.jl")

using GLMakie


backend = CPUBackend{Int64,Float64}

Params = SolverParameters(backend, 2, (:ord4, :ord4, :ord4), :tord2_1_1)
Grid = PeriodicGrid(backend,PDE,0.01,(0.16,0.16,0.16))
Method, Memory, Stats = M1(PDE,Params,Grid)
Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)

for j in 1:10
    step!(Problem)
end
using LinearAlgebra
Fb = lu()