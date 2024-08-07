include("examples/PDE_ex2D_1.jl")
using GLMakie
using CUDA
CUDA.allowscalar(false)
CUDA.device!(1)

backend = GPUBackend{Int64,Float64}

Params = SolverParameters(backend, 2, (:ord4, :ord4), :tord2_1_1)
Grid = PeriodicGrid(backend,PDE,0.01,(0.16,0.16))
Method, Memory, Stats = M1(PDE,Params,Grid)
Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)

for i in 1:200
    step!(Problem)
end