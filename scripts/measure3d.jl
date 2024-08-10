include("./examples/3D/PDE_ex3D_1.jl")

using GLMakie
using CUDA
CUDA.allowscalar(false)
CUDA.device!(1)
npoints = [200,400,600,800,1000,1200]
grid_size = map(x->16/x, npoints)
backend = GPUBackend{Int64,Float64}

for gsize in grid_size
    Params = SolverParameters(backend, 2, (:ord4, :ord4, :ord4), :tord2_1_1)
    Grid = PeriodicGrid(backend,PDE,0.01,(gsize,gsize,gsize))
    Method, Memory, Stats = M1(PDE,Params,Grid)
    Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)
    for j in 1:10
        step!(Problem)
    end
end
