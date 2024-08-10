include("examples/PDE_ex2D_1.jl")
using GLMakie
using CUDA
CUDA.allowscalar(false)
CUDA.device!(1)

backend = GPUBackend{Int64,Float64}

Params = SolverParameters(backend, 2, (:ord4, :ord4), :tord2_1_1)
Grid = PeriodicGrid(backend,PDE,0.01,(0.016,0.016))
Method, Memory, Stats = M4(PDE,Params,Grid)
Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)

fig = Figure()
fig3daxis = Axis3(fig[1, 1])
hidespines!(fig3daxis)
hidedecorations!(fig3daxis)
intxy = -8..8

obs = Observable(reshape(abs2.(current_state(Memory)[:,1]|>Array),(1600,1600)))

surface!(fig3daxis, intxy, intxy, obs, colormap = :balance)
GLMakie.record(fig,"2Dgpu.gif",1:300;framerate=15) do t
    obs[] = reshape(abs2.(current_state!(Memory)[:,1]|>Array),(1600,1600))
    step!(Problem)
end

for j in 1:400
    step!(Problem)
end 