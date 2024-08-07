include("examples/PDE_ex2D_1Gross.jl")
using GLMakie
using CUDA
CUDA.allowscalar(false)
CUDA.device!(1)

backend = GPUBackend{Int64,Float64}

Params = SolverParameters(backend, 2, (:ord4, :ord4), :tord2_1_1)
Grid = PeriodicGrid(backend, PDE, 0.01, (0.16, 0.16))
Method, Memory, Stats = M2(PDE, Params, Grid)

Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)
g = Figure(; size = (800, 600), fontsize = 25)
ax3 = Axis3(g[1, 1])
ax3.zlabel = ""
ax3.yticks = WilkinsonTicks(6; k_min = 5)
ax3.xticks = WilkinsonTicks(6; k_min = 5)
for i in 1:250
    for i in 1:50
        step!(Problem)
    end
    empty!(ax3)

    ax3.title = "iteration $(current_iteration(Stats)) time $(round(Grid.τ*current_iteration(Stats),digits=2))"
    systemnd!(ax3, Memory, Grid, 2)
    sleep(0.5)
end



systemnd!(ax3, Memory, Grid, 1)

g

ax4 = Axis3(g[1, 2])
ax4.zlabel = ""
ax4.yticks = WilkinsonTicks(6; k_min = 5)
ax4.xticks = WilkinsonTicks(6; k_min = 5)


for i in 1:250
    for i in 1:15
        step!(Problem)
    end
    empty!(ax3)
    empty!(ax4)

    ax3.title = "iteration $(current_iteration(Stats)) time $(round(Grid.τ*current_iteration(Stats),digits=2)) c: 1"
    ax4.title = "iteration $(current_iteration(Stats)) time $(round(Grid.τ*current_iteration(Stats),digits=2)) c: 2"

    systemnd!(ax3, Memory, Grid, 1)
    systemnd!(ax4, Memory, Grid, 2)
    sleep(0.5)
end









for i in 1:200
    step!(Problem)
end