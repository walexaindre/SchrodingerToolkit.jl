include("../examples/PDE_ex2D_1.jl")
using GLMakie
using CUDA
CUDA.allowscalar(false)
CUDA.device!(1)

backend = GPUBackend{Int64,Float64}

Params = SolverParameters(backend, 2, (:ord4, :ord4), :tord2_1_1; stats=ConfigRuntimeStatsOptions(2; log_component_update_steps=true) )
Grid = PeriodicGrid(backend, PDE, 0.01, (0.16, 0.16))
Method, Memory, Stats = M1(PDE, Params, Grid)

Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)

for i in 1:340
    step!(Problem)
end

f,ax = plot_systemnd(Memory, Grid, 1)
f,ax = plot_solver_iterations(Stats, Grid)
f,ax = plot_mass_per_component(Stats, Grid,1)
f,ax = plot_preprocessing_time(Stats, Grid)
f,ax = plot_absolute_error_mass_per_component(Stats, Grid,1)
f,ax = plot_execution_time(Stats, Grid)
f,ax = plot_total_mass(Stats, Grid)
f,ax = plot_absolute_error_total_mass(Stats, Grid)
f,ax = plot_absolute_error_system_energy(Stats, Grid)
f,ax = plot_solver_time(Stats, Grid)
f,ax = plot_component_update_steps(Stats,Grid)

save("a.png",f)