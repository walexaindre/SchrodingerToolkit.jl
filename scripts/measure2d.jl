include("./examples/PDE_ex2D_1.jl")
using GLMakie
using CUDA
using JSON
CUDA.allowscalar(false)
CUDA.device!(1)
npoints = [200, 400, 600, 800, 1000, 1200]
grid_size = map(x -> 16 / x, npoints)

log_freq = ConfigRuntimeStatsOptions(Int(1);log_component_update_steps=true)

for gsize in grid_size
    backend = CPUBackend{Int64,Float64}
    Params = SolverParameters(backend, 2, (:ord4, :ord4), :tord2_1_1; stats = log_freq)
    Grid = PeriodicGrid(backend, PDE, 0.01, (gsize, gsize))
    Method, Memory, Stats = M1(PDE, Params, Grid)
    Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)
    for j in 1:18
        step!(Problem)
    end

    svec1 = step_time(Stats)
    svec2 = solver_time(Stats)
    sint1 = component_update_calls(Stats)
    svec3 = component_update_steps(Stats)
    a_time = assembly_time(Method)

    println("Grid size: ", gsize)
    println("Backend: ", backend)
    println("Method: ", typeof(Method))
    println("Step time avg: ", sum(svec1) / length(svec1))
    println("Solver time avg: ", sum(svec2) / length(svec2))
    println("Component update steps avg: ", sum(svec3) / sint1)
    println("Component update calls: ", sint1)
    println("Asembly time: ", a_time)

    def = Dict("Grid size" => gsize, "Backend" => backend, "Method" => typeof(Method),
               "Step time avg" => sum(svec1) / length(svec1),
               "Solver time avg" => sum(svec2) / length(svec2),
               "Component update steps avg" => sum(svec3) / sint1,
               "Component update calls" => sint1, "Asembly time" => a_time)

    open("M_1_$(backend)_$(round(16/gsize)).json", "a") do f
        JSON.print(f, def)
    end

end
