include("./examples/3D/PDE_ex3D_1.jl")

using GLMakie
using CUDA
using JSON
CUDA.allowscalar(false)
CUDA.device!(1)
npoints = [30,40,50,80,130]
grid_size = map(x->16/x, npoints)
backend = GPUBackend{Int64,Float64}
log_freq = ConfigRuntimeStatsOptions(Int(1);log_solver_info=true, log_component_update_steps=true)
for gsize in grid_size
    Params = SolverParameters(backend, 2, (:ord4, :ord4, :ord4), :tord2_1_1;stats = log_freq)
    Grid = PeriodicGrid(backend,PDE,0.01,(gsize,gsize,gsize))
    Method, Memory, Stats = M1(PDE,Params,Grid)
    @show size(current_state!(Memory))

    Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)
    for j in 1:18
        step!(Problem)
    end

    svec1 = step_time(Stats)
    svec2 = solver_time(Stats)
    sint1 = component_update_calls(Stats)
    svec3 = component_update_steps(Stats)
    total_component_update_steps = sum(svec3)
    a_time = assembly_time(Method)
    iterations = solver_iterations(Stats)
    iterations_avg = sum(iterations) / (total_component_update_steps)

    size(current_state!(Memory))

    println("Grid size: ", gsize)
    println("Backend: ", backend)
    println("Method: ", typeof(Method))
    println("Step time avg: ", sum(svec1) / length(svec1))
    println("Solver time avg: ", sum(svec2) / length(svec2))
    println("Solver iterations avg: ", iterations_avg)
    println("Component update steps avg: ", sum(svec3) / sint1)
    println("Component update calls: ", sint1)
    println("Asembly time: ", a_time)

    def = Dict("Grid size" => round(gsize,digits=2), "Backend" => backend, "Method" => "$(Base.typename(typeof(Method)))",
               "Step time avg" => sum(svec1) / length(svec1),
               "Solver time avg" => sum(svec2) / length(svec2),
               "Solver iterations avg" => iterations_avg,
               "Component update steps avg" => sum(svec3) / sint1,
               "Component update calls" => sint1, "Asembly time" => a_time)


    open("$(Base.typename(typeof(Method)))_$(backend)_$(Int(16/gsize))_$(ndims(PDE))D.json", "a") do f
        JSON.print(f, def)
    end
end
