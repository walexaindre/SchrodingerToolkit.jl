include("examples/PDE_ex2D_1.jl")

backend = GPUBackend{Int64,Float64}

Params = SolverParameters(backend, 2, (:ord4, :ord4), :tord2_1_1)
Grid = PeriodicGrid(backend,PDE,0.01,(0.16,0.16))

Stats = initialize_stats(Vector{Float64},Vector{Int}, 3, 2,100)

