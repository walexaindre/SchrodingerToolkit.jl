include("examples/PDE_ex2D_1.jl")

using CUDA
using CUDA.CUSPARSE
using CUDSS

using LinearAlgebra

using SparseArrays
device!(1)
backend = GPUBackend{Int64,Float64}


Params = SolverParameters(backend, 2, (:ord4, :ord4), :tord2_1_1)
Grid = PeriodicGrid(backend,PDE,0.01,(0.16,0.16))
Method, Memory, Stats = M1(PDE,Params,Grid)
Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)

opB = B.opB

Fr = lu(opB)

a = rand(ComplexF64, 10000)|>CuVector

CUDA.@time Fr \ a


opBcpu = opB|>SparseMatrixCSC

Frcpu = lu(opBcpu)

acpu = rand(ComplexF64, 10000)

@time Frcpu \ acpu