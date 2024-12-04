#include("examples/PDE_ex2D_1.jl")
include("examples/PDE_ex2D_NR.jl")
using LinearAlgebra
using GLMakie
using CUDA
CUDA.allowscalar(false)
CUDA.device!(1)

backend = GPUBackend{Int64,Float64}

hxhy = 0.5
TF = 0.5
τst = 0.01

space_refinement(idx) = (hxhy/(2^idx), hxhy/(2^idx))
time_refinement(idx) = τst/(2^idx)

PDET = SchrodingerPDEPolynomial((Ω,Ω),(C1,C2),F,N,TF)
Params = SolverParameters(backend, 2, (:ord2, :ord2), :tord2_1_1)
Grid1 = PeriodicGrid(backend, PDE, time_refinement(0), space_refinement(0))
Grid2 = PeriodicGrid(backend, PDE, time_refinement(1), space_refinement(1))
Grid3 = PeriodicGrid(backend, PDE, time_refinement(2), space_refinement(2))
Grid4 = PeriodicGrid(backend, PDE, time_refinement(3), space_refinement(3))
Grid5 = PeriodicGrid(backend, PDE, time_refinement(4), space_refinement(4))

Method1, Memory1, Stats1 = M4(PDET, Params, Grid1)
Problem1 = SchrodingerProblem(Method1, Memory1, Stats1, PDET, Params)

Method2, Memory2, Stats2 = M4(PDET, Params, Grid2)
Problem2 = SchrodingerProblem(Method2, Memory2, Stats2, PDET, Params)

Method3, Memory3, Stats3 = M1(PDET, Params, Grid3)
Problem3 = SchrodingerProblem(Method3, Memory3, Stats3, PDET, Params)

Method4, Memory4, Stats4 = M4(PDET, Params, Grid4)
Problem4 = SchrodingerProblem(Method4, Memory4, Stats4, PDET, Params)

Method5, Memory5, Stats5 = M4(PDET, Params, Grid5)
Problem5 = SchrodingerProblem(Method5, Memory5, Stats5, PDET, Params)

ctrf=1
for i in 0:time_refinement(0):TF
    if ctrf%200 == 0
        println("Time: ", i)
    end 
    ctrf+=1
    step!(Problem1)
end
println("Finished 1")
ctrf=1
for i in 0:time_refinement(1):TF
    if ctrf%200 == 0
        println("Time: ", i)
    end 
    ctrf+=1
    step!(Problem2)
end
println("Finished 2")
ctrf=1
for i in 0:time_refinement(2):TF
    if ctrf%200 == 0
        println("Time: ", i)
    end 
    ctrf+=1
    step!(Problem3)
end
ctrf=1
println("Finished 3")
for i in 0:time_refinement(3):TF
    if ctrf%200 == 0
        println("Time: ", i)
    end 
    ctrf+=1
    step!(Problem4)
end

println("Finished 4")
ctrf=0
for i in 0:time_refinement(4):TF
    if ctrf%200 == 0
        println("Time: ", i)
    end 
    ctrf+=1
    step!(Problem5)
end

println("Finished 5")


#Convergence Error Calculation


#Calculating convergence error

function L2(memory,grid)
    curr= current_state!(memory)

    refsize = size(Grid5)
    currsize = size(grid)

    ratio = (div(refsize[1],currsize[1]), div(refsize[2] , currsize[2]))

    @show currsize
    @show refsize
    @show ratio
    @show size(curr)

    cpu_copy_curr = curr|>Array
    cpu_copy_ref = current_state!(Memory5)|>Array

    cpu_copy_curr1 = reshape(cpu_copy_curr[:,1], currsize)
    cpu_copy_ref1 = reshape(cpu_copy_ref[:,1], refsize)
    
    cpu_copy_curr2 = reshape(cpu_copy_curr[:,2], currsize)
    cpu_copy_ref2 = reshape(cpu_copy_ref[:,2], refsize)

    cpu_ref_down_scale1 = cpu_copy_ref1[1:ratio[1]:end, 1:ratio[2]:end]
    cpu_ref_down_scale2 = cpu_copy_ref2[1:ratio[1]:end, 1:ratio[2]:end]

    norm( cpu_ref_down_scale1 - cpu_copy_curr1), norm(cpu_ref_down_scale2 - cpu_copy_curr2)

end


function L2exact(memory,grid,time)
    curr = current_state!(memory)|>Array
    g1 = exp.(im*(2.0*grid[:,1]+    grid[:,2].-4.5*time))
    g2 = exp.(im*(    grid[:,1]+2.0*grid[:,2].-4.5*time))

    println("Norms: ", norm(curr[:,1]), " ",norm(curr[:,2]))
    println("Norms: ", norm(g1), " ",norm(g2))
    println("Diff: ", (curr[:,1]-g1)[1], " ",(curr[:,2]-g2)[1])
    println("Theoretical: ",g1[1], " ",g2[1])
    println("Value: ",curr[1], " ",curr[2])
    norm(curr[:,1]-g1),norm(curr[:,2]-g2)
end

function Lorder(memory_low,grid_low,memory_high,grid_high,tlow,thigh,time)
    low = L2exact(memory_low,grid_low,time)
    high = L2exact(memory_high,grid_high,time)

    order = log.(low./high)./log(tlow/thigh)
end


order(low,high,tlow,thigh) =  log.(low./high)./log(tlow/thigh)
