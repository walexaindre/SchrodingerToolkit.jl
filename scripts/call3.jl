include("examples/3D/PDE_ex3D_1.jl")

using GLMakie
using CUDA
CUDA.allowscalar(false)
CUDA.device!(1)

backend = GPUBackend{Int64,Float64}

Params = SolverParameters(backend, 2, (:ord4, :ord4, :ord4), :tord2_1_1)
Grid = PeriodicGrid(backend,PDE,0.01,(0.16,0.16,0.16))
Method, Memory, Stats = M1(PDE,Params,Grid)
Problem = SchrodingerProblem(Method, Memory, Stats, PDE, Params)


colors = to_colormap(:balance)
n = length(colors)
alpha = [ones(n÷3);zeros(n-2*(n÷3));ones(n÷3)]
res  = -8..8
resv = collect(res)
rs = map(x->(x.r,x.g,x.b),colors)
cmap_alpha = map(x->RGBAf(x[2][1],x[2][2],x[2][3],alpha[x[1]]),enumerate(rs))

current_mem = current_state(Memory)|> Array
v1 = reshape(abs2.(current_mem[:,1]),(100,100,100))

volume!(ax,-8..8,-8..8,-8..8,obs,algorithm=:mip)

obs = Observable(v1)
fig = volume(res,res,res,obs,algorithm=:mip,colormap=cmap_alpha)
hidespines!(fig)
hidedecorations!(fig.axis)
hidespines!(ax)
hidedecorations!(ax)

fig

for i in 1:30
    for j in 1:10
        step!(Problem)
    end
    

end

GLMakie.record(f,"3D.gif",1:400;framerate=20) do t
    obs[] = reshape(abs2.(current_state!(Memory)[:,1]|>Array),(100,100,100))
    step!(Problem)
end

f = Figure()

ax = Axis3(f[1,1])

f

save("./3dplt.png", f)