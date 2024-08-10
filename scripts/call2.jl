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





















@recipe(_ExecutionTime,Stats) do scene
    Attributes(
        color = theme(scene,:color),
        colormap = theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end





function Makie.plot!(Sys::_ExecutionTime)
    Stats = Sys[1]
    τ = Stats[].τ
    log_freq = Stats[].log_frequency

    points = Observable(Float64[])    
    timerange = Observable(range(τ,step=τ,length=4))
    avg = Observable(0.0)

    function update_plot(Stats)
        empty!(points[])
        maxindex = Stats[].store_index-1
        timerange[] = range(Stats[].τ,step=Stats[].τ*log_freq,length=maxindex)
        append!(points[],Stats[].step_time[1:maxindex])
        avg[] = sum(points[])/maxindex
    end

    update_plot(Stats)
    
    Makie.scatterlines!(Sys,timerange,points,  color = RGBf(0.0039,0.239216,0.5333) , colormap = Sys.colormap, inspectable = Sys.inspectable, visible = Sys.visible,marker=:hexagon,strokewidth=1,linewidth=2)
    Makie.ablines!(Sys,avg[],0,color= RGBf(0.698,0.168,0.0745),linestyle=:dash,linewidth=2,fontsize=Sys.fontsize)
    Sys

end

function executiontime(Stats,ounit)
    p=_executiontime(Stats,ounit,fontsize=50)

    p.axis.xlabel = rich("t",subscript("n"))
    p.axis.ygridvisible = false
    p.axis.xgridvisible = false
    p.axis.yminorgridvisible = false
    p.axis.yticks = WilkinsonTicks(6,k_min=5)
    p.axis.xticks = WilkinsonTicks(6,k_min=5)
    p.axis.xlabelsize = p.axis.ylabelsize = 24
    p.axis.xticklabelsize=p.axis.yticklabelsize=24
    p.axis.ytickformat = xs -> [string(transform(v,ounit)) for v in xs]
    p    
end


fig = Figure(size = (800, 600), fontsize = 40,dpi=300)

Ax = Axis(fig[1,1],xlabel = "Time (ms)", ylabel = "Step time (ms)")
Ax.xlabel = rich("t",subscript("n"),"(ms)")
Ax.ygridvisible = false
Ax.xgridvisible = false
Ax.yminorgridvisible = false
Ax.yticks = WilkinsonTicks(6,k_min=5)
Ax.xticks = WilkinsonTicks(6,k_min=5)
Ax.xlabelsize = 24
Ax.xticklabelsize=24

scatterlines!(Ax,0..3,1000*rmx,  color = RGBf(0.0039,0.239216,0.5333) , colormap = :viridis, inspectable = true, visible = true,marker=:hexagon,strokewidth=1,linewidth=2)
fig

save("./nonsolvingtime.png", fig)