using GLMakie

fig = diagram_startup()
ax2 = Axis(fig[1,1])
diagram_axis_conf(ax2)
samplefig = lines!(ax2, 0..10, sin, color = :red)
function diagram_startup(; fig_height::Int = 1000, fig_ratio = 4 / 3)
    fig_width = ceil(Int, fig_height * fig_ratio)
    fig_size = (fig_width, fig_height)
    fig = Figure(; size = fig_size)
    return fig
end


function diagram_axis_conf(axis)

    is3d = hasproperty(axis,:zticks)

    ndims = is3d ? 3 : 2
    prefix = is3d ? [:x,:y,:z] : [:x,:y]   

    property = Dict{Symbol,Any}()

    property[:labelsize] = 24
    property[:ticklabelsize] = 18
    property[:ticks] = WilkinsonTicks(7; k_min = 6)
  
    for dim in prefix
        for (prop,value) in property
            setproperty!(axis,Symbol(dim,prop),value)
        end
    end
end

macro diagbase(figure)
    return quote
        println($figure)
    end
end
macro define_var(varname, value)
    quote
        $(esc(varname)) = $(esc(value))
    end
end