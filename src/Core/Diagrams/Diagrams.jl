# [TODO]: Makie isn't at release 1.0 yet. Some code may need to be updated in the future.

function diagram_startup(; fig_height::Int = 1000, fig_ratio = 4 / 3)
    fig_width = ceil(Int, fig_height * fig_ratio)
    fig_size = (fig_width, fig_height)
    fig = Figure(; size = fig_size)
    return fig
end

function diagram_hide_conf(axis)
    is3d = hasproperty(axis, :zticks)

    ndims = is3d ? 3 : 2
    prefix = is3d ? [:x, :y, :z] : [:x, :y]

    property = Dict{Symbol,Any}()

    property[:gridvisible] = false
    if ndims == 2
        property[:minorgridvisible] = false
    end

    for dim in prefix
        for (prop, value) in property
            setproperty!(axis, Symbol(dim, prop), value)
        end
    end
end

function diagram_hide_env(axis)
    hidespines!(axis)
    hidedecorations!(axis)
end

function diagram_axis_conf(axis; labelsize = 24, ticklabelsize = 18)
    is3d = hasproperty(axis, :zticks)

    ndims = is3d ? 3 : 2
    prefix = is3d ? [:x, :y, :z] : [:x, :y]

    property = Dict{Symbol,Any}()

    property[:labelsize] = labelsize
    property[:ticklabelsize] = ticklabelsize
    property[:ticks] = WilkinsonTicks(7; k_min = 6)

    for dim in prefix
        for (prop, value) in property
            setproperty!(axis, Symbol(dim, prop), value)
        end
    end
end

function diagram_label_conf(axis; xlabel = "x-label", ylabel = "y-label",
                            zlabel = "z-label")
    is3d = hasproperty(axis, :zticks)

    ndims = is3d ? 3 : 2

    property = Dict{Symbol,Any}()

    if ndims >= 2
        property[:xlabel] = xlabel
        property[:ylabel] = ylabel
    end

    if ndims == 3
        property[:zlabel] = zlabel
    end

    for (prop, value) in property
        setproperty!(axis, prop, value)
    end
end