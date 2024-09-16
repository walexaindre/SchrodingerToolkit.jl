# [TODO]: Makie isn't at release 1.0 yet. Some code may need to be updated in the future.


"""
    diagram_startup(; fig_height::Int = 1000, fig_ratio = 4 / 3)

Create a figure with the specified height and ratio.

# Keyword Arguments
- `fig_height::Int`: The height of the figure.
- `fig_ratio`: The ratio of the figure.

# Returns
- `fig`: The figure.
"""
function diagram_startup(; fig_height::Int = 1000, fig_ratio = 4 / 3)
    fig_width = ceil(Int, fig_height * fig_ratio)
    fig_size = (fig_width, fig_height)
    fig = Figure(; size = fig_size)
    return fig
end

"""
    diagram_hide_conf(axis)

Hide the grid and minor grid lines.

# Arguments
- `axis`: The axis.
"""
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

"""
    diagram_hide_env(axis)

Hide the spines and decorations.

# Arguments
- `axis`: The axis.
"""
function diagram_hide_env(axis)
    hidespines!(axis)
    hidedecorations!(axis)
end

"""
    diagram_axis_conf(axis; labelsize = 48, ticklabelsize = 42)

Configure the axis.

# Arguments
- `axis`: The axis.

# Keyword Arguments
- `labelsize::Int`: The size of the label.
- `ticklabelsize::Int`: The size of the tick label.
"""
function diagram_axis_conf(axis; labelsize = 48, ticklabelsize = 42)
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

"""
    diagram_label_conf(axis; xlabel = "x-label", ylabel = "y-label", zlabel = "z-label")

Configure the labels.

# Arguments
- `axis`: The axis.

# Keyword Arguments
- `xlabel::String`: The x-label.
- `ylabel::String`: The y-label.
- `zlabel::String`: The z-label.
"""
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


"""
    diagram_base_2d(xval, yval; xlabel = "x - label", ylabel = "y - label", kwargs...)

Create a 2D diagram with average line.

# Arguments
- `xval`: The x-values.
- `yval`: The y-values.

# Keyword Arguments
- `xlabel::String`: The x-label.
- `ylabel::String`: The y-label.
- `kwargs...`: Additional keyword arguments passed to `Makie`.

# Returns
- `fig`: The figure.
- `axis`: The axis.
"""
function diagram_base_2d(xval, yval, avg; xlabel = "x - label", ylabel = "y - label", kwargs...)
    linewidth = haskey(kwargs, :linewidth) ? kwargs[:linewidth] : 2.5
    linestyle = haskey(kwargs, :linestyle) ? kwargs[:linestyle] : :dash
    color = haskey(kwargs, :color) ? kwargs[:color] : RGBf(0.0039, 0.239216, 0.5333)
    colormap = haskey(kwargs, :colormap) ? kwargs[:colormap] : :default

    fig = diagram_startup()
    axis = Axis(fig[1, 1])
    diagram_axis_conf(axis)
    diagram_label_conf(axis; xlabel = xlabel, ylabel = ylabel)
    diagram_hide_conf(axis)

    scatterlines!(axis, xval, yval; color = color, linewidth = linewidth,
                  linestyle = linestyle, colormap = colormap, kwargs...)

    ablines!(axis, avg, 0; color = :red,linestyle = linestyle, linewidth = linewidth, kwargs...)

    fig, axis
end

"""
    diagram_2d(xval, yval; xlabel = "x - label", ylabel = "y - label", xscale = identity, yscale = identity, kwargs...)

Create a 2D diagram.

# Arguments
- `xval`: The x-values.
- `yval`: The y-values.

# Keyword Arguments
- `xlabel::String`: The x-label.
- `ylabel::String`: The y-label.

# Returns
- `fig`: The figure.
- `axis`: The axis.
"""
function diagram_2d(xval,yval; xlabel= "x - label",ylabel= "y - label",xscale = identity, yscale = identity, kwargs...)
    linewidth = haskey(kwargs, :linewidth) ? kwargs[:linewidth] : 2.5
    linestyle = haskey(kwargs, :linestyle) ? kwargs[:linestyle] : :dash
    color = haskey(kwargs, :color) ? kwargs[:color] : RGBf(0.0039, 0.239216, 0.5333)
    colormap = haskey(kwargs, :colormap) ? kwargs[:colormap] : :default

    fig = diagram_startup()
    axis = Axis(fig[1, 1])
    axis.xscale = xscale
    axis.yscale = yscale

    diagram_axis_conf(axis)
    diagram_label_conf(axis; xlabel = xlabel, ylabel = ylabel)
    diagram_hide_conf(axis)
    
    scatterlines!(axis, xval, yval; color = color, linewidth = linewidth,
                  linestyle = linestyle, colormap = colormap, kwargs...)              

    fig,axis
    
end