#Unstable... Makie.jl isn't at release 1.0 yet.

function plot_systemnd(Memory, Grid, idx; kwargs...)
    ndimsv = ndims(Grid)
    x = get_range(Grid, 1)
    y = get_range(Grid, 2)

    xel = length(x)
    yel = length(y)

    fig_height = haskey(kwargs, :fig_height) ? kwargs[:fig_height] : 1000
    fig_ratio = haskey(kwargs, :fig_ratio) ? kwargs[:fig_ratio] : 4 / 3
    x_label = haskey(kwargs, :x_label) ? kwargs[:x_label] : "x  ($xel points)"
    y_label = haskey(kwargs, :y_label) ? kwargs[:y_label] : "y  ($yel points)"
    z_label = haskey(kwargs, :z_label) ? kwargs[:z_label] :
              L"\left| \psi \right|^2_{%$idx}" #Interpolation for LaTeX %$ variable
    colormap = haskey(kwargs, :colormap) ? kwargs[:colormap] : :default

    fig = diagram_startup(; fig_height, fig_ratio)

    axis = Axis3(fig[1, 1])
    diagram_axis_conf(axis)
    diagram_label_conf(axis; xlabel = x_label, ylabel = y_label, zlabel = z_label)
    diagram_hide_conf(axis)

    current_state = current_state!(Memory)

    if ndimsv == 2
        points = Array(abs2.(view(current_state, :, idx)))
        z = reshape(points, size(Grid))
        surf_scene = surface!(axis, x, y, z; colormap = cgrad(colormap; rev = true),
                              kwargs...)
        Colorbar(fig[1, 2], surf_scene; ticklabelsize = 18)

    else
        Error("3D plotting not supported yet.")
    end
    fig, axis
end

function plot_execution_time(Stats, Grid; kwargs...)
    linewidth = haskey(kwargs, :linewidth) ? kwargs[:linewidth] : 2.5
    linestyle = haskey(kwargs, :linestyle) ? kwargs[:linestyle] : :dash
    color = haskey(kwargs, :color) ? kwargs[:color] : RGBf(0.0039, 0.239216, 0.5333)
    colormap = haskey(kwargs, :colormap) ? kwargs[:colormap] : :default

    fig = diagram_startup()
    axis = Axis(fig[1, 1])
    diagram_axis_conf(axis)
    diagram_label_conf(axis; xlabel = rich("t", subscript("n")),
                       ylabel = "Step time (ms)")
    diagram_hide_conf(axis)

    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    steps = step_time(Stats) * 1000.0
    avg = mean(steps) * 1000.0

    scatterlines!(axis, time_range, steps; color = color, linewidth = linewidth,
                  linestyle = linestyle, colormap = colormap, kwargs...)
    ablines!(axis, avg, 0; linestyle = linestyle, linewidth = linewidth, kwargs...)

    fig, axis
end

function plot_solver_time(Stats,Grid; kwargs...)
    linewidth = haskey(kwargs, :linewidth) ? kwargs[:linewidth] : 2.5
    linestyle = haskey(kwargs, :linestyle) ? kwargs[:linestyle] : :dash
    color = haskey(kwargs, :color) ? kwargs[:color] : RGBf(0.0039, 0.239216, 0.5333)
    colormap = haskey(kwargs, :colormap) ? kwargs[:colormap] : :default

    fig = diagram_startup()
    axis = Axis(fig[1, 1])
    diagram_axis_conf(axis)
    diagram_label_conf(axis; xlabel = rich("t", subscript("n")),
                       ylabel = "Solver time (ms)")
    diagram_hide_conf(axis)

    sol_time = solver_time(Stats) * 1000.0
    avg = mean(sol_time) * 1000.0

    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    scatterlines!(axis, time_range, sol_time; color = color, linewidth = linewidth,
                  linestyle = linestyle, colormap = colormap, kwargs...)

    ablines!(axis, avg, 0; linestyle = linestyle, linewidth = linewidth, kwargs...)

    fig, axis
end

function plot_mass_per_component(Stats, index)
end

export plot_systemnd, plot_execution_time, plot_mass_per_component, plot_solver_time
