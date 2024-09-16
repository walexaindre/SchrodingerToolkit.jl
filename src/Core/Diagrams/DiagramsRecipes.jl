#Unstable... Makie.jl isn't at release 1.0 yet.

"""
    plot_systemnd(Memory, Grid, idx; kwargs...)

Plot the system at current time step.

# Arguments
- `Memory`: The memory struct.
- `Grid`: The grid.
- `idx`: The index of the component.
- `kwargs...`: Additional keyword arguments passed to `Makie`.
"""
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
        Colorbar(fig[1, 2], surf_scene; ticklabelsize = 42)

    else
        Error("3D plotting not supported yet.")
    end
    fig, axis
end

"""
    plot_execution_time(Stats::RTStats, Grid; kwargs...)

Plot the execution time at each time step.

# Arguments
- `Stats::RTStats`: The runtime statistics.
- `Grid`: The grid.
- `kwargs...`: Additional keyword arguments passed to `diagram_base_2d`.
"""
function plot_execution_time(Stats::RTStats, Grid; kwargs...) where {RTStats <: AbstractRuntimeStats}
    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    steps = step_time(Stats) * 1000.0
    avg = mean(steps) * 1000.0

    diagram_base_2d(time_range, steps, avg; xlabel = rich("t", subscript("n")),
                    ylabel = "Step time (ms)",
                    kwargs...)
end

"""
    plot_solver_time(Stats::RTStats, Grid; kwargs...)

Plot the solver time at each time step.

# Arguments
- `Stats::RTStats`: The runtime statistics.
- `Grid`: The grid.
- `kwargs...`: Additional keyword arguments passed to `diagram_base_2d`.
"""
function plot_solver_time(Stats::RTStats, Grid; kwargs...) where {RTStats <: AbstractRuntimeStats}
    sol_time = solver_time(Stats) * 1000.0
    avg = mean(sol_time) * 1000.0

    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    diagram_base_2d(time_range, sol_time, avg; xlabel = rich("t", subscript("n")),
                    ylabel = "Solver time (ms)", kwargs...)
end

function plot_preprocessing_time(Stats::RTStats, Grid; kwargs...) where {RTStats <: AbstractRuntimeStats}
    diffv = calculate_diff_step_and_solver_time(Stats)
    avg = mean(diffv) * 1000.0

    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    diagram_base_2d(time_range, diffv, avg; xlabel = rich("t", subscript("n")),
                    ylabel = "Preprocessing time (ms)", kwargs...)
end

"""
    plot_absolute_error_mass_per_component(Stats::RTStats, Grid, index; kwargs...)

Plot the absolute error in the mass at each component at each time step.

# Arguments
- `Stats::RTStats`: The runtime statistics.
- `Grid`: The grid.
- `index`: The index of the component.
- `kwargs...`: Additional keyword arguments passed to `diagram_2d`.
"""
function plot_absolute_error_mass_per_component(Stats::RTStats, Grid, index; kwargs...) where {RTStats <: AbstractRuntimeStats}
    sys_mass_diff = calculate_diff_system_mass(Stats, index)
    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    diagram_2d(time_range, sys_mass_diff; xlabel = rich("t", subscript("n")),
                    ylabel = "Absolute Error",yscale=Makie.pseudolog10, kwargs...)
end


"""
    plot_mass_per_component(Stats::RTStats, Grid, index; kwargs...)

Plot the mass at each component at each time step.

# Arguments
- `Stats::RTStats`: The runtime statistics.
- `Grid`: The grid.
- `index`: The index of the component.
- `kwargs...`: Additional keyword arguments passed to `diagram_2d`.
"""
function plot_mass_per_component(Stats::RTStats, Grid, index; kwargs...) where {RTStats <: AbstractRuntimeStats}
    sys_mass_diff = system_mass(Stats, index)
    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    diagram_2d(time_range, sys_mass_diff; xlabel = rich("t", subscript("n")),
                    ylabel = "Mass at component $index", kwargs...)
end

"""
    plot_absolute_error_total_mass(Stats::RTStats, Grid)

Plot the absolute error in the total mass at each time step.

# Arguments
- `Stats::RTStats`: The runtime statistics.
- `Grid`: The grid.
"""
function plot_absolute_error_total_mass(Stats::RTStats, Grid) where {RTStats <: AbstractRuntimeStats}
    sys_mass_diff = calculate_diff_system_total_mass(Stats)
    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    diagram_2d(time_range, sys_mass_diff; xlabel = rich("t", subscript("n")),
                    ylabel = "Absolute error",yscale = Makie.pseudolog10)
end

"""
    plot_total_mass(Stats::RTStats, Grid)

Plot the total mass at each time step.

# Arguments
- `Stats::RTStats`: The runtime statistics.
- `Grid`: The grid.
"""
function plot_total_mass(Stats::RTStats, Grid) where {RTStats <: AbstractRuntimeStats}
    sys_mass_diff = system_total_mass(Stats)
    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    diagram_2d(time_range, sys_mass_diff; xlabel = rich("t", subscript("n")),
                    ylabel = "Total mass")
end


"""
    plot_solver_iterations(Stats::RTStats, Grid; kwargs...)

Plot the number of iterations taken by the solver at each time step.

# Arguments
- `Stats::RTStats`: The runtime statistics.
- `Grid`: The grid.
- `kwargs...`: Additional keyword arguments passed to `diagram_base_2d`.
"""
function plot_solver_iterations(Stats::RTStats , Grid; kwargs...) where {RTStats <: AbstractRuntimeStats}
    iter = solver_iterations(Stats)
    avg = mean(iter)

    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    diagram_base_2d(time_range, iter, avg; xlabel = rich("t", subscript("n")),
                    ylabel = "Solver iterations", kwargs...)
end

"""
    plot_absolute_error_system_energy(Stats::RTStats, Grid; kwargs...)

Plot the absolute error in the system energy at each time step.

# Arguments
- `Stats::RTStats`: The runtime statistics.
- `Grid`: The grid.
- `kwargs...`: Additional keyword arguments passed to `diagram_2d`.
"""
function plot_absolute_error_system_energy(Stats::RTStats, Grid; kwargs...) where {RTStats <: AbstractRuntimeStats}
    sys_energy = calculate_diff_system_energy(Stats)
    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)
    diagram_2d(time_range, sys_energy; xlabel = rich("t", subscript("n"),yscale = Makie.pseudolog10),
                    ylabel = "System Energy", kwargs...)
    
end

"""
    plot_component_update_steps(Stats::RTStats, Grid; kwargs...)

Plot the number of steps taken by the component update function at each time step.

# Arguments
- `Stats::RTStats`: The runtime statistics.
- `Grid`: The grid.
- `kwargs...`: Additional keyword arguments passed to `diagram_base_2d`.
"""
function plot_component_update_steps(Stats::RTStats, Grid; kwargs...) where {RTStats <: AbstractRuntimeStats}
    steps = component_update_steps(Stats)
    avg = mean(steps)
    log_freq = Stats.log_frequency
    τ = Grid.τ * log_freq

    len = length(Stats)

    time_range = range(τ; step = τ, length = len)

    diagram_base_2d(time_range, steps, avg; xlabel = rich("t", subscript("n")),
                    ylabel = "Component update steps", kwargs...)
    
end

export plot_systemnd, plot_execution_time, plot_mass_per_component, plot_solver_time,
       plot_preprocessing_time, plot_absolute_error_mass_per_component,
       plot_absolute_error_total_mass, plot_total_mass, plot_solver_iterations, plot_absolute_error_system_energy, plot_component_update_steps
