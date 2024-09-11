Base.length(R::RuntimeStats{IntType,FloatType}) where {IntType,FloatType} = R.store_index -
                                                                            1
Base.size(C::ComponentMass{FloatType,ArrayType}) where {FloatType,ArrayType} = size(C.mass)
Base.length(C::ComponentMass{FloatType,ArrayType}) where {FloatType,ArrayType} = length(C.mass)

@inline current_iteration(R::RuntimeStats) = R.current_iteration

@inline start_energy(R::RuntimeStats) = R.system_energy[end]
@inline system_energy(R::RuntimeStats) = R.system_energy
@inline system_energy(R::RuntimeStats, index) = R.system_energy[index]
@inline system_mass(R::RuntimeStats, component) = R.system_mass[component]
@inline system_mass(R::RuntimeStats, component, index) = R.system_mass[component][index]
@inline system_total_mass(R::RuntimeStats) = R.system_total_mass
@inline start_mass(R::RuntimeStats, component) = R.system_mass[component][end]
@inline start_mass(R::RuntimeStats) = map(x -> x[end], R.system_mass)
@inline start_total_mass(R::RuntimeStats) = R.system_total_mass[end]

@inline component_update_steps(R::RuntimeStats) = R.component_update_steps[1:length(R)]
@inline component_update_calls(R::RuntimeStats) = R.component_update_call_count

"""
    isregistering_stats(R::RuntimeStats)

Check if the stats are being registered.
"""
@inline isregistering_stats(R::RuntimeStats) = R.config[5]

"""
    islog(R::RuntimeStats)

Check if the stats are being logged. Note that this is different from `isregistering_stats`. This is used to check if the stats are being logged at the current time step.
 
"""
@inline islog(R::RuntimeStats) = R.log_data

"""
    islocked(R::RuntimeStats)

Check if the stats are locked. This is used to check if the stats are locked and cannot be modified in any way.

## Note

Any attempt to use or modify the `Stats` structure if it is locked will throw an error.
"""
@inline islocked(R::RuntimeStats) = R.config[4]

@inline islog_system_total_mass(R::RuntimeStats) = R.config[2]
@inline islog_solver_info(R::RuntimeStats) = R.config[1]
@inline islog_component_update_steps(R::RuntimeStats) = R.config[3]

@inline function Base.getindex(C::ComponentMass{FloatType,
                                                 ArrayType},
                               index) where {FloatType,ArrayType}
    C.mass[index]
end

@inline function Base.setindex!(C::ComponentMass{FloatType,ArrayType},
                                value::FloatType,
                                index) where {FloatType,ArrayType}
    C.mass[index] = value
end

"""
    initialize_stats(FloatVectorType, IntVectorType, ncomponents, log_freq, time_steps;
                     log_stats = true, log_solver_info = true,
                     log_component_update_steps = false,
                     log_system_total_mass = false, islocked = false)

Initialize the stats structure.

## Arguments

- `FloatVectorType`: The type of the floating point vector.
- `IntVectorType`: The type of the integer vector.
- `ncomponents`: The number of components.
- `log_freq`: The frequency to log the stats.
- `time_steps`: The total number of time steps.

## Optional arguments

- `log_stats`: Log the stats.
- `log_solver_info`: Log the solver information.
- `log_component_update_steps`: Log the component update steps.
- `log_system_total_mass`: Log the system total mass.
- `islocked`: Lock the stats.

## Returns

- `RuntimeStats`: The stats structure.

"""
function initialize_stats(::Type{FloatVectorType}, ::Type{IntVectorType},
                          ncomponents::IntType, log_freq::IntType,
                          time_steps::IntType; log_stats::Bool = true,
                          log_solver_info::Bool = true,
                          log_component_update_steps::Bool = false,
                          log_system_total_mass::Bool = false,
                          islocked::Bool = false) where {IntType,
                                                         IntVectorType,
                                                         FloatVectorType}
    config = (log_solver_info, log_system_total_mass, log_component_update_steps,
              islocked, log_stats, false, false, false)

    seq_size = div(time_steps, log_freq) + 1
    seq_solver = log_solver_info ? seq_size : 0
    seq_total_mass = log_system_total_mass ? seq_size : 0
    seq_component_update_steps = log_component_update_steps ? seq_size : 0

    sys_energy = vundef(FloatVectorType, seq_size)
    sys_mass = ntuple(Returns(ComponentMass(similar(sys_energy))), ncomponents)
    sys_total_mass = vundef(FloatVectorType, seq_total_mass)
    sys_time = vundef(FloatVectorType, seq_size)
    solver_time = vzeros(FloatVectorType, seq_solver)
    solver_iterations = vzeros(IntVectorType, seq_solver)
    component_update_steps = vzeros(IntVectorType, seq_component_update_steps)

    RuntimeStats(sys_energy, sys_mass, sys_total_mass, sys_time, solver_time,
                 solver_iterations, component_update_steps,
                 log_freq, config, false, zero(IntType), zero(IntType), one(IntType))
end

function initialize_stats(::Type{FloatVectorType}, ::Type{IntVectorType}, PDE::SPDE,
                          conf::SolverConfig,
                          grid::PGrid) where {FloatVectorType,IntVectorType,
                                              SPDE<:SchrodingerPDE,
                                              PGrid<:PeriodicGrid}
    log_freq = stats_logfreq(conf)
    log_component_update_steps = stats_log_component_update_steps(conf)
    log_system_total_mass = stats_log_system_total_mass(conf)
    log_solver_info = stats_log_solver_info(conf)
    islocked = stats_locked(conf)
    log_stats = stats_log_stats(conf)

    ncomp = ncomponents(PDE)
    tsteps = estimate_timesteps(PDE, grid)
    initialize_stats(FloatVectorType, IntVectorType, ncomp, log_freq, tsteps;
                     log_stats = log_stats,
                     log_solver_info = log_solver_info,
                     log_component_update_steps = log_component_update_steps,
                     log_system_total_mass = log_system_total_mass,
                     islocked = islocked)
end

"""
    startup_stats!(stats::Stats, start_mass, start_energy)

Initialize the stats structure with the mass and energy at time step `0`.
"""
function startup_stats!(stats::Stats, start_mass,
                        start_energy) where {Stats<:RuntimeStats}
    last_index = lastindex(stats.step_time)
    update_mass!(stats, start_mass, last_index)
    update_system_energy!(stats, start_energy, last_index)
    advance_iteration!(stats)
end

"""
    startup_stats!(stats::Stats, start_mass, start_energy, total_mass)

Initialize the stats structure with the mass, energy and total mass at time step `0`.
"""
function startup_stats!(stats::Stats, start_mass, start_energy,
                        total_mass) where {Stats<:RuntimeStats}
    last_index = lastindex(stats.step_time)
    update_mass!(stats, start_mass, last_index)
    update_system_energy!(stats, start_energy, last_index)
    update_system_total_mass!(stats, total_mass, last_index)
    advance_iteration!(stats)
end

"""
    advance_iteration!(R::RuntimeStats{IntType,FloatType}) where {IntType<:Integer,
                                                                 FloatType<:AbstractFloat}

Advance the iteration in the stats structure.
"""
@inline function advance_iteration!(R::RuntimeStats{IntType,FloatType}) where {IntType<:Integer,
                                                                               FloatType<:AbstractFloat}
    if !islocked(R) && isregistering_stats(R)
        R.current_iteration += 1
        if mod(R.current_iteration, R.log_frequency) == 0
            R.log_data = true
        end
        return true
    elseif islocked(R)
        error("Stats are locked")
    end
    return false
end
"""
    update_mass!(stats, mass, idx)

Update the mass per component in the stats structure.

## Arguments
- `stats::Stats`: The stats structure.
- `mass::Vector{FloatType}`: The mass per component.
- `idx::IntType`: The index to update.
"""
function update_mass!(stats::Stats, mass,
                       idx::IntType) where {IntType<:Integer,
                                            FloatType<:AbstractFloat,
                                            Stats<:RuntimeStats{IntType,FloatType}}
    mass = Vector(mass)

    if !islocked(stats) && isregistering_stats(stats)
        for comp in 1:length(stats.system_mass)
            stats.system_mass[comp][idx] = mass[comp]
        end
        return true
    elseif islocked(stats)
        error("Stats are locked")
    end
    false
end

"""
    update_system_energy!(stats, energy, idx)

Update the system energy in the stats structure.

## Arguments
- `stats::Stats`: The stats structure.
- `energy::FloatType`: The energy.
- `idx::IntType`: The index to update.
"""
function update_system_energy!(stats::Stats, energy::FloatType,
                               idx::IntType) where {IntType<:Integer,
                                                    FloatType<:AbstractFloat,
                                                    Stats<:RuntimeStats{IntType,
                                                                        FloatType}}
    if !islocked(stats) && isregistering_stats(stats)
        stats.system_energy[idx] = energy
        return true
    elseif islocked(stats)
        error("Stats are locked")
    end
    false
end

"""
    update_system_total_mass!(stats, total_mass, idx)

Update the system total mass in the stats structure.

## Arguments
- `stats::Stats`: The stats structure.
- `total_mass::FloatType`: The total mass.
- `idx::IntType`: The index to update.
"""
function update_system_total_mass!(stats::Stats, total_mass::FloatType,
                                    idx::IntType) where {IntType<:Integer,
                                                         FloatType<:AbstractFloat,
                                                         Stats<:RuntimeStats{IntType,
                                                                             FloatType}}
    if !islocked(stats) && islog_system_total_mass(stats) && isregistering_stats(stats)
        stats.system_total_mass[idx] = total_mass
        return true
    elseif islocked(stats)
        error("Stats are locked")
    end
    return false
end

"""
    update_solver_info!(stats, time, iterations)

Update the linear solver information in the stats structure.

## Arguments
- `stats::Stats`: The stats structure.
- `time::FloatType`: The time spent in the solver.
- `iterations::IntType`: The number of iterations performed in the solver.
"""
function update_solver_info!(stats::Stats, time,
                             iterations) where {Stats<:RuntimeStats}
    if !islocked(stats) && islog_solver_info(stats) && isregistering_stats(stats)
        if islog(stats)
            idx = stats.store_index
            stats.solver_time[idx] += time
            stats.solver_iterations[idx] += iterations
        end
    elseif islocked(stats)
        error("Stats are locked")
    end
end

"""
    update_component_update_steps!(stats, steps)

Update the component update steps in the stats structure.

## Arguments
- `stats::Stats`: The stats structure.
- `steps::IntType`: The number of steps performed in the component update.
"""
function update_component_update_steps!(stats::Stats, steps) where {Stats<:RuntimeStats}
    if !islocked(stats) && islog_component_update_steps(stats) &&
       isregistering_stats(stats)
        if islog(stats)
            idx = stats.store_index
            stats.component_update_steps[idx] += steps
            stats.component_update_call_count += 1
            return true
        end
        return false
    elseif islocked(stats)
        error("Stats are locked")
    end
    return false
end

"""
    calculate_diff_system_energy(stats)

Calculate the absolute value of the difference between the system energy and the energy at time step `0`.

## Arguments
- `stats::Stats`: The stats structure.
"""
function calculate_diff_system_energy(stats::Stats) where {Stats<:RuntimeStats}
    startup_energy = start_energy(stats)

    abs.(system_energy(stats)[1:length(stats)] .- startup_energy)
end

"""
    calculate_diff_system_mass(stats)

Calculate the absolute value of the difference between the system mass and the mass at time step `0`.

## Arguments
- `stats::Stats`: The stats structure.
"""
function calculate_diff_system_mass(stats::Stats,
                                     index) where {Stats<:RuntimeStats}
    startup_mass = start_mass(stats, index)

    return abs.(system_mass(stats, index)[1:length(stats)] .- startup_mass)
end

"""
    calculate_diff_system_total_mass(stats)

Calculate the absolute value of the difference between the system total mass and the total mass at time step `0`.

## Arguments
- `stats::Stats`: The stats structure.
"""
function calculate_diff_system_total_mass(stats::Stats) where {Stats<:RuntimeStats}
    startup_mass = start_total_mass(stats)

    abs.(system_total_mass(stats)[1:length(stats)] .- startup_mass)
end

"""
    step_time(stats)

Get the time spent in the current step.

## Arguments
- `stats::Stats`: The stats structure.
"""
step_time(stats::Stats) where {Stats<:RuntimeStats} = stats.step_time[1:length(stats)]

"""
    solver_time(stats)

Get the time spent in the solver.

## Arguments
- `stats::Stats`: The stats structure.
"""
solver_time(stats::Stats) where {Stats<:RuntimeStats} = stats.solver_time[1:length(stats)]

"""
    calculate_diff_step_and_solver_time(stats)

Calculate the absolute value of the difference between the time spent in the current step and the time spent in the solver.

## Arguments
- `stats::Stats`: The stats structure.
"""
function calculate_diff_step_and_solver_time(stats::Stats) where {Stats<:RuntimeStats}
    step_time = step_time(stats)
    solver_time = solver_time(stats)
    return abs.(step_time .- solver_time)
end

"""
    update_stats!(stats, step_time, mass_per_component, sys_energy)

Update the stats structure with the current time, mass per component and system energy.

## Arguments
- `stats::Stats`: The stats structure.
- `step_time::FloatType`: Amount of time spent in the current step.
- `mass_per_component`: The mass per component.
- `sys_energy::FloatType`: The system energy.
"""
function update_stats!(stats::Stats, step_time, mass_per_component,
                       sys_energy) where {Stats<:RuntimeStats}
    if !islocked(stats) && isregistering_stats(stats)
        if islog(stats)
            idx = stats.store_index
            stats.step_time[idx] = step_time
            update_mass!(stats, mass_per_component, idx)
            update_system_energy!(stats, sys_energy, idx)
            stats.store_index += 1
            stats.log_data = false
        end
        advance_iteration!(stats)
    elseif islocked(stats)
        error("Stats are locked")
    end
end

"""
    update_stats!(stats, step_time, mass_per_component, sys_energy, total_mass)

Update the stats structure with the current time, mass per component, system energy and total mass.

## Arguments
- `stats::Stats`: The stats structure.
- `step_time::FloatType`: Amount of time spent in the current step.
- `mass_per_component`: The mass per component.
- `sys_energy::FloatType`: The system energy.
- `total_mass::FloatType`: The total mass.
"""
function update_stats!(stats::Stats, step_time, mass_per_component, sys_energy,
                       total_mass) where {Stats<:RuntimeStats}
    if !islocked(stats) && isregistering_stats(stats)
        if islog(stats)
            idx = stats.store_index
            stats.step_time[idx] = step_time
            update_mass!(stats, mass_per_component, idx)
            update_system_energy!(stats, sys_energy, idx)
            update_system_total_mass!(stats, total_mass, idx)
            stats.store_index += 1
            stats.log_data = false
        end
        advance_iteration!(stats)
    elseif islocked(stats)
        error("Stats are locked")
    end
end

"""
    update_stats!(stats, memory, grid, PDE, work_timer)

Update the stats structure with the current time, mass per component, system energy and total mass.

## Arguments
- `stats::Stats`: The stats structure.
- `memory`: The memory structure.
- `grid`: The grid structure.
- `PDE`: The PDE structure.
- `work_timer`: The timer for the work done in the current step.

## Side effects
- Advances the iteration in the stats structure no matter if the stats are being logged or not. 
This is to keep the stats structure consistent with the simulation. 
For that reason is important to call this function at the end of the simulation step.
"""
function update_stats!(stats::Stats, memory::MemType, grid, PDE,
                       work_timer) where {Stats<:RuntimeStats,MemType<:AbstractMemory}
    if isregistering_stats(stats)
        if islog(stats)
            mass_per_component = system_mass(memory, grid)
            energy = system_energy(memory, PDE, grid)
            if islog_system_total_mass(stats)
                total_mass = system_total_mass(memory, PDE, grid,
                                                 mass_per_component)
                update_stats!(stats, work_timer, mass_per_component, energy, total_mass)
            end
            update_stats!(stats, work_timer, mass_per_component, energy)
        else
            advance_iteration!(stats)
        end
    end
end

function serialize(array, len)
    output = vundef(typeof(array), len + 1)
    output[2:(len + 1)] .= array[1:len]
    output[1] = array[end]
    output
end

function deserialize(::Type{VectorType}, array) where {VectorType<:AbstractVector}
    sz = length(array)
    output = vundef(VectorType, sz)
    output[1:(sz - 1)] .= array[2:end]
    output[end] = array[1]
    output
end

function serialize(stats::Stats, path) where {Stats<:RuntimeStats}
    max_elem = length(stats)

    serialized_energy = serialize(stats.system_energy, max_elem)
    serialized_power = map(x -> serialize(x.power, max_elem), stats.system_power)
    serialized_time = serialize(stats.step_time, max_elem)

    solver_time = serialize(stats.solver_time, max_elem)
    solver_iterations = serialize(stats.solver_iterations, max_elem)

    logfreq = stats.log_frequency
    current_iter = stats.current_iteration
    store_index = stats.store_index

    parameters = Dict()
    parameters["log_frequency"] = logfreq
    parameters["current_iteration"] = current_iter
    parameters["store_index"] = store_index
    parameters["date"] = Dates.now(UTC)
    parameters["struct"] = "RuntimeStats"
    parameters["float_type"] = string(eltype(solver_time))
    parameters["int_type"] = string(eltype(max_elem))
    parameters["float_vector_type"] = string(typeof(stats.system_energy))
    parameters["int_vector_type"] = string(typeof(stats.solver_iterations))
    parameters["ncomponents"] = length(stats.system_power)
    parameters["samples"] = max_elem

    tojson = Dict()
    tojson["system_energy"] = serialized_energy
    tojson["system_power"] = serialized_power
    tojson["step_time"] = serialized_time
    tojson["solver_time"] = solver_time
    tojson["solver_iterations"] = solver_iterations
    tojson["metadata"] = parameters

    open(path, "w") do io
        JSON.print(io, tojson, 4)
    end
end

function deserialize(::Type{RuntimeStats}, path::String)
    valid_float_types = Dict("Float64" => Float64, "Float32" => Float32)
    valid_int_types = Dict("Int64" => Int64, "Int32" => Int32)
    valid_float_vector_types = Dict("Vector{Float64}" => Vector{Float64},
                                    "Vector{Float32}" => Vector{Float32})
    valid_int_vector_types = Dict("Vector{Int64}" => Vector{Int64},
                                  "Vector{Int32}" => Vector{Int32})

    fromjson = JSON.parsefile(path)

    confparams = fromjson["metadata"]

    if confparams["struct"] != "RuntimeStats"
        error("The file does not contain a RuntimeStats structure")
    end

    FloatVectorType = valid_float_vector_types[confparams["float_vector_type"]]
    IntVectorType = valid_int_vector_types[confparams["int_vector_type"]]
    FloatType = valid_float_types[confparams["float_type"]]
    IntType = valid_int_types[confparams["int_type"]]

    system_energy = deserialize(FloatVectorType, fromjson["system_energy"])

    system_power = Tuple(map(x -> ComponentMass(deserialize(FloatVectorType, x)),
                             fromjson["system_power"]))

    step_time = deserialize(FloatVectorType, fromjson["step_time"])

    solver_time = deserialize(FloatVectorType, fromjson["solver_time"])
    solver_iterations = deserialize(IntVectorType, fromjson["solver_iterations"])

    parameters = fromjson["parameters"]

    logfreq = IntType(parameters["log_frequency"])
    current_iter = IntType(parameters["current_iteration"])
    store_index = IntType(parameters["store_index"])

    RuntimeStats(system_energy, system_power, step_time, solver_time,
                 solver_iterations, logfreq, true, false, current_iter,
                 store_index)
end

export initialize_stats, update_stats!, update_mass!, update_system_energy!,
       system_mass, system_energy, deserialize, serialize, calculate_diff_system_energy,
       calculate_diff_system_mass, start_mass, start_energy, islocked,
       current_iteration, islog, update_component_update_steps!,
       update_system_total_mass!, component_update_calls, component_update_steps,
       calculate_diff_system_total_mass, solver_time, step_time