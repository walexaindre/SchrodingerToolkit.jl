Base.length(R::RuntimeStats{IntType,FloatType}) where {IntType,FloatType} = R.store_index -
                                                                            1
Base.size(C::ComponentPower{FloatType,ArrayType}) where {FloatType,ArrayType} = size(C.power)
Base.length(C::ComponentPower{FloatType,ArrayType}) where {FloatType,ArrayType} = length(C.power)

@inline current_iteration(R::RuntimeStats) = R.current_iteration

@inline start_energy(R::RuntimeStats) = R.system_energy[end]
@inline system_energy(R::RuntimeStats) = R.system_energy
@inline system_energy(R::RuntimeStats, index) = R.system_energy[index]
@inline system_power(R::RuntimeStats, component) = R.system_power[component]
@inline system_power(R::RuntimeStats, component, index) = R.system_power[component][index]
@inline system_total_power(R::RuntimeStats) = R.system_total_power
@inline start_power(R::RuntimeStats, component) = R.system_power[component][end]
@inline start_power(R::RuntimeStats) = map(x -> x[end], R.system_power)
@inline start_total_power(R::RuntimeStats) = R.system_total_power[end]

@inline component_update_steps(R::RuntimeStats) = R.component_update_steps
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

@inline islog_system_total_power(R::RuntimeStats) = R.config[2]
@inline islog_solver_info(R::RuntimeStats) = R.config[1]
@inline islog_component_update_steps(R::RuntimeStats) = R.config[3]

@inline function Base.getindex(C::ComponentPower{FloatType,
                                                 ArrayType},
                               index) where {FloatType,ArrayType}
    C.power[index]
end

@inline function Base.setindex!(C::ComponentPower{FloatType,ArrayType},
                                value::FloatType,
                                index) where {FloatType,ArrayType}
    C.power[index] = value
end

"""
    initialize_stats(FloatVectorType, IntVectorType, ncomponents, log_freq, time_steps;
                     log_stats = true, log_solver_info = true,
                     log_component_update_steps = false,
                     log_system_total_power = false, islocked = false)

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
- `log_system_total_power`: Log the system total power.
- `islocked`: Lock the stats.

## Returns

- `RuntimeStats`: The stats structure.

"""
function initialize_stats(::Type{FloatVectorType}, ::Type{IntVectorType},
                          ncomponents::IntType, log_freq::IntType,
                          time_steps::IntType; log_stats::Bool = true,
                          log_solver_info::Bool = true,
                          log_component_update_steps::Bool = false,
                          log_system_total_power::Bool = false,
                          islocked::Bool = false) where {IntType,
                                                         IntVectorType,
                                                         FloatVectorType}
    config = (log_solver_info, log_system_total_power, log_component_update_steps,
              islocked, log_stats, false, false, false)

    seq_size = div(time_steps, log_freq) + 1
    seq_solver = log_solver_info ? seq_size : 0
    seq_total_power = log_system_total_power ? seq_size : 0
    seq_component_update_steps = log_component_update_steps ? seq_size : 0

    sys_energy = vundef(FloatVectorType, seq_size)
    sys_power = ntuple(Returns(ComponentPower(similar(sys_energy))), ncomponents)
    sys_total_power = vundef(FloatVectorType, seq_total_power)
    sys_time = vundef(FloatVectorType, seq_size)
    solver_time = vzeros(FloatVectorType, seq_solver)
    solver_iterations = vzeros(IntVectorType, seq_solver)
    component_update_steps = vzeros(IntVectorType, seq_component_update_steps)

    RuntimeStats(sys_energy, sys_power, sys_total_power, sys_time, solver_time,
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
    log_system_total_power = stats_log_system_total_power(conf)
    log_solver_info = stats_log_solver_info(conf)
    islocked = stats_locked(conf)
    log_stats = stats_log_stats(conf)

    ncomp = ncomponents(PDE)
    tsteps = estimate_timesteps(PDE, grid)
    initialize_stats(FloatVectorType, IntVectorType, ncomp, log_freq, tsteps;
                     log_stats = log_stats,
                     log_solver_info = log_solver_info,
                     log_component_update_steps = log_component_update_steps,
                     log_system_total_power = log_system_total_power,
                     islocked = islocked)
end

"""
    startup_stats!(stats::Stats, start_power, start_energy)

Initialize the stats structure with the power and energy at time step `0`.
"""
function startup_stats!(stats::Stats, start_power,
                        start_energy) where {Stats<:RuntimeStats}
    last_index = lastindex(stats.step_time)
    update_power!(stats, start_power, last_index)
    update_system_energy!(stats, start_energy, last_index)
    advance_iteration!(stats)
end

"""
    startup_stats!(stats::Stats, start_power, start_energy, total_power)

Initialize the stats structure with the power, energy and total power at time step `0`.
"""
function startup_stats!(stats::Stats, start_power, start_energy,
                        total_power) where {Stats<:RuntimeStats}
    last_index = lastindex(stats.step_time)
    update_power!(stats, start_power, last_index)
    update_system_energy!(stats, start_energy, last_index)
    update_system_total_power!(stats, total_power, last_index)
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
    update_power!(stats, power, idx)

Update the power per component in the stats structure.

## Arguments
- `stats::Stats`: The stats structure.
- `power::Vector{FloatType}`: The power per component.
- `idx::IntType`: The index to update.
"""
function update_power!(stats::Stats, power,
                       idx::IntType) where {IntType<:Integer,
                                            FloatType<:AbstractFloat,
                                            Stats<:RuntimeStats{IntType,FloatType}}
    power = Vector(power)

    if !islocked(stats) && isregistering_stats(stats)
        for comp in 1:length(stats.system_power)
            stats.system_power[comp][idx] = power[comp]
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
    update_system_total_power!(stats, total_power, idx)

Update the system total power in the stats structure.

## Arguments
- `stats::Stats`: The stats structure.
- `total_power::FloatType`: The total power.
- `idx::IntType`: The index to update.
"""
function update_system_total_power!(stats::Stats, total_power::FloatType,
                                    idx::IntType) where {IntType<:Integer,
                                                         FloatType<:AbstractFloat,
                                                         Stats<:RuntimeStats{IntType,
                                                                             FloatType}}
    if !islocked(stats) && islog_system_total_power(stats) && isregistering_stats(stats)
        stats.system_total_power[idx] = total_power
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
    calculate_diff_system_power(stats)

Calculate the absolute value of the difference between the system power and the power at time step `0`.

## Arguments
- `stats::Stats`: The stats structure.
"""
function calculate_diff_system_power(stats::Stats,
                                     index) where {Stats<:RuntimeStats}
    startup_power = start_power(stats, index)

    return abs.(system_power(stats, index)[1:length(stats)] .- startup_power)
end

"""
    calculate_diff_system_total_power(stats)

Calculate the absolute value of the difference between the system total power and the total power at time step `0`.

## Arguments
- `stats::Stats`: The stats structure.
"""
function calculate_diff_system_total_power(stats::Stats) where {Stats<:RuntimeStats}
    startup_power = start_total_power(stats)

    abs.(system_total_power(stats)[1:length(stats)] .- startup_power)
end

"""
    update_stats!(stats, step_time, power_per_component, sys_energy)

Update the stats structure with the current time, power per component and system energy.

## Arguments
- `stats::Stats`: The stats structure.
- `step_time::FloatType`: Amount of time spent in the current step.
- `power_per_component`: The power per component.
- `sys_energy::FloatType`: The system energy.
"""
function update_stats!(stats::Stats, step_time, power_per_component,
                       sys_energy) where {Stats<:RuntimeStats}
    if !islocked(stats) && isregistering_stats(stats)
        if islog(stats)
            idx = stats.store_index
            stats.step_time[idx] = step_time
            update_power!(stats, power_per_component, idx)
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
    update_stats!(stats, step_time, power_per_component, sys_energy, total_power)

Update the stats structure with the current time, power per component, system energy and total power.

## Arguments
- `stats::Stats`: The stats structure.
- `step_time::FloatType`: Amount of time spent in the current step.
- `power_per_component`: The power per component.
- `sys_energy::FloatType`: The system energy.
- `total_power::FloatType`: The total power.
"""
function update_stats!(stats::Stats, step_time, power_per_component, sys_energy,
                       total_power) where {Stats<:RuntimeStats}
    if !islocked(stats) && isregistering_stats(stats)
        if islog(stats)
            idx = stats.store_index
            stats.step_time[idx] = step_time
            update_power!(stats, power_per_component, idx)
            update_system_energy!(stats, sys_energy, idx)
            update_system_total_power!(stats, total_power, idx)
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

Update the stats structure with the current time, power per component, system energy and total power.

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
function update_stats!(stats::Stats, memory, grid, PDE,
                       work_timer) where {Stats<:RuntimeStats}
    if isregistering_stats(stats)
        if islog(stats)
            power_per_component = system_power(memory, grid)
            energy = system_energy(memory, PDE, grid)
            if islog_system_total_power(conf)
                total_power = system_total_power(memory, PDE, grid,
                                                 power_per_component)
                update_stats!(stats, work_timer, power_per_component, energy, total_power)
            end
            update_stats!(stats, work_timer, power_per_component, energy)
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

    system_power = Tuple(map(x -> ComponentPower(deserialize(FloatVectorType, x)),
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

export initialize_stats, update_stats!, update_power!, update_system_energy!,
       system_power, system_energy, deserialize, serialize, calculate_diff_system_energy,
       calculate_diff_system_power, start_power, start_energy, islocked,
       current_iteration, islog, update_component_update_steps!,
       update_system_total_power!, component_update_calls, component_update_steps,
       calculate_diff_system_total_power