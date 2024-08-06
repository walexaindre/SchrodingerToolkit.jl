Base.length(R::RuntimeStats{IntType,FloatType}) where {IntType,FloatType} = R.store_index -
                                                                            1
Base.size(C::ComponentPower{FloatType,ArrayType}) where {FloatType,ArrayType} = size(C.power)
Base.length(C::ComponentPower{FloatType,ArrayType}) where {FloatType,ArrayType} = length(C.power)

@inline islocked(R::RuntimeStats) = R.locked
@inline islog(R::RuntimeStats) = R.log_data

@inline current_iteration(R::RuntimeStats) = R.current_iteration

@inline start_energy(R::RuntimeStats) = R.system_energy[end]
@inline system_energy(R::RuntimeStats) = R.system_energy
@inline system_energy(R::RuntimeStats, index) = R.system_energy[index]
@inline system_power(R::RuntimeStats, component) = R.system_power[component]
@inline system_power(R::RuntimeStats, component, index) = R.system_power[component][index]
@inline start_power(R::RuntimeStats, component) = R.system_power[component][end]
@inline start_power(R::RuntimeStats) = map(x -> x[end], R.system_power)

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

function initialize_stats(::Type{FloatVectorType}, ::Type{IntVectorType},
                          ncomponents::IntType, log_freq::IntType,
                          time_steps::IntType,
                          islog_solver_info::Bool = true) where {IntType,IntVectorType,
                                                                 FloatVectorType}
    seq_size = div(time_steps, log_freq) + 1
    seq_solver = islog_solver_info ? seq_size : 0
    sys_energy = vundef(FloatVectorType, seq_size)
    sys_power = ntuple(Returns(ComponentPower(similar(sys_energy))), ncomponents)
    sys_time = vundef(FloatVectorType, seq_size)
    solver_time = vzeros(FloatVectorType, seq_solver)
    solver_iterations = vzeros(IntVectorType, seq_solver)

    RuntimeStats(sys_energy, sys_power, sys_time, solver_time, solver_iterations,
                 log_freq, false, false, zero(IntType), one(IntType))
end

function initialize_stats(::Type{FloatVectorType}, ::Type{IntVectorType}, PDE::SPDE,
                          conf::SolverConfig,
                          grid::PGrid) where {FloatVectorType,IntVectorType,
                                              SPDE<:SchrodingerPDE,
                                              PGrid<:PeriodicGrid}
    log_freq = stats_logfreq(conf)
    islog_solver_info = stats_logsol(conf)
    ncomp = ncomponents(PDE)
    tsteps = estimate_timesteps(PDE, grid)
    initialize_stats(FloatVectorType, IntVectorType, ncomp, log_freq, tsteps,
                     islog_solver_info)
end

function startup_stats!(stats::Stats, start_power,
                        start_energy) where {Stats<:RuntimeStats}
    last_index = lastindex(stats.step_time)
    update_power!(stats, start_power, last_index)
    update_system_energy!(stats, start_energy, last_index)
    advance_iteration!(stats)
end

@inline function advance_iteration!(R::RuntimeStats{IntType,FloatType}) where {IntType<:Integer,
                                                                               FloatType<:AbstractFloat}
    if !islocked(R)
        R.current_iteration += 1
        if mod(R.current_iteration, R.log_frequency) == 0
            R.log_data = true
        end
    else
        error("Stats are locked")
    end
end

function update_power!(stats::Stats, power,
                       idx::IntType) where {IntType<:Integer,
                                            FloatType<:AbstractFloat,
                                            Stats<:RuntimeStats{IntType,FloatType}}
    power = Vector(power)

    if !islocked(stats)
        for comp in 1:length(stats.system_power)
            stats.system_power[comp][idx] = power[comp]
        end
    else
        error("Stats are locked")
    end
end

function update_system_energy!(stats::Stats, energy::FloatType,
                               idx::IntType) where {IntType<:Integer,
                                                    FloatType<:AbstractFloat,
                                                    Stats<:RuntimeStats{IntType,
                                                                        FloatType}}
    if !islocked(stats)
        stats.system_energy[idx] = energy
    else
        error("Stats are locked")
    end
end

"""
    update_solver_info!(stats, time, iterations)

    Update the linear solver information in the stats structure.

    # Arguments
    - `stats::Stats`: The stats structure.
    - `time::FloatType`: The time spent in the solver.
    - `iterations::IntType`: The number of iterations performed in the solver.
"""
function update_solver_info!(stats::Stats, time,
                             iterations) where {Stats<:RuntimeStats}
    if !islocked(stats)
        if stats.log_data && length(stats.solver_time) > 0
            idx = stats.store_index
            stats.solver_time[idx] += time
            stats.solver_iterations[idx] += iterations
        end
    else
        error("Stats are locked")
    end
end

"""
    calculate_diff_system_energy(stats)

    Calculate the absolute value of the difference between the system energy and the energy at time step `0`.

    # Arguments
    - `stats::Stats`: The stats structure.
"""
function calculate_diff_system_energy(stats::Stats) where {Stats<:RuntimeStats}
    startup_energy = start_energy(stats)

    abs.(system_energy(stats)[1:length(stats)] .- startup_energy)
end

"""
    calculate_diff_system_power(stats)

    Calculate the absolute value of the difference between the system power and the power at time step `0`.

    # Arguments
    - `stats::Stats`: The stats structure.
"""
function calculate_diff_system_power(stats::Stats,
                                     index) where {Stats<:RuntimeStats}
    startup_power = start_power(stats, index)

    return abs.(system_power(stats, index)[1:length(stats)] .- startup_power)
end

"""
    update_stats!(stats, step_time, power_per_component, sys_energy)

    Update the stats structure with the current time, power per component and system energy.

    # Arguments
    - `stats::Stats`: The stats structure.
    - `step_time::FloatType`: Amount of time spent in the current step.
    - `power_per_component`: The power per component.
    - `sys_energy::FloatType`: The system energy.
"""
function update_stats!(stats::Stats, step_time, power_per_component,
                       sys_energy) where {Stats<:RuntimeStats}
    if !islocked(stats)
        if stats.log_data
            idx = stats.store_index
            stats.step_time[idx] = step_time
            update_power!(stats, power_per_component, idx)
            update_system_energy!(stats, sys_energy, idx)
            stats.store_index += 1
            stats.log_data = false
        end
        advance_iteration!(stats)
    else
        error("Stats are locked")
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
       calculate_diff_system_power, start_power, start_energy, islocked,current_iteration,islog