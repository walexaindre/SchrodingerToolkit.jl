"""
    AbstractRuntimeStats{IntType,FloatType,ArrayType}

This abstract type provides basic tracking for runtime analysis of PDE simulation.

# Functions that you need to implement your own `Stats` type:

- `length(Stats)`: Current number of taken samples.

- `islog(Stats)`: Returns `true` if the `Stats` is logging problem related data.
- `islocked(Stats)`: Returns `true` if the `Stats` is locked meaning that you can not modify it. This is used to prevent data corruption after deserialization.
- `islog_solver_info(Stats)`: Returns `true` if the `Stats` is logging solver related data.
- `islog_system_total_power(Stats)`: Returns `true` if the `Stats` is logging the total power of the system.
- `islog_component_update_steps(Stats)`: Returns `true` if the `Stats` is logging the number of times that the components are updated.

- `initialize_stats(::Type{FloatVectorType}, ::Type{IntVectorType},ncomponents,log_freq,time_steps,islog_solver_info::Bool)`: Initialize the memory underlying `Stats` type. 
- `startup_stats!(Stats,start_power,start_energy)`: Initialize the power and energy at time step `0`. This is usally stored at the last index of the memory.
- `update_power!(Stats,power,index)`: Update the power at time step `index`.
- `update_energy!(Stats,energy,index)`: Update the energy at time step `index`.
- `advance_iteration!(Stats)`: Move from sample  `n`  to  `n+1`.

- `serialize(Stats,path)`: Serialize the `Stats` type to a file.
- `deserialize(::Type{Stats},path)`: Deserialize the `Stats` type from a file. The type here is for dispatching purposes.

- `system_power(Stats,component_idx,index)`: Get the power at `component_idx` in time step `index`.
- `system_power(Stats,component_idx)`: Get the full internal power vector at `component_idx`.
- `system_energy(Stats,index)`: Get the energy at time step `index`.
- `system_energy(Stats)`: Get the full internal energy vector.
- `start_power(Stats)`: Get the power for all components at time step `0`.
- `start_power(Stats,component_idx)`: Get the power at `component_idx` at time step `0`.
- `start_energy(Stats)`: Get the energy at time step `0`.
- `get_solver_time(Stats,index)`: Get the solver time at time step `index`.
- `get_solver_iterations(Stats,index)`: Get the solver iterations at time step `index`.

Note `1`: Any attempt to modify the `Stats` type after it is locked must throw an error or a warning.

Note `2`: In general is preferred to use CPU vectors for `Stats` type and taking care of the conversion to other types when needed. This is because indexing is slow and disallowed in GPU arrays.

Note `3`: As an API decision, to reduce boilerplate code is expected that when you intialize any method, then this structure is initialized too.
This is because the `Stats` type is used to store incoming data from the simulation and data types are known at that moment. 

# Optional functions:

- `update_stats!(Stats,time,PDE,Grid,Mem,ItStop)`: For iterative solvers
- `update_stats!(Stats,time,PDE,Grid,Mem)`: For direct solvers

# Examples
"""
abstract type AbstractRuntimeStats{IntType,FloatType,IntArrayType,FloatArrayType} end

struct ComponentPower{FloatType,ArrayType<:AbstractArray{FloatType}} <:
       AbstractVector{FloatType}
    power::ArrayType
end

mutable struct RuntimeStats{IntType,FloatType,IntVectorType,
                            FloatVectorType<:AbstractArray{FloatType},
                            Power<:Tuple{Vararg{ComponentPower{FloatType,FloatVectorType}}}} <:
               AbstractRuntimeStats{IntType,FloatType,IntVectorType,FloatVectorType}
    const system_energy::FloatVectorType
    const system_power::Power

    const system_total_power::FloatVectorType

    const step_time::FloatVectorType

    const solver_time::FloatVectorType
    const solver_iterations::IntVectorType

    const component_update_steps::IntVectorType
    const log_frequency::IntType
    
    # Index: 1 = log_solver_info, 2 = log_system_total_power, 3 = log_component_update_steps, 4 = locked, 5 = log_stats, 6 - 8 = unused
    const config::NTuple{8,Bool}

    log_data::Bool

    component_update_call_count::IntType
    current_iteration::IntType
    store_index::IntType
end