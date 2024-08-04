"""
    AbstractRuntimeStats{IntType,FloatType,ArrayType}

This abstract type provides basic tracking for runtime of PDE simulation.

# Functions that you need to implement your own `Stats` type:

- `length(Stats)`: Current number of taken samples.
- `islocked(Stats)`: Returns `true` if the `Stats`` is locked meaning that you can not move from sample `n` to `n+1`. However some fields can be updated depending on your design.
- `initialize_stats(::Type{VectorType},ncomponents,log_freq,time_steps,log_solver_info)`: Initialize the memory underlying `Stats` type. 
- `startup_stats!(Stats,start_power,start_energy)`: Initialize the power and energy at time step `0`. This is usally stored at the last index of the memory.
- `update_power!(Stats,power,index)`: Update the power at time step `index`.
- `advance_iteration!(Stats)`: Move from sample  `n`  to  `n+1`.
- `serialize(Stats,path)`: Serialize the `Stats` type to a file.
- `deserialize(Stats,path)`: Deserialize the `Stats` type from a file.


# Optional functions:

- `update_stats!(Stats,time,PDE,Grid,Mem,ItStop)`: For iterative solvers
- `update_stats!(Stats,time,PDE,Grid,Mem)`: For direct solvers

# Examples
"""
abstract type AbstractRuntimeStats{IntType,FloatType,IntArrayType,FloatArrayType} end

struct ComponentPower{FloatType,ArrayType<:AbstractArray{FloatType}} <: AbstractVector{FloatType}
    power::ArrayType
end

mutable struct RuntimeStats{IntType,FloatType,IntArrayType,FloatArrayType<:AbstractArray{FloatType},Power<:Tuple{Vararg{ComponentPower{FloatType,FloatArrayType}}}} <: AbstractRuntimeStats{IntType,FloatType,IntArrayType,FloatArrayType} 
    const system_energy::FloatArrayType
    const system_power::Power
    const step_time::FloatArrayType
    const solver_time::FloatArrayType
    const solver_iterations::IntArrayType
    const log_frequency::IntType
    const locked::Bool
    
    log_data::Bool
    current_iteration::IntType
    store_index::IntType
end