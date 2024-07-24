"""
    AbstractStats{IntType,FloatType,ArrayType}

This abstract type provides basic tracking for runtime of PDE simulation.

# Functions that you need to implement your own `Stats` type:

- `length(Stats)`: Current number of taken samples.
- `islocked(Stats)`: Returns `true` if the `Stats`` is locked meaning that you can not move from sample `n` to `n+1`. However some fields can be updated depending on your design.
- `initialize_stats(Stats,ncomponents,log_freq,τ,timesteps,log_solverinfo)`: Here you should initialize the `Stats` type. `ncomponents` is the number of components in the PDE, `log_freq` is the frequency of logging, `τ` is the time step, `timesteps` is the total number of time steps, and `log_solverinfo` is a boolean that indicates whether to log solver information.
- `update_system_energy(Stats,energy,index)`: Update the system energy at time step `index`.
- `update_power!(Stats,power,index)`: Update the power at time step `index`.
- `advance!(Stats)`: Move from sample `n` to `n+1`.
- `serialize(Stats,path)`: Serialize the `Stats` type to a file.
- `deserialize(Stats,path)`: Deserialize the `Stats` type from a file.


# Optional functions:

- `update_stats!(Stats,time,PDE,Grid,Mem,ItStop)`

# Examples
"""
abstract type AbstractStats{IntType,FloatType,ArrayType} end

