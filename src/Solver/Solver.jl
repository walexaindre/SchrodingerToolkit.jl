Base.show(io::IO, conf::SolverConfig) = print(io,
                                              "SolverConfig:\n\
                                              \n\ttime order:        $(conf.time_order),\
                                              \n\tspace order:       $(conf.space_order),\
                                              \n\tbackend type:      $(conf.backend_type),\
                                              \n\tverbose:           $(conf.verbose),\
                                              \n\tshow_progress:     $(conf.show_progress),\
                                              \n\toutput_folder:     $(conf.output_folder_path),\
                                              \n\tstats:             $(conf.stats),\
                                              \n\tdebug:             $(conf.debug))")

function SolverParameters(::Type{ComputeBackend}, dim::IntType,
                          spaceorder = ntuple(Returns(:ord2), dim),
                          timeorder = :tord2_1_1;
                          verbose::IntType = zero(IntType),
                          stats::Tuple{Bool,IntType,Bool} = (true, IntType(5),true),
                          debug::Tuple{Bool,IntType,IntType} = (false, zero(IntType),
                                                                zero(IntType)),
                          show_progress::Bool = false,
                          output_folder::String = "/data") where {IntType,
                                                                  ComputeBackend<:AbstractBackend{IntType}}
    SolverConfig(timeorder,
                 spaceorder,
                 ComputeBackend,
                 verbose,
                 show_progress,
                 output_folder,
                 stats,
                 debug)
end


Base.ndims(conf::SolverConfig) = length(conf.space_order)

backend_type(conf::SolverConfig) = conf.backend_type
debug(conf::SolverConfig) = conf.debug[1]
debug_level(conf::SolverConfig) = conf.debug[2]
debug_sublevel(conf::SolverConfig) = conf.debug[3]
stats(conf::SolverConfig) = conf.stats[1]
stats_logfreq(conf::SolverConfig) = conf.stats[2]
stats_logsol(conf::SolverConfig) = conf.stats[3]
verbose(conf::SolverConfig) = conf.verbose
show_progress(conf::SolverConfig) = conf.show_progress
output_folder(conf::SolverConfig) = conf.output_folder_path
time_order(conf::SolverConfig) = conf.time_order
space_order(conf::SolverConfig) = conf.space_order

method(problem::SchrodingerProblem) = problem.Method
memory(problem::SchrodingerProblem) = problem.Memory
stats(problem::SchrodingerProblem) = problem.Stats
PDE(problem::SchrodingerProblem) = problem.PDE
config(problem::SchrodingerProblem) = problem.Config

"""
    With this function you can advance one step in time.

        The implementation of this function is mandatory.

        From the user perspective, this function should be called to advance the solution one step in time.

        The expected parameters are:
            - `method`: The method to be used to advance the solution.
            - `memory`: The memory of the method.
                Expectations about memory are: 
                    - The memory must be allocated by the method at initialization.
                    - The memory must be updated by the method.
                    - You must provide a method called `current_state!(memory)` that allows the caller to a reference to the current state (time step).
                    - You must provide a method called `current_state(memory)` that allows the caller to get a copy of the current state of (time step).
            - `stats`: The statistics of the method.
                    - If you want to log statistics and you perform some kind of component wise operation then:
                        - You must call `update_solver_info!(stats, time, niter)` to store every time spent at linear solvers.
                    - In general 
            - `PDE`: The PDE to be solved.
            - `config`: The configuration of the solver.
        
        Here dispatch is used to call the correct version of step!.

        The return of this function must be the time taken to advance one step in time.
"""
function step! end

@inline step!(P::SchrodingerProblem) = step!(method(P), memory(P),stats(P), PDE(P), config(P))


"""

    With this function you can update a component of the problem.

        The implementation of this function is suggested but not mandatory.

"""
function update_component! end
    

export  SolverParameters, step!, update_component!