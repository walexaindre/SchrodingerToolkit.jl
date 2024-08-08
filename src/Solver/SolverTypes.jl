struct SolverConfig{N,IntType,Backend}
    time_order::Symbol
    space_order::NTuple{N,Symbol}

    backend_type::Type{Backend}

    verbose::IntType

    show_progress::Bool
    output_folder_path::String

    #stopping_criteria::StoppingCriteria

    # Index: 1 = log_frequency, 2 = log_solver_info, 3 = log_system_total_power, 4 = log_component_update_steps, 5 = locked, 6 = log_stats,
    stats::Tuple{IntType,Vararg{Bool,5}}
    debug::Tuple{Bool,IntType,IntType}
end

struct SchrodingerProblem{Method,Storage,Statistics,SolverConf,PDEq}
    Method::Method
    Memory::Storage
    Stats::Statistics
    PDE::PDEq
    Config::SolverConf
end

export SchrodingerProblem, SolverConfig