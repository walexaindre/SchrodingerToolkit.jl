struct SolverConfig{N,IntType,Backend}
    time_order::Symbol
    space_order::NTuple{N,Symbol}

    backend_type::Type{Backend}
    
    verbose::IntType
    
    show_progress::Bool
    output_folder_path::String

    #stopping_criteria::StoppingCriteria

    stats::Tuple{Bool,IntType}
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