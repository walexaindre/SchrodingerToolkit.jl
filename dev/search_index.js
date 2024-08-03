var documenterSearchIndex = {"docs":
[{"location":"Theory/DFP/","page":"DFP","title":"DFP","text":"In the journey of explore how to improve parallelism we need to learn the theory behind how to discretize the equation.","category":"page"},{"location":"Theory/Problem/","page":"Problem","title":"Problem","text":"Let Omega be a Cartesian region in mathbbR^m with m=123. We are interested in get an approximation of the complex field Psi Omega times (0T longrightarrow mathbbC^N satisfying the following system of non linear equations in dimensionless form:","category":"page"},{"location":"Theory/Problem/","page":"Problem","title":"Problem","text":"    i Psi_t = D Delta Psi + mathcalV(mathbfx) odot Psi+mathcalJ(Psi)","category":"page"},{"location":"Theory/Problem/","page":"Problem","title":"Problem","text":"where mathbfx = (x_1cdotsx_m), Psi^2=left( psi_1^2cdotspsi_2^2 right), D=diag(d_i) is a mathbbR^Ntimes N diagonal matrix with positive dispersive coefficients, odot is the Hadamard or element-wise product, F mathbbR^N longrightarrow mathbbR, mathcalV Omega longrightarrow mathbbR^N represents an external trapping potential. Finally, the linear functional mathcalJ mathbbC^N longrightarrow mathbbC^N models the internal atomic Josephson junction, and takes the form mathcalJ(Psi)_n = Gammasum_j=1jneq n^N psi_j, where Gamma is a real coefficient","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SchrodingerToolkit","category":"page"},{"location":"#SchrodingerToolkit","page":"Home","title":"SchrodingerToolkit","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SchrodingerToolkit.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SchrodingerToolkit]","category":"page"},{"location":"#SchrodingerToolkit.AbstractBackend","page":"Home","title":"SchrodingerToolkit.AbstractBackend","text":"Abstract type for our available backends\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.AbstractRuntimeStats","page":"Home","title":"SchrodingerToolkit.AbstractRuntimeStats","text":"AbstractRuntimeStats{IntType,FloatType,ArrayType}\n\nThis abstract type provides basic tracking for runtime of PDE simulation.\n\nFunctions that you need to implement your own Stats type:\n\nlength(Stats): Current number of taken samples.\nislocked(Stats): Returns true if the Statsis locked meaning that you can not move from samplenton+1`. However some fields can be updated depending on your design.\ninitialize_stats(Stats,ncomponents,log_freq,τ,timesteps,log_solverinfo): Here you should initialize the Stats type. ncomponents is the number of components in the PDE, log_freq is the frequency of logging, τ is the time step, timesteps is the total number of time steps, and log_solverinfo is a boolean that indicates whether to log solver information.\nupdate_system_energy!(Stats,energy,index): Update the system energy at time step index.\nupdate_power!(Stats,power,index): Update the power at time step index.\nadvance!(Stats): Move from sample n to n+1.\nserialize(Stats,path): Serialize the Stats type to a file.\ndeserialize(Stats,path): Deserialize the Stats type from a file.\n\nOptional functions:\n\nupdate_stats!(Stats,time,PDE,Grid,Mem,ItStop): For iterative solvers\nupdate_stats!(Stats,time,PDE,Grid,Mem): For direct solvers\n\nExamples\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.CPUBackend","page":"Home","title":"SchrodingerToolkit.CPUBackend","text":"This represents a CPU based backend\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.GPUBackend","page":"Home","title":"SchrodingerToolkit.GPUBackend","text":"This represents a GPU based backend\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.OffsetAllZero","page":"Home","title":"SchrodingerToolkit.OffsetAllZero","text":"No dimension has the zero offset\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.OffsetNonZero","page":"Home","title":"SchrodingerToolkit.OffsetNonZero","text":"Every dimension has the zero offset\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.OffsetUniqueZero","page":"Home","title":"SchrodingerToolkit.OffsetUniqueZero","text":"Only the first dimension has the zero offset\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.SchrodingerPDE","page":"Home","title":"SchrodingerToolkit.SchrodingerPDE","text":"SchrodingerPDE\n\nThis is an abstract type that represents a Schrodinger PDE. This type is used to define a Schrodinger PDE structure.\n\nParameters:\n\n- `N::Int`: Number of dimensions\n- `RealType`: Type for Real numbers. This is used to define the internal algorithms and data structures.\n\nIs important to note, that the derived structures must handle boundary conditions, initial conditions and important parts of the PDE.\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.SchrodingerPDEComponent","page":"Home","title":"SchrodingerToolkit.SchrodingerPDEComponent","text":"SchrodingerPDEComponent\n\nThis structure represents a component of the Schrodinger PDE. This is a simple structure that contains the necessary information\nto solve the Schrodinger PDE for a single component.\n\nFields:\n\n- `σ::Tv`: Dispersion coefficient => Must be non negative\n- `f::Fn`: ∂∇F/∂xᵢ(x) = f(x)\n- `ψ::InitialCondition`: Initial condition \n- `V::TrappingPotential`: External trapping potential\n- `Γ::Γv`: Josephson Junction Coefficient => In some documents this coefficient is the same for all components but I guess it can be different under some circumstances.\n\nBy default `V` and `Γ` are set to `nothing` but they can be set to a trapping potential and a Josephson junction coefficient respectively if you want to work with a Gross Pitaevsky PDE.\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.SchrodingerPDENonPolynomic","page":"Home","title":"SchrodingerToolkit.SchrodingerPDENonPolynomic","text":"SchrodingerPDENonPolynomic\n\nThis structure represents a Schrodinger PDE with non-polynomical potentials.\n\nFields:\n\n- `boundaries::NTuple{N,NTuple{2,Tv}}`: Boundaries of the cartesian domain where you want to solve the PDE.\n- `components::MComp`: Components of the Schrodinger PDE\n- `F::Potential`: Non-polynomical potential\n- `T::Tv`: Final time\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.SchrodingerPDEPolynomic","page":"Home","title":"SchrodingerToolkit.SchrodingerPDEPolynomic","text":"SchrodingerPDEPolynomic\n\nThis structure represents a Schrodinger PDE with polynomical potentials.\n\nFields:\n\n- `boundaries::NTuple{N,NTuple{2,Tv}}`: Boundaries of the cartesian domain where you want to solve the PDE.\n- `components::MComp`: Components of the Schrodinger PDE\n- `F::Potential`: Polynomical potential\n- `N::Optimized`: Optimized structure for polynomical potentials. Here is stored an auxiliary function without the catastrophic cancellation problem.\n- `T::Tv`: Final time\n\n\n\n\n\n","category":"type"},{"location":"#SchrodingerToolkit.xPUBackend","page":"Home","title":"SchrodingerToolkit.xPUBackend","text":"This represents a heterogeneous based backend (Can use a mixture of CPU and GPU)\n\n\n\n\n\n","category":"type"},{"location":"#Base.eltype-Union{Tuple{Type{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}}}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential}","page":"Home","title":"Base.eltype","text":"Type of the elements for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#Base.eltype-Union{Tuple{Type{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}}}, Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"Base.eltype","text":"Type of the elements for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#Base.ndims-Union{Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential}","page":"Home","title":"Base.ndims","text":"Number of dimensions for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#Base.ndims-Union{Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}}, Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"Base.ndims","text":"Number of dimensions for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.estimate_timesteps-Union{Tuple{PGrid}, Tuple{SPDE}, Tuple{SPDE, PGrid}} where {SPDE<:SchrodingerPDE, PGrid<:PeriodicGrid}","page":"Home","title":"SchrodingerToolkit.estimate_timesteps","text":"Estimate the number of timesteps for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.evaluate_ψ!-Union{Tuple{PGrid}, Tuple{SPDE}, Tuple{SPDE, PGrid, Any}} where {SPDE<:SchrodingerPDE, PGrid<:PeriodicGrid}","page":"Home","title":"SchrodingerToolkit.evaluate_ψ!","text":"Evaluation of the initial condition for the Schrodinger PDE where the output is stored in Container\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.evaluate_ψ-Union{Tuple{PGrid}, Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}, PGrid, Any}} where {N, Tv, MComp, Potential, Optimized, PGrid<:PeriodicGrid}","page":"Home","title":"SchrodingerToolkit.evaluate_ψ","text":"Evaluation of the initial condition for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.evaluate_ψ-Union{Tuple{PGrid}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}, PGrid, Any}} where {N, Tv, MComp, Potential, PGrid<:PeriodicGrid}","page":"Home","title":"SchrodingerToolkit.evaluate_ψ","text":"Evaluation of the initial condition for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_V-Union{Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}, Int64}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_V","text":"Trapping potential for the i-th Schrodinger PDE component\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_V-Union{Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}, Int64}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_V","text":"Trapping potential for the i-th Schrodinger PDE component\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_boundary-Union{Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}, Int64}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_boundary","text":"Boundary start and end point at index dimension\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_boundary-Union{Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}, Int64}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_boundary","text":"Boundary start and end point at index dimension\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_component-Union{Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}, Int64}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_component","text":"Obtain the component of the Schrodinger PDE at index\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_component-Union{Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}, Int64}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_component","text":"Obtain the component of the Schrodinger PDE at index\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_f-Union{Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}, Int64}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_f","text":"Function for the Schrodinger PDE component\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_f-Union{Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}, Int64}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_f","text":"Function for the Schrodinger PDE component\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_field-Union{Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_field","text":"Potential for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_field-Union{Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}}, Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_field","text":"Potential for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_optimized-Union{Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_optimized","text":"Non polynomial potentials can't be optimized so return nothing\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_optimized-Union{Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}}, Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_optimized","text":"Simplified polynomial to evade catastrophic cancellation\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_time_boundary-Union{Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_time_boundary","text":"Finish time for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_time_boundary-Union{Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}}, Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_time_boundary","text":"Finish time for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_Γ-Union{Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}, Int64}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_Γ","text":"Josephson Junction Coefficient for the i-th Schrodinger PDE component\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_Γ-Union{Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}, Int64}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_Γ","text":"Josephson Junction Coefficient for the i-th Schrodinger PDE component\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_σ-Union{Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}, Int64}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_σ","text":"Dispersion coefficient for the Schrodinger PDE component at index\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_σ-Union{Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}, Int64}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_σ","text":"Dispersion coefficient for the Schrodinger PDE component at index\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_ψ-Union{Tuple{Optimized}, Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDEPolynomic{N, Tv, MComp, Potential, Optimized}, Int64}} where {N, Tv, MComp, Potential, Optimized}","page":"Home","title":"SchrodingerToolkit.get_ψ","text":"Initial condition for the Schrodinger PDE component\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.get_ψ-Union{Tuple{Potential}, Tuple{MComp}, Tuple{Tv}, Tuple{N}, Tuple{SchrodingerPDENonPolynomic{N, Tv, MComp, Potential}, Int64}} where {N, Tv, MComp, Potential}","page":"Home","title":"SchrodingerToolkit.get_ψ","text":"Initial condition for the Schrodinger PDE component\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.iscpu","page":"Home","title":"SchrodingerToolkit.iscpu","text":"Function to check if a backend is a CPU based backend\n\n\n\n\n\n","category":"function"},{"location":"#SchrodingerToolkit.isgpu","page":"Home","title":"SchrodingerToolkit.isgpu","text":"Function to check if a backend is a GPU based backend\n\n\n\n\n\n","category":"function"},{"location":"#SchrodingerToolkit.isxpu","page":"Home","title":"SchrodingerToolkit.isxpu","text":"Function to check if a backend is a heterogeneous based backend\n\n\n\n\n\n","category":"function"},{"location":"#SchrodingerToolkit.matrix_to_vector-Union{Tuple{Type{M}}, Tuple{M}} where M<:(DenseMatrix)","page":"Home","title":"SchrodingerToolkit.matrix_to_vector","text":"S = matrix_to_vector(M)\n\nReturn the dense vector storage type S related to the dense matrix storage type M.\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.measure-Union{Tuple{PeriodicGrid{V, T, R, N}}, Tuple{N}, Tuple{R}, Tuple{T}, Tuple{V}} where {V<:Integer, T<:Real, R<:AbstractRange{T}, N}","page":"Home","title":"SchrodingerToolkit.measure","text":"measure(A::PeriodicGrid{V,T,R,N}) where {V<:Integer,T<:Real,R<:AbstractRange{T},N}\n\nReturn the measure of the grid.\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.ncomponents-Tuple{PDEeq} where PDEeq<:SchrodingerPDE","page":"Home","title":"SchrodingerToolkit.ncomponents","text":"Number of components for the Schrodinger PDE\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.step!","page":"Home","title":"SchrodingerToolkit.step!","text":"With this function you can advance one step in time.\n\n    The implementation of this function is mandatory.\n\n    From the user perspective, this function should be called to advance the solution one step in time.\n\n    The expected parameters are:\n        - `method`: The method to be used to advance the solution.\n        - `memory`: The memory of the method.\n            Expectations about memory are: \n                - The memory must be allocated by the method at initialization.\n                - The memory must be updated by the method.\n                - You must provide a method called `current_state!(memory)` that allows the caller to a reference to the current state (time step).\n                - You must provide a method called `current_state(memory)` that allows the caller to get a copy of the current state of (time step).\n        - `stats`: The statistics of the method.\n                - If you want to log statistics and you perform some kind of component wise operation then:\n                    - You must call `update_solver_info!(stats, time, niter)` to store every time spent at linear solvers.\n        - `PDE`: The PDE to be solved.\n        - `config`: The configuration of the solver.\n    \n    Here dispatch is used to call the correct version of step!.\n\n    The return of this function must be the time taken to advance one step in time.\n\n\n\n\n\n","category":"function"},{"location":"#SchrodingerToolkit.test_if_zero","page":"Home","title":"SchrodingerToolkit.test_if_zero","text":"v = test_if_zero(V,idx)\n\nReturn true if the element at position idx of the vector V is zero.\n\n\n\n\n\n","category":"function"},{"location":"#SchrodingerToolkit.update_component!","page":"Home","title":"SchrodingerToolkit.update_component!","text":"With this function you can update a component of the problem.\n\n    The implementation of this function is suggested but not mandatory.\n\n\n\n\n\n","category":"function"},{"location":"#SchrodingerToolkit.update_solver_info!-Union{Tuple{Stats}, Tuple{Stats, Any, Any}} where Stats<:SchrodingerToolkit.RuntimeStats","page":"Home","title":"SchrodingerToolkit.update_solver_info!","text":"update_solver_info!(stats, time, iterations)\n\nUpdate the linear solver information in the stats structure.\n\n# Arguments\n- `stats::Stats`: The stats structure.\n- `time::FloatType`: The time spent in the solver.\n- `iterations::IntType`: The number of iterations performed in the solver.\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.vcanonical","page":"Home","title":"SchrodingerToolkit.vcanonical","text":"v = vcanonical(S, n, idx)\n\nCreate a vector of storage type S of length n where the position idx is one and any other position is zero.\n\n\n\n\n\n","category":"function"},{"location":"#SchrodingerToolkit.vfirst","page":"Home","title":"SchrodingerToolkit.vfirst","text":"v = vfirst(V)\n\nReturn the first element of the vector V.\n\n\n\n\n\n","category":"function"},{"location":"#SchrodingerToolkit.vlast","page":"Home","title":"SchrodingerToolkit.vlast","text":"v = vlast(V)\n\nReturn the last element of the vector V.\n\n\n\n\n\n","category":"function"},{"location":"#SchrodingerToolkit.vones-Tuple{Any, Any}","page":"Home","title":"SchrodingerToolkit.vones","text":"v = vones(S, n)\n\nCreate a vector of storage type S of length n only composed of one.\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.vseq-Tuple{Any, Any}","page":"Home","title":"SchrodingerToolkit.vseq","text":"v = vseq(S, n)\n\nCreate a vector of storage type S of length n where the position i is i.\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.vundef-Tuple{Any, Any}","page":"Home","title":"SchrodingerToolkit.vundef","text":"v = vundef(S, n)\n\nCreate an uninitialized vector of storage type S of length n.\n\n\n\n\n\n","category":"method"},{"location":"#SchrodingerToolkit.vzeros-Tuple{Any, Any}","page":"Home","title":"SchrodingerToolkit.vzeros","text":"v = vzeros(S, n)\n\nCreate a vector of storage type S of length n only composed of zero.\n\n\n\n\n\n","category":"method"},{"location":"Usage/QuickStart/","page":"-","title":"-","text":"To use this project you need to define your problem or model in a understandable way for this package. To do so, you need to follow the next steps:","category":"page"}]
}
