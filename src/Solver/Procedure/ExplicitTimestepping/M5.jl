include("M5Memory.jl")

struct M5{RealType,Grid,TKernel,ItSolver,TTime,StoppingCriterion} <:
       AbstractSolverMethod{RealType}
    grid::Grid
    Kernel::TKernel

    Δτ::RealType

    linear_solve_params::ItSolver
    time_collection::TTime
    stopping_criteria::StoppingCriterion
    assembly_time::RealType
end

@inline function update_component!(method::M5, memory, stats, PDE, τ, σ, component_index)

    #Method to be applied here: U_τ+Δτ = U_t + Δτ f(Uₜ) - Δτ/Δt (Uτ - Uₜ)
    #Where Uₜ = f(U) and τ is fictional time and t is real time.
    Δτ = method.Δτ
    Δt = get_τ(PDE)
    Δtm = Δτ / Δt
    ncheck = method.ncheck

    solved = false
    exit_iterations = 0
    stopping_criteria_m5 = stopping_criteria(method)

    for l in 1:get_max_iterations(stopping_criteria_m5)
        update = l % ncheck == 0 # Check if we should calculate the norm 

        if update
            stage1 .= zl
        end

        UτΔτ = 



        if update
            #Calculate the norm
            znorm = grid_measure * norm(stage1)
            solved = znorm <=
                     get_atol(stopping_criteria_m2) +
                     get_rtol(stopping_criteria_m2) * znorm
            if solved
                exit_iterations = l
                break
            end
        end
    end

    if !solved
        @warn "Convergence not reached in $(get_max_iterations(stopping_criteria_m5)) iterations..."
    end
    exit_iterations
end

function step!(method::M5, memory, stats, PDE, conf::SolverConfig)
    start_timer = time()
    grid = method.grid

    σ_forward = get_σ(PDE)
    σ_backward = reverse(σ_forward)

    for τ in time_collection(method)
        #Forward
        for (component_index, σ) in enumerate(σ_forward)
            steps = update_component!(method, memory, stats, PDE, τ, σ, component_index)
            update_component_update_steps!(stats, steps)
        end

        #Backward
        for (component_index, σ) in zip(length(σ_backward):-1:1, σ_backward)
            steps = update_component!(method, memory, stats, PDE, τ, σ, component_index)
            update_component_update_steps!(stats, steps)
        end
    end
    work_timer = time() - start_timer

    update_stats!(stats, memory, grid, PDE,
                  work_timer)
    work_timer
end

export M5, step!