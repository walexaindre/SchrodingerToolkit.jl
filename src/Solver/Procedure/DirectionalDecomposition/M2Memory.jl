struct M2Memory{RealType,ComplexType<:Complex{RealType},
                ComplexVectorType<:AbstractVector{ComplexType},
                ComplexArrayType<:AbstractArray{ComplexType,2},SolverStorage,ItSolver,
                LinOp1,
                LinOp2}
    current_state::ComplexArrayType
    current_state_abs2::ComplexArrayType
    temp_state_abs2::ComplexVectorType
    stage1::ComplexVectorType
    stage2::ComplexVectorType
    b0_temp::ComplexVectorType
    b_temp::ComplexVectorType

    component_temp::ComplexVectorType

    opA::LinOp1
    opD::LinOp2

    solver_mem::SolverStorage
    energy_solver_params::ItSolver
end

@inline current_state(memory::M2Memory) = copy(memory.current_state)
@inline current_state(memory::M2Memory, component_idx) = memory.current_state[:,
                                                                              component_idx]
@inline current_state!(memory::M2Memory) = memory.current_state
@inline current_state!(memory::M2Memory, component_idx) = view(memory.current_state, :,
                                                               component_idx)

@inline solver_memory!(memory::M2Memory) = memory.solver_mem

@inline energy_solver_params(memory::M2Memory) = memory.energy_solver_params

function system_perturbation_of_total_mass(M::M2Memory, PDE, grid)
    curr_state = current_state!(M)
    coef = junction_coefficient(PDE)
    τ = grid.τ

    josephson_energy = zero(eltype(curr_state))

    for comp_i in 1:ncomponents(PDE)
        for comp_j in (comp_i + 1):ncomponents(PDE)
            josephson_energy += dot(view(curr_state, :, comp_i),
                                    view(curr_state, :, comp_j))
        end
    end

    josephson_energy *= (τ * coef)

    measure(grid) * sum(sum(abs2.(curr_state); dims = 1)) - imag(josephson_energy)
end

function system_power(M::M2Memory, grid)
    curr_state = current_state!(M)
    measure(grid) * vec(sum(abs2.(curr_state); dims = 1)) - imag(josephson_energy)
end

function aux_sys_energy_josephson(Γ, curr_state)
    partial_res = zero(eltype(curr_state))
    for comp_i in 1:ncomponents(PDE)
        for comp_j in 1:ncomponents(PDE)
            if comp_i != comp_j
                partial_res += dot(view(curr_state, :, comp_i),
                                   view(curr_state, :, comp_j))
            end
        end
    end

    Γ * real(partial_res)
end

function aux_sys_energy_trapping_potential(PDE,curr_state)
    V = get_trapping_potential(PDE)
    grid_points = collect_points(grid)|> typeof(curr_state)



    for comp_i in 1:ncomponents(PDE)
        trapping_potential(PDE,comp_i)
        sum(V.() .* abs2.(curr_state))

    end

end

function system_energy(M::M2Memory, PDE, grid)
    A = M.opA
    D = M.opD

    state_abs2 = M.current_state_abs2
    stage1 = M.stage1
    components = current_state!(M)
    sol_mem = solver_memory!(M)
    sol_params = energy_solver_params(M)

    @. state_abs2 = abs2(components)
    F = get_field(PDE)

    energy = zero(eltype(components))

    σ_collection = get_σ(PDE)

    for (idx, σ) in enumerate(σ_collection)
        @views comp = components[:, idx]
        b = D * comp

        gmres!(sol_mem,
               A,
               b;
               restart = true,
               atol = get_atol(sol_params),
               rtol = get_rtol(sol_params),
               itmax = get_max_iterations(sol_params))

        energy -= σ * dot(comp, sol_mem.x)
    end

    stage1 .= F(state_abs2)
    vecenergy = Vector(sum(stage1; dims = 1))
    energy += vecenergy[1]

    real(energy) * measure(grid)
end