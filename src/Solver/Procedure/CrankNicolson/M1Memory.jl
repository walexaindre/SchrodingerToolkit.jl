struct M1Memory{RealType,ComplexType<:Complex{RealType},
                ComplexVectorType<:AbstractVector{ComplexType},
                ComplexArrayType<:AbstractArray{ComplexType,2},SolverStorage,ItSolver,
                LinOp1,
                LinOp2} <: AbstractMemory
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

@inline current_state(memory::M1Memory) = copy(memory.current_state)
@inline current_state(memory::M1Memory, component_idx) = memory.current_state[:,
                                                                              component_idx]
@inline current_state!(memory::M1Memory) = memory.current_state
@inline current_state!(memory::M1Memory, component_idx) = view(memory.current_state, :,
                                                               component_idx)

@inline solver_memory!(memory::M1Memory) = memory.solver_mem

@inline energy_solver_params(memory::M1Memory) = memory.energy_solver_params

@inline initialize_gmres_kylov_solver(vcomplex, memory_size = 20) = GmresSolver(length(vcomplex),
                                                                                length(vcomplex),
                                                                                memory_size,
                                                                                typeof(vcomplex))

function M1Memory(::Type{ComplexVectorType}, ::Type{ComplexArrayType}, PDE, conf, grid,
                  opA,
                  opD; krylov_gmres_memory = 20,
                  energy_solver_params = IterativeLinearSolver(backend_type(conf))) where {ComplexVectorType,
                                                                                           ComplexArrayType}
    ncomp = ncomponents(PDE)
    mesh_elems = length(grid)

    vecmem = vzeros(ComplexVectorType, mesh_elems)
    arrmem = vzeros(ComplexArrayType, (mesh_elems, ncomp))

    M1Memory(arrmem, copy(arrmem), vecmem, copy(vecmem), copy(vecmem), copy(vecmem),
             copy(vecmem),
             copy(vecmem),
             opA, opD, initialize_gmres_kylov_solver(vecmem, krylov_gmres_memory),
             energy_solver_params)
end

@inline function system_power(M::M1Memory, grid)
    curr_state = current_state!(M)
    measure(grid) * vec(sum(abs2.(curr_state); dims = 1))
end

@inline function system_total_power(M::M1Memory, grid::AG) where{AG<:AbstractPDEGrid}
    sum(system_power(M, grid))
end

@inline function system_total_power(M::M1Memory, power_per_component::Power) where {Power<:AbstractVector}
    sum(power_per_component)
end

function system_energy(M::M1Memory, PDE, grid)
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
    vecenergy = sum(stage1)
    energy += vecenergy

    real(energy) * measure(grid)
end

export current_state,current_state!