#Temporary workaround before this is patched at CUDA.jl
function SparseArrays.spdiagm(v::CuArray{Tv}) where {Tv}
    nzVal = v
    N = Int32(length(v))

    colPtr = CuArray(one(Int32):(N + one(Int32)))
    rowVal = CuArray(one(Int32):N)
    dims = (N, N)
    CuSparseMatrixCSC(colPtr, rowVal, nzVal, dims)
end

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

function M2Memory(::Type{ComplexVectorType}, ::Type{ComplexArrayType}, PDE, conf, grid,
                  opA,
                  opD; krylov_gmres_memory = 20,
                  energy_solver_params = IterativeLinearSolver(backend_type(conf))) where {ComplexVectorType,
                                                                                           ComplexArrayType}
    ncomp = ncomponents(PDE)
    mesh_elems = length(grid)

    vecmem = vzeros(ComplexVectorType, mesh_elems)
    arrmem = vzeros(ComplexArrayType, (mesh_elems, ncomp))

    M2Memory(arrmem, copy(arrmem), vecmem, copy(vecmem), copy(vecmem), copy(vecmem),
             copy(vecmem),
             copy(vecmem),
             opA, opD, initialize_gmres_kylov_solver(vecmem, krylov_gmres_memory),
             energy_solver_params)
end

@inline function aux_josephson_power(M::M2Memory, τ,
                             Γ)
    curr_state = current_state!(M)

    josephson_energy = zero(eltype(curr_state))

    for comp_i in 1:ncomponents(PDE)
        for comp_j in (comp_i + 1):ncomponents(PDE)
            josephson_energy += dot(view(curr_state, :, comp_i),
                                    view(curr_state, :, comp_j))
        end
    end

    josephson_energy *= (τ * Γ)

    josephson_energy
end

@inline function system_power(M::M2Memory, grid)
    curr_state = current_state!(M)

    measure(grid) * vec(sum(abs2.(curr_state); dims = 1))
end

@inline function system_total_power(M::M2Memory, PDE::PDEq, grid) where {PDEq<:SchrodingerPDE}
    curr_state = current_state!(M)
    aux_josephson_power(M, grid.τ, junction_coefficient(PDE))
    measure(grid) * sum(sum(abs2.(curr_state); dims = 1)) - imag(josephson_energy)
end

@inline function system_total_power(M::M2Memory, PDE::PDEq, grid,
                            power_per_component) where {PDEq<:SchrodingerPDE}
    aux_josephson_power(M, grid.τ, junction_coefficient(PDE))
    sum(power_per_component) - imag(josephson_energy)
end

function aux_sys_energy_josephson(PDE, curr_state)
    Γ = junction_coefficient(PDE)

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

function aux_sys_energy_trapping_potential(PDE, grid, curr_state)
    grid_points = typeof(curr_state)(collect_points(grid))

    partial_res = zero(eltype(curr_state))

    for comp_i in 1:ncomponents(PDE)
        V = trapping_potential(PDE, comp_i)

        partial_res += sum(V(grid_points) .* abs2.(curr_state))
    end
    real(partial_res)
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

    aux_josephson = aux_sys_energy_josephson(PDE, components)
    aux_trapping_potential = aux_sys_energy_trapping_potential(PDE, grid, components)
    println("aux_josephson: ", aux_josephson)
    println("aux_trapping_potential: ", aux_trapping_potential)
    ren = real(energy) * measure(grid) + aux_josephson + aux_trapping_potential
    println("sys_energy: ", ren)
    ren
end