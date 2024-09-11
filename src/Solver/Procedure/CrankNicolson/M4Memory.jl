struct M4Memory{RealType,ComplexType<:Complex{RealType},
                ComplexVectorType<:AbstractVector{ComplexType},
                ComplexArrayType<:AbstractArray{ComplexType,2},
                LinOp1,
                LinOp2,
                LinOp3} <: AbstractMemory
    current_state::ComplexArrayType
    current_state_abs2::ComplexArrayType
    temp_state_abs2::ComplexVectorType
    stage1::ComplexVectorType
    stage2::ComplexVectorType
    b0_temp::ComplexVectorType
    b_temp::ComplexVectorType

    component_temp::ComplexVectorType

    opA::LinOp1
    factored_opA::LinOp2
    opD::LinOp3
end

@inline current_state(memory::M4Memory) = copy(memory.current_state)
@inline current_state(memory::M4Memory, component_idx) = memory.current_state[:,
                                                                              component_idx]
@inline current_state!(memory::M4Memory) = memory.current_state
@inline current_state!(memory::M4Memory, component_idx) = view(memory.current_state, :,
                                                               component_idx)

function M4Memory(::Type{ComplexVectorType}, ::Type{ComplexArrayType}, PDE, conf, grid,
                  opA,
                  factored_opA,
                  opD) where {ComplexVectorType,
                              ComplexArrayType}
    ncomp = ncomponents(PDE)
    mesh_elems = length(grid)
    vecmem = vzeros(ComplexVectorType, mesh_elems)
    arrmem = vzeros(ComplexArrayType, (mesh_elems, ncomp))

    M4Memory(arrmem, copy(arrmem), vecmem, copy(vecmem), copy(vecmem), copy(vecmem),
             copy(vecmem),
             copy(vecmem),
             opA, factored_opA, opD)
end

@inline function system_mass(M::M4Memory, grid)
    curr_state = current_state!(M)
    measure(grid) * vec(sum(abs2.(curr_state); dims = 1))
end

@inline function system_total_mass(M::M4Memory, grid::AG) where {AG<:AbstractPDEGrid}
    sum(system_mass(M, grid))
end

@inline function system_total_mass(M::M4Memory,
                                    power_per_component::Power) where {Power<:AbstractVector}
    sum(power_per_component)
end

function system_energy(M::M4Memory, PDE, grid)
    factored_opA = M.factored_opA
    D = M.opD

    state_abs2 = M.current_state_abs2
    stage1 = M.stage1
    stage2 = M.stage2
    components = current_state!(M)

    @. state_abs2 = abs2(components)
    F = get_field(PDE)

    energy = zero(eltype(components))

    σ_collection = get_σ(PDE)

    for (idx, σ) in enumerate(σ_collection)
        @views comp = components[:, idx]
        stage1 .= D * comp

        stage2 .= factored_opA \ stage1

        energy -= σ * dot(comp, stage2)
    end

    stage1 .= F(state_abs2)
    vecenergy = sum(stage1)
    energy += vecenergy

    real(energy) * measure(grid)
end