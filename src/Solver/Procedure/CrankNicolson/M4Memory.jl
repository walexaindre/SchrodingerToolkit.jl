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