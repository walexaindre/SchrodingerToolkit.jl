SchrodingerPDEComponent(σ::Tv,f::Fn,ψ::IC,V::TP = nothing,Γ::Γv = nothing) where {Tv,Fn,IC,TP,Γv} = SchrodingerPDEComponent(σ,f,ψ,V,Γ)

@inline has_trapping_potential(Component::SchrodingerPDEComponent{Tv,Fn,InitialCondition,TrappingPotential}) where {Tv,Fn,InitialCondition,TrappingPotential} = if TrappingPotential == Nothing return false else return true end

@inline has_josephson_junction(Component::SchrodingerPDEComponent{Tv,Γv,Fn,InitialCondition,TrappingPotential}) where {Tv,Γv,Fn,InitialCondition,TrappingPotential} = if Γv == Nothing return false else return true end

@inline has_josephson_junction(PDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = any(has_josephson_junction, PDE.components)
@inline has_josephson_junction(PDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = any(has_josephson_junction, PDE.components)
@inline has_trapping_potential(PDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = any(has_trapping_potential, PDE.components)
@inline has_trapping_potential(PDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = any(has_trapping_potential, PDE.components)



"Number of dimensions for the Schrodinger PDE"
@inline Base.ndims(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = N

"Number of dimensions for the Schrodinger PDE"
@inline Base.ndims(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = N

#"Number of components for the Schrodinger PDE"
#@inline Base.length(SPDE::SchrodingerPDEPolynomic{N,Tv,Comp,Potential,Optimized}) where {N,Tv,Comp,Potential,Optimized} = length(Comp)

#"Number of components for the Schrodinger PDE"
#@inline Base.length(SPDE::SchrodingerPDENonPolynomic{N,Tv,Comp,Potential}) where {N,Tv,Comp,Potential} = length(Comp)



"Type of the elements for the Schrodinger PDE"
@inline Base.eltype(::Type{SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}}) where {N,Tv,MComp,Potential,Optimized} = Tv

"Type of the elements for the Schrodinger PDE"
@inline Base.eltype(::Type{SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}}) where {N,Tv,MComp,Potential} = Tv

"Boundary start and end point at index dimension"
@inline get_boundary(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.boundaries[index]

"Boundary start and end point at index dimension"
@inline get_boundary(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.boundaries[index]

"Obtain the component of the Schrodinger PDE at index"
@inline get_component(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index]

"Obtain the component of the Schrodinger PDE at index"
@inline get_component(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index]

"Dispersion coefficient for the Schrodinger PDE component at index"
@inline get_σ(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].σ

@inline get_σ(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = ntuple(idx -> SPDE.components[idx].σ,
                                                                                                                             length(SPDE.components))

"Dispersion coefficient for the Schrodinger PDE component at index"
@inline get_σ(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].σ

@inline get_σ(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = ntuple(idx -> SPDE.components[idx].σ,
                                                                                                            length(SPDE.components))

"Function for the Schrodinger PDE component"
@inline get_f(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].f

"Function for the Schrodinger PDE component"
@inline get_f(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].f

"Initial condition for the Schrodinger PDE component"
@inline get_ψ(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].ψ

"Initial condition for the Schrodinger PDE component"
@inline get_ψ(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].ψ

"Trapping potential for the i-th Schrodinger PDE component"
@inline get_V(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].V

"Trapping potential for the i-th Schrodinger PDE component"
@inline get_V(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].V

"Josephson Junction Coefficient for the i-th Schrodinger PDE component"
@inline get_Γ(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].Γ

"Josephson Junction Coefficient for the i-th Schrodinger PDE component"
@inline get_Γ(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].Γ

"Potential for the Schrodinger PDE"
@inline get_field(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = SPDE.F

"Potential for the Schrodinger PDE"
@inline get_field(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = SPDE.F

"Simplified polynomial to evade catastrophic cancellation"
@inline get_optimized(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = SPDE.N

"Non polynomial potentials can't be optimized so return nothing"
@inline get_optimized(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = nothing

"Finish time for the Schrodinger PDE"
@inline get_time_boundary(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = SPDE.T

"Finish time for the Schrodinger PDE"
@inline get_time_boundary(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = SPDE.T

"Evaluation of the initial condition for the Schrodinger PDE"
@inline function evaluate_ψ(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,
                                                          Optimized}, P::PGrid,
                            Memory) where {N,Tv,MComp,Potential,Optimized,
                                           PGrid<:PeriodicGrid}
    cstate = Memory.current_state
    points = typeof(cstate)(collect_points(P))
    for i in 1:ncomponents(SPDE)
        ψ_ = view(cstate, :, i)
        func = get_ψ(SPDE, i)
        ψ_ .= func(points)
    end
    nothing
end

"Evaluation of the initial condition for the Schrodinger PDE"
@inline function evaluate_ψ(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential},
                            P::PGrid,
                            Memory) where {N,Tv,MComp,Potential,PGrid<:PeriodicGrid}
    cstate = Memory.current_state
    points = typeof(cstate)(collect_points(P))
    for i in 1:ncomponents(SPDE)
        ψ_ = view(cstate, :, i)
        func = get_ψ(SPDE, i)
        ψ_ .= func(points)
    end
    nothing
end

"Number of components for the Schrodinger PDE"
@inline ncomponents(SPDE::PDEeq) where {PDEeq<:SchrodingerPDE} = length(SPDE.components)

"Estimate the number of timesteps for the Schrodinger PDE"
@inline function estimate_timesteps(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,
                                                                  Potential,
                                                                  Optimized},
                                    P::PGrid) where {N,Tv,MComp,Potential,Optimized,
                                                     PGrid<:PeriodicGrid}
    τ = P.τ
    T = get_time_boundary(SPDE)
    rank = 0:τ:T
    return length(rank) - 1
end

"Estimate the number of timesteps for the Schrodinger PDE"
@inline function estimate_timesteps(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,
                                                                     Potential},
                                    P::PGrid) where {N,Tv,MComp,Potential,
                                                     PGrid<:PeriodicGrid}
    τ = P.τ
    T = get_time_boundary(SPDE)
    rank = 0:τ:T
    return length(rank) - 1
end

export get_boundary, get_component, get_σ, get_f, get_ψ, get_field,
       get_time_boundary, ncomponents, get_optimized, evaluate_ψ, estimate_timesteps