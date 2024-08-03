abstract type SchrodingerPDE{N,RealType} end

"PDE component"
struct SchrodingerPDEComponent{Tv,Γv,Fn,InitialCondition,TrappingPotential} 
    σ::Tv #Dispersion coefficient => Must be non negative
    f::Fn #∂∇F/∂xᵢ(x) = f(x)
    ψ::InitialCondition #Initial condition 
    V::TrappingPotential # External trapping potential
    Γ::Γv #Josephson Junction Coefficient => In some documents this coefficient is the same for all components but I guess it can be different under some circumstances.
    SchrodingerPDEComponent(σ::Tv, f::Fn, ψ::InitialCondition, V::TrappingPotential = nothing, Γ::Γv = nothing) where {Tv,Fn,InitialCondition,TrappingPotential,Γv} = new(σ,f,ψ,V,Γ)
end

"Generic handler for non polynomic potentials"
struct SchrodingerPDENonPolynomic{N,Tv,MComp,Potential} <: SchrodingerPDE{N,Tv}
    boundaries::NTuple{N,NTuple{2,Tv}}
    components::MComp
    F::Potential
    T::Tv
end

"Optimized structure for polynomical potentials"
struct SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized} <: SchrodingerPDE{N,Tv}
    boundaries::NTuple{N,NTuple{2,Tv}}
    components::MComp
    F::Potential
    N::Optimized
    T::Tv
end