"""
    SchrodingerPDE

    This is an abstract type that represents a Schrodinger PDE. This type is used to define a Schrodinger PDE structure.

    Parameters:

    - `N::Int`: Number of dimensions
    - `RealType`: Type for Real numbers. This is used to define the internal algorithms and data structures.

    Is important to note, that the derived structures must handle boundary conditions, initial conditions and important parts of the PDE.

"""
abstract type SchrodingerPDE{N,RealType} end

"""
    SchrodingerPDEComponent

    This structure represents a component of the Schrodinger PDE. This is a simple structure that contains the necessary information
    to solve the Schrodinger PDE for a single component.

    Fields:

    - `σ::Tv`: Dispersion coefficient => Must be non negative
    - `f::Fn`: ∂∇F/∂xᵢ(x) = f(x)
    - `ψ::InitialCondition`: Initial condition 
    - `V::TrappingPotential`: External trapping potential
    - `Γ::Γv`: Josephson Junction Coefficient => In some documents this coefficient is the same for all components but I guess it can be different under some circumstances.

    By default `V` and `Γ` are set to `nothing` but they can be set to a trapping potential and a Josephson junction coefficient respectively if you want to work with a Gross Pitaevsky PDE.

"""
struct SchrodingerPDEComponent{Tv,Γv,Fn,InitialCondition,TrappingPotential} 
    σ::Tv
    f::Fn
    ψ::InitialCondition
    V::TrappingPotential
    Γ::Γv
    SchrodingerPDEComponent(σ::Tv, f::Fn, ψ::InitialCondition, V::TrappingPotential = nothing, Γ::Γv = nothing) where {Tv,Fn,InitialCondition,TrappingPotential,Γv} = new{Tv,Γv,Fn,InitialCondition,TrappingPotential}(σ,f,ψ,V,Γ)
end

"""
    SchrodingerPDENonPolynomic

    This structure represents a Schrodinger PDE with non-polynomical potentials.

    Fields:

    - `boundaries::NTuple{N,NTuple{2,Tv}}`: Boundaries of the cartesian domain where you want to solve the PDE.
    - `components::MComp`: Components of the Schrodinger PDE
    - `F::Potential`: Non-polynomical potential
    - `T::Tv`: Final time
"""
struct SchrodingerPDENonPolynomic{N,Tv,MComp,Potential} <: SchrodingerPDE{N,Tv}
    boundaries::NTuple{N,NTuple{2,Tv}}
    components::MComp
    F::Potential
    T::Tv
end

"""
    SchrodingerPDEPolynomic

    This structure represents a Schrodinger PDE with polynomical potentials.

    Fields:

    - `boundaries::NTuple{N,NTuple{2,Tv}}`: Boundaries of the cartesian domain where you want to solve the PDE.
    - `components::MComp`: Components of the Schrodinger PDE
    - `F::Potential`: Polynomical potential
    - `N::Optimized`: Optimized structure for polynomical potentials. Here is stored an auxiliary function without the catastrophic cancellation problem.
    - `T::Tv`: Final time
"""
struct SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized} <: SchrodingerPDE{N,Tv}
    boundaries::NTuple{N,NTuple{2,Tv}}
    components::MComp
    F::Potential
    N::Optimized
    T::Tv
end

export SchrodingerPDE, SchrodingerPDEComponent, SchrodingerPDEPolynomic, SchrodingerPDENonPolynomic