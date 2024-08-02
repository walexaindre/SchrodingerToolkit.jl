
@inline has_trapping_potential(Component::SchrodingerPDEComponent{Tv,Fn,InitialCondition,TrappingPotential}) where {Tv,Fn,InitialCondition,TrappingPotential} = if TrappingPotential == Nothing return false else return true end

@inline has_josephson_junction(Component::SchrodingerPDEComponent{Tv,Γv,Fn,InitialCondition,TrappingPotential}) where {Tv,Γv,Fn,InitialCondition,TrappingPotential} = if Γv == Nothing return false else return true end

@inline has_josephson_junction(PDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = any(has_josephson_junction, PDE.components)
@inline has_josephson_junction(PDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = any(has_josephson_junction, PDE.components)
@inline has_trapping_potential(PDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = any(has_trapping_potential, PDE.components)
@inline has_trapping_potential(PDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = any(has_trapping_potential, PDE.components)