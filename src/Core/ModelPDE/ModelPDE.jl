function Base.show(io::IO,
                   Component::SchrodingerPDEComponent{Tv,Γv,Fn,
                                                      InitialCondition,
                                                      TrappingPotential}) where {Tv,
                                                                                 Γv,
                                                                                 Fn,
                                                                                 InitialCondition,
                                                                                 TrappingPotential}
    println("SchrodingerPDEComponent:")
    println("    σ: $(Component.σ)")
    println("    f: $(Component.f |> typeof)")
    println("    ψ: $(Component.ψ |> typeof)")
    println("    V: $(Component.V |> typeof)")
    println("    Γ: $(Component.Γ |> typeof)")
end

function Base.show(io::IO,
                   PDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}) where {N,Tv,
                                                                                  MComp,
                                                                                  Potential}
    println("SchrodingerPDENonPolynomial{$N, $Tv}: ")
    println("    Boundaries: $(PDE.boundaries)")
    println("    Components: $(PDE.components)")
    println("    Potential : $(PDE.F |> typeof)")
    println("Final time: $(PDE.T)")
end

function Base.show(io::IO,
                   PDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}) where {N,
                                                                                         Tv,
                                                                                         MComp,
                                                                                         Potential,
                                                                                         Optimized}
    println("SchrodingerPDEPolynomial{$N, $Tv}: ")
    println("    Boundaries: $(PDE.boundaries)")
    println("    Components: $(PDE.components)")
    println("    Potential : $(PDE.F |> typeof)")
    println("    Optimized : $(PDE.N |> typeof)")
    println("Final time: $(PDE.T)")
end

@inline has_trapping_potential(Component::SchrodingerPDEComponent{Tv,Γv,Fn,
                                                                  InitialCondition,
                                                                  TrappingPotential}) where {Tv,
                                                                                             Γv,
                                                                                             Fn,
                                                                                             InitialCondition,
                                                                                             TrappingPotential} =
    if TrappingPotential == Nothing
        return false
    else
        return true
    end

@inline has_josephson_junction(Component::SchrodingerPDEComponent{Tv,Γv,Fn,
                                                                  InitialCondition,
                                                                  TrappingPotential}) where {Tv,
                                                                                             Γv,
                                                                                             Fn,
                                                                                             InitialCondition,
                                                                                             TrappingPotential} =
    if Γv == Nothing
        return false
    else
        return true
    end

@inline has_josephson_junction(PDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = any(has_josephson_junction,
                                                                                                                          PDE.components)
@inline has_josephson_junction(PDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = any(has_josephson_junction,
                                                                                                                                           PDE.components)
@inline has_trapping_potential(PDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = any(has_trapping_potential,
                                                                                                                          PDE.components)
@inline has_trapping_potential(PDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = any(has_trapping_potential,
                                                                                                                                           PDE.components)

"Number of dimensions for the Schrodinger PDE"
@inline Base.ndims(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = N

"Number of dimensions for the Schrodinger PDE"
@inline Base.ndims(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = N

#"Number of components for the Schrodinger PDE"
#@inline Base.length(SPDE::SchrodingerPDEPolynomic{N,Tv,Comp,Potential,Optimized}) where {N,Tv,Comp,Potential,Optimized} = length(Comp)

#"Number of components for the Schrodinger PDE"
#@inline Base.length(SPDE::SchrodingerPDENonPolynomic{N,Tv,Comp,Potential}) where {N,Tv,Comp,Potential} = length(Comp)

"Type of the elements for the Schrodinger PDE"
@inline Base.eltype(::Type{SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}}) where {N,Tv,MComp,Potential,Optimized} = Tv

"Type of the elements for the Schrodinger PDE"
@inline Base.eltype(::Type{SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}}) where {N,Tv,MComp,Potential} = Tv

"Boundary start and end point at index dimension"
@inline get_boundary(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.boundaries[index]

"Boundary start and end point at index dimension"
@inline get_boundary(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.boundaries[index]

"Obtain the component of the Schrodinger PDE at index"
@inline get_component(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index]

"Obtain the component of the Schrodinger PDE at index"
@inline get_component(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index]

"Obtain all the components of the Schrodinger PDE"
@inline get_components(SPDE::PDEeq) where {PDEeq<:SchrodingerPDE} = SPDE.components

"Dispersion coefficient for the Schrodinger PDE component at index"
@inline get_σ(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].σ

@inline get_σ(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = ntuple(idx -> SPDE.components[idx].σ,
                                                                                                                              length(SPDE.components))

"Dispersion coefficient for the Schrodinger PDE component at index"
@inline get_σ(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].σ

@inline get_σ(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = ntuple(idx -> SPDE.components[idx].σ,
                                                                                                             length(SPDE.components))

"Function for the Schrodinger PDE component"
@inline get_f(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].f

"Function for the Schrodinger PDE component"
@inline get_f(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].f

"Initial condition for the Schrodinger PDE component"
@inline get_ψ(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].ψ

"Initial condition for the Schrodinger PDE component"
@inline get_ψ(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].ψ

"Trapping potential for the i-th Schrodinger PDE component"
@inline get_V(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].V

"Trapping potential for the i-th Schrodinger PDE component"
@inline get_V(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].V

"Josephson Junction Coefficient for the i-th Schrodinger PDE component"
@inline get_Γ(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].Γ

"Josephson Junction Coefficient for the i-th Schrodinger PDE component"
@inline get_Γ(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].Γ

"Potential for the Schrodinger PDE"
@inline get_field(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = SPDE.F

"Potential for the Schrodinger PDE"
@inline get_field(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = SPDE.F

"Simplified polynomial to evade catastrophic cancellation"
@inline get_optimized(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = SPDE.N

"Non polynomial potentials can't be optimized so return nothing"
@inline get_optimized(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = nothing

"Finish time for the Schrodinger PDE"
@inline get_time_boundary(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = SPDE.T

"Finish time for the Schrodinger PDE"
@inline get_time_boundary(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = SPDE.T

"Evaluation of the initial condition for the Schrodinger PDE"
@inline function evaluate_ψ(SPDE::SchrodingerPDEPolynomial{N,Tv,MComp,Potential,
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
@inline function evaluate_ψ(SPDE::SchrodingerPDENonPolynomial{N,Tv,MComp,Potential},
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

"Evaluation of the initial condition for the Schrodinger PDE where the output is stored in `Container`"
@inline function evaluate_ψ!(PDE::SPDE, P::PGrid,
                             Container) where {SPDE<:SchrodingerPDE,
                                               PGrid<:PeriodicGrid}
    points = typeof(Container)(collect_points(P))
    for i in 1:ncomponents(PDE)
        ψ_ = view(Container, :, i)
        func = get_ψ(PDE, i)
        ψ_ .= func(points)
    end
    nothing
end

"Number of components for the Schrodinger PDE"
@inline ncomponents(SPDE::PDEeq) where {PDEeq<:SchrodingerPDE} = length(SPDE.components)

"Enumerate each different V for the Schrodinger PDE. Here equalities are tested by reference"
@inline function enumerate_V(SPDE::PDEeq) where {PDEeq<:SchrodingerPDE}
    dict_keys = collect(1:ncomponents(SPDE))
    marked = zeros(Bool, ncomponents(SPDE))

    dict_values = Vector{Int}(undef, ncomponents(SPDE))

    for (idx, comp) in enumerate(SPDE.components)
        if marked[idx]
            continue
        elseif !has_trapping_potential(comp)
            dict_values[idx] = -1
            marked[idx] = true
        end

        for (idx2, comp2) in zip((idx + 1):ncomponents(SPDE), SPDE.components)
            if marked[idx2] || !has_trapping_potential(comp2)
                continue
            else
                if comp.V == comp2.V
                    dict_values[idx2] = idx
                    marked[idx2] = true
                end
            end
        end
    end

    all_marked = all(marked)
    if !all_marked
        throw(ArgumentError("Some trapping potentials are not enumerated"))
    end

    todelete = findall(x->x == -1, Iterators.reverse(dict_values)) 

    deleteat!(dict_keys, todelete)
    deleteat!(dict_values, todelete)

    return Dict(dict_keys, dict_values)
end

"""
    junction!(SPDE::PDEeq, memory, output)

Check if the Schrodinger PDE has a Josephson Junction and return the sum of the current state if it has one.
"""
@inline function junction!(SPDE::PDEeq, memory, output,
                           index) where {PDEeq<:SchrodingerPDE}
    curr_state = current_state!(memory)

    if has_josephson_junction(SPDE)
        if ncomponents(SPDE) == 2
            output .= curr_state[:, 3 - index]
        else
            valid_rank = collect(1:ncomponents(SPDE)) #Get the valid ranks
            deleteat!(valid_rank, index) #Delete the current index
            rank_tosum = view(curr_state, :, valid_rank) #Get a view of the current state without the current index
            output .= sum(rank_tosum; dims = 2) #Reduce
        end

    else
        throw(ArgumentError("Unreachable reached: The Schrodinger PDE doesn't have a Josephson junction"))
    end
end

@inline junction_coefficient(SPDE::PDEeq) where {PDEeq<:SchrodingerPDE} = get_Γ(SPDE, 1)

@inline trapping_potential(SPDE::PDEeq, index) where {PDEeq<:SchrodingerPDE} = get_V(SPDE,
                                                                                     index)

"Estimate the number of timesteps for the Schrodinger PDE"
@inline function estimate_timesteps(PDE::SPDE,
                                    P::PGrid) where {SPDE<:SchrodingerPDE,
                                                     PGrid<:PeriodicGrid}
    τ = P.τ
    T = get_time_boundary(PDE)
    rank = 0:τ:T
    return length(rank) - 1
end

export get_boundary, get_component, get_σ, get_f, get_ψ, get_field,
       get_time_boundary, ncomponents, get_optimized, evaluate_ψ, estimate_timesteps