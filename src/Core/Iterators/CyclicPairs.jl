# Code from ffevotte => https://discourse.julialang.org/t/cyclic-pair-iterator/23156/4

struct CyclicPairs{T}
    subiter::T
end

@inline function Base.iterate(cp::CyclicPairs)
    i = iterate(cp.subiter)
    i === nothing && return nothing
    first, substate = i

    iterate(cp, (substate, first, first, #=finished=#false))
end

@inline function Base.iterate(cp::CyclicPairs, state)
    (substate, latest, first, finished) = state

    i = iterate(cp.subiter, substate)
    if i === nothing
        if finished
            return nothing
        else
            return ((latest, first), (substate, latest, first, #=finished=#true))
        end
    end
    current, substate = i

    return ((latest, current), (substate, current, first, #=finished=#false))
end

@inline Base.length(cp::CyclicPairs) = length(cp.subiter)
@inline Base.eltype(cp::CyclicPairs) = Tuple{eltype(cp.subiter), eltype(cp.subiter)}