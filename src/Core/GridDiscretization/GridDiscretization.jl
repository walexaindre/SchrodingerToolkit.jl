function validate_constraints(a::T, b::T, c::T, α::T, β::T,
                              order::V) where {V<:Integer,
                                               T<:AbstractFloatOrRational{V}}
    if order < 2
        throw(ArgumentError("Order must be greater than 1 and even at space discretization => [ order = $order ]"))
    end

    if mod(order, 2) == 1
        throw(ArgumentError("Order must be even at space discretization => [ order = $order ]"))
    end

    if !isapprox(a + b + c, 1 + 2 * (α + β))
        throw(ArgumentError("Constraint not satisfied => [ a + b + c - 1 - 2 * ( α + β ) = $(a + b + c - 1 - 2 * α - 2 * β) ] != 0"))
    end

    for idx in 2:2:(order - 2)
        pow2 = 2^idx
        pow3 = 3^idx

        fact = (idx + 1) * (idx + 2)

        left_side = a + pow2 * b + pow3 * c
        right_side = fact * (α + pow2 * β)

        if !isapprox(left_side, right_side)
            throw(ArgumentError("Constraint not satisfied => [ a + 2 ^ $idx * b + 3 ^ $idx * c -  ($(idx+2)!/$idx!)( α - 2 ^ $idx * β ) = $(abs(left_side-right_side)) ] !≈ 0"))
        end
    end
end

@inline function validate_constraints(SpaceDiscretization::SecondDerivativeCoefficients{V,T}) where {V<:Integer,
                                                                                             T<:AbstractFloatOrRational{V}}
    return validate_constraints(SpaceDiscretization.a, SpaceDiscretization.b,
                                SpaceDiscretization.c, SpaceDiscretization.α,
                                SpaceDiscretization.β, SpaceDiscretization.order)
end

@inline function validate_positive_definite(α::T,
                                    β::T) where {T<:AbstractFloatOrRational{Int}}

    #By Gershgorin's theorem the generated matrix A is PSD
    #If |λ-1|≤ 2 ( α + β ) => -2 ( α + β ) ≤ λ - 1 ≤ 2 ( α + β ) => 1 - 2 ( α + β ) ≤ λ ≤ 1 + 2 ( α + β )
    #Then 1 - 2 ( α + β ) > 0 is required...
    if 1 - 2 * (α + β) <= 0
        throw(ArgumentError("Generated A matrix will not be positive definite..."))
    end
end

@inline function validate_positive_definite(SpaceDiscretization::SecondDerivativeCoefficients{V,
                                                                                      T}) where {V<:Integer,
                                                                                                 T<:AbstractFloatOrRational{V}}
    validate_positive_definite(SpaceDiscretization.α, SpaceDiscretization.β)
end

@inline function check_validity(SpaceDiscretization::SecondDerivativeCoefficients{V,T}) where {T<:AbstractFloat,
                                                                                       V<:Integer}
    validate_constraints(SpaceDiscretization)
    validate_positive_definite(SpaceDiscretization)
end

@inline function check_validity(a::T, b::T, c::T, α::T, β::T,
                        order::V) where {V<:Integer,
                                         T<:AbstractFloatOrRational{V}}
    validate_constraints(a, b, c, α, β, order)
    validate_positive_definite(α, β)
end

@inline function SecondDerivativeFiniteDifferenceSchemeAssembly(a::T, b::T, c::T, α::T, β::T,
                             order::V) where {V<:Integer,
                                              T<:AbstractFloatOrRational{V}}
    check_validity(a, b, c, α, β, order)
    return SecondDerivativeCoefficients{V,T}(a, b, c, α, β, order)
end

export SecondDerivativeFiniteDifferenceSchemeAssembly, check_validity

include("GridDiscretizationDefaults.jl")
include("GridDiscretizationAssembly.jl")