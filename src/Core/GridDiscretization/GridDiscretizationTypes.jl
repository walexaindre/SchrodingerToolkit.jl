AbstractFloatOrRational{V} = Union{AbstractFloat,Rational{V}} where {V<:Integer}

abstract type FiniteDifferenceScheme end

abstract type SecondDerivativeFiniteDifferenceScheme <: FiniteDifferenceScheme end

struct SecondDerivativeCoefficients{V <: Integer,T <: AbstractFloatOrRational{V}} <: SecondDerivativeFiniteDifferenceScheme
    a::T
    b::T
    c::T
    α::T
    β::T
    order::V
end

export SecondDerivativeCoefficients, SecondDerivativeFiniteDifferenceScheme