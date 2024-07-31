abstract type AbstractGrid{T<:Integer,N} <: AbstractArray{T,N} end

# An implementation of an abstract N dimensional mesh.
# Here we can map from the N dimensional index to the 1D index and vice versa.

struct PeriodicAbstractGrid{T<:Integer,N} <: AbstractGrid{T,N}
    dims::NTuple{N,T}
    multiplied_dims::NTuple{N,T}
end

struct GhostPeriodicAbstractGrid{T<:Integer,N} <: AbstractGrid{T,N}
    dims::NTuple{N,T}
    multiplied_dims::NTuple{N,T}
    depth::T
end

export AbstractGrid, PeriodicAbstractGrid, GhostPeriodicAbstractGrid