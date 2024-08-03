################################################################################################
##                                                                                         
##     Generic Methods for AbstractGrid                           
##
################################################################################################

@inline Base.size(A::T) where {T<:AbstractPDEGrid} = A.dims
@inline Base.length(A::T) where {T<:AbstractPDEGrid} = prod(A.dims)

@inline range_start_step_stop(start, step, stop) = start:step:stop
@inline range_start_stop_step(start, stop, step) = start:step:stop
@inline range_start_stop_step(start_stop, step) = start_stop[1]:step:start_stop[2]
@inline range_start_stop_length(start, stop, length) = range(start, stop;
                                                             length = length)
@inline range_start_stop_length(start_stop, length) = range(start_stop[1],
                                                            start_stop[2];
                                                            length = length)
@inline range_length(rank) = Base.step(rank)

################################################################################################
##                                                                                         
##     Periodic Grid Methods                           
##
################################################################################################

#TODO Perform integrity checks...

@inline PeriodicGrid(::Type{V}, ::Type{T}, τ::T, ranges::NTuple{N,R}) where {V<:Integer,T<:Real,R<:AbstractRange{T},N} = PeriodicGrid(ranges,
                                                                                                                                      size.(ranges,
                                                                                                                                            1),
                                                                                                                                      range_length.(ranges),
                                                                                                                                      τ)
@inline PeriodicGrid(::Type{V}, ::Type{T}, τ::T, ranges::R...) where {V<:Integer,T<:Real,R<:AbstractRange{T}} = PeriodicGrid(V,
                                                                                                                             T,
                                                                                                                             τ,
                                                                                                                             ranges)

@inline function PeriodicGrid(::Type{V}, ::Type{T}, τ::T, N::NTuple{_N,V},
                              boundaries::NTuple{_N,Tuple{T,T}}) where {V<:Integer,
                                                                        T<:Real,_N}
    ranges = Tuple(range_start_stop_length(boundaries[idx], N[idx]) for idx in 1:_N)
    dims = Tuple(size(ranges[idx], 1) for idx in 1:_N)

    return PeriodicGrid(ranges,
                        dims,
                        range_length.(ranges),
                        τ)
end

@inline function PeriodicGrid(::Type{V}, ::Type{T}, τ::T, h::NTuple{N,T},
                              boundaries::NTuple{N,Tuple{T,T}}) where {V<:Integer,
                                                                       T<:Real,N}
    ranges = Tuple(range_start_stop_step(boundaries[idx], h[idx]) for idx in 1:N)
    dims = Tuple(size(ranges[idx], 1) for idx in 1:N)

    return PeriodicGrid(ranges,
                        dims,
                        h,
                        τ)
end

@inline PeriodicGrid(::Type{ComputeBackend}, PDE::PDEeq, τ::RealType, h::NTuple{N,RealType}) where {N,IntType,RealType,PDEeq<:SchrodingerPDE{N,RealType},ComputeBackend<:AbstractBackend{IntType,RealType}} = PeriodicGrid(IntType,
                                                                                                                                                                                                                                RealType,
                                                                                                                                                                                                                                τ,
                                                                                                                                                                                                                                h,
                                                                                                                                                                                                                                PDE.boundaries)

@inline PeriodicGrid(::Type{ComputeBackend}, PDE::PDEeq, τ::RealType, N::NTuple{_N,IntType}) where {_N,IntType,RealType,PDEeq<:SchrodingerPDE{_N,RealType},ComputeBackend<:AbstractBackend{IntType,RealType}} = PeriodicGrid(IntType,
                                                                                                                                                                                                                                 RealType,
                                                                                                                                                                                                                                 τ,
                                                                                                                                                                                                                                 N,
                                                                                                                                                                                                                                 PDE.boundaries)

@inline PeriodicGrid(::Type{V}, ::Type{T}, τ::T, h::NTuple{N,T}, boundaries::Tuple{T,T}...) where {V<:Integer,T<:Real,N} = PeriodicGrid(V,
                                                                                                                                        T,
                                                                                                                                        τ,
                                                                                                                                        h,
                                                                                                                                        boundaries)

@inline function Base.getindex(A::PeriodicGrid{V,T,R,N},
                               I::Vararg{V,N}) where {V<:Integer,T<:Real,
                                                      R<:AbstractRange{T},N}
    #println("IS: ", I)
    @inbounds getindex.(A.ranges, mod1.(I, A.dims))
end

@inline function Base.getindex(A::PeriodicGrid{V,T,R,N}, ::Colon,
                               col::V) where {V<:Integer,T<:Real,R<:AbstractRange{T},
                                              N}
    @boundscheck begin
        if !(1 <= col <= N)
            throw(BoundsError(A, (:, col)))
        end
    end

    #total_product = prod(size(A))
    left_slice = 1:(col - 1)
    right_slice = (col + 1):N

    #This product is handy because prod([]) = 1 (product of an empty collection is 1) 
    #and a range has dimension 1
    inner_elems = prod(size(A)[left_slice])
    outer_elems = prod(size(A)[right_slice])

    return repeat(A.ranges[col]; inner = inner_elems, outer = outer_elems)
end

function get_range(A::PeriodicGrid{V,T,R,N},
                   dim::V) where {V<:Integer,T<:Real,R<:AbstractRange{T},N}
    A.ranges[dim]
end

function collect_points(A::PeriodicGrid{V,T,R,N}) where {V<:Integer,T<:Real,
                                                         R<:AbstractRange{T},N}
    foldl(hcat, (A[:, idx] for idx in 1:N))
end

@inline Base.convert(::Type{PeriodicAbstractGrid{V,N}}, A::PeriodicGrid{V,T,R,N}) where {V<:Integer,T<:Real,R<:AbstractRange{T},N} = PeriodicAbstractMesh(V,
                                                                                                                                                          A.dims)

"""
    measure(A::PeriodicGrid{V,T,R,N}) where {V<:Integer,T<:Real,R<:AbstractRange{T},N}

Return the measure of the grid.
"""
@inline measure(A::PeriodicGrid{V,T,R,N}) where {V<:Integer,T<:Real,R<:AbstractRange{T},N} = prod(A.h)