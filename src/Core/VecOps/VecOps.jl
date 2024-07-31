# Here are defined a simple set of functions to handle vectors and matrices. The functions are defined for both CPU and GPU storage.

"""
    S = matrix_to_vector(M)

Return the dense vector storage type `S` related to the dense matrix storage type `M`.
"""
function matrix_to_vector(::Type{M}) where {M<:DenseMatrix}
    T = hasproperty(M, :body) ? M.body : M
    par = T.parameters
    npar = length(par)
    (2 ≤ npar ≤ 3) || error("Type $M is not supported.")
    if npar == 2
        S = T.name.wrapper{par[1],1}
    else
        S = T.name.wrapper{par[1],1,par[3]}
    end
    return S
end

matrix_to_vector(::Type{M}) where {M<:AbstractCuSparseMatrix} = M.types[3]

matrix_to_vector(::Type{M}) where {M<:SparseMatrixCSC} = M.types[5]

#matrix_to_vector(::Type{M}) where {M<:SparseMatrixCSR} = M.types[5]

"""
    v = vundef(S, n)

Create an uninitialized vector of storage type `S` of length `n`.
"""
vundef(S, n) = S(undef, n)

"""
    v = vzeros(S, n)

Create a vector of storage type `S` of length `n` only composed of zero.
"""
vzeros(S, n) = fill!(S(undef, n), zero(eltype(S)))

"""
    v = vones(S, n)

Create a vector of storage type `S` of length `n` only composed of one.
"""
vones(S, n) = fill!(S(undef, n), one(eltype(S)))

"""
    v = vcanonical(S, n, idx)

Create a vector of storage type `S` of length `n` where the position `idx` is one and any other position is zero.
"""
function vcanonical(S, n, idx = 1)
    v = vzeros(S, n)
    v[idx:idx] .= 1
    return v
end

"""
    v = vseq(S, n)

Create a vector of storage type `S` of length `n` where the position `i` is `i`.
"""
vseq(S, n) = cumsum(vones(S, n))

"""
    v = test_if_zero(V,idx)

Return `true` if the element at position `idx` of the vector `V` is zero.
"""
function test_if_zero end

function test_if_zero(Vec::GPUV, idx) where {GPUV<:CuVector}
    CUDA.@allowscalar if Vec[idx] == zero(eltype(Vec))
        return true
    end
    return false
end

function test_if_zero(Vec::GenV, idx) where {GenV<:AbstractVector}
    if Vec[idx] == zero(eltype(Vec))
        return true
    end
    return false
end

"""
    v = vfirst(V)

Return the first element of the vector `V`.
"""
function vfirst end

vfirst(Vec::GenV) where {GenV<:AbstractVector} = first(Vec)

vfirst(Vec::GPUV) where {GPUV<:CuVector} = CUDA.@allowscalar first(Vec)

"""
    v = vlast(V)

Return the last element of the vector `V`.
"""
function vlast end

vlast(Vec::GenV) where {GenV<:AbstractVector} = last(Vec)

vlast(Vec::GPUV) where {GPUV<:CuVector} = CUDA.@allowscalar last(Vec)