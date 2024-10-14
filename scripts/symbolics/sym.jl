using Symbolics
using Symbolics: variables
using LinearAlgebra

ar = variables(:x, 1:4, 1:4, 1:4)

repr2d = permutedims(reshape(ar, 4, :), (2, 1))
stencil = [2 -1 0 -1; -1 2 -1 0; 0 -1 2 -1; -1 0 -1 2]
expanded_stencil = kron(I(4), stencil)

repr2d_off = repr2d * stencil

repr2d_off_front = expanded_stencil * repr2d