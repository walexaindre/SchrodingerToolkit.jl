#The daunting part here is what a backend should be. My thoughts about this is that a backend esentially provide some types
#and functions that are used by some other part of the code. I can't armonize this easily and I think that the best way to do this
#is esentially defining a backend in general like CPU, GPU, etc. The types and functions that are used by the backend must be defined
#by the caller of the backend.

#The sane expectations to me are simply ask for integer and real/floating point types. The method defined and the related structures
#must use these types to create the necessary structures. In my mind, reusing code for backend was a good idea, but I think that
#this is not the case in general.

"Abstract type for our available backends"
abstract type AbstractBackend{IntType,RealType} end

"This represents a CPU based backend"
struct CPUBackend{IntType,RealType} <: AbstractBackend{IntType,RealType}
end

"This represents a GPU based backend"
struct GPUBackend{IntType,RealType} <: AbstractBackend{IntType,RealType}
end

"This represents a heterogeneous based backend (Can use a mixture of CPU and GPU)"
struct xPUBackend{IntType,RealType} <: AbstractBackend{IntType,RealType}
end

"Function to check if a backend is a CPU based backend"
function iscpu end

"Function to check if a backend is a GPU based backend"
function isgpu end

"Function to check if a backend is a heterogeneous based backend"
function isxpu end

@inline BackendReal(::Type{BackendType}) where {IntType,RealType,BackendType<:AbstractBackend{IntType,RealType}} = RealType
@inline BackendInt(::Type{BackendType}) where {IntType,RealType,BackendType<:AbstractBackend{IntType,RealType}} = IntType

@inline iscpu(::Type{CPUBackend{I,F}}) where {I,F} = true
@inline iscpu(::Type{GPUBackend{I,F}}) where {I,F} = false
@inline iscpu(::Type{xPUBackend{I,F}}) where {I,F} = false

@inline isgpu(::Type{GPUBackend{I,F}}) where {I,F} = true
@inline isgpu(::Type{CPUBackend{I,F}}) where {I,F} = false
@inline isgpu(::Type{xPUBackend{I,F}}) where {I,F} = false

@inline isxpu(::Type{xPUBackend{I,F}}) where {I,F} = true
@inline isxpu(::Type{CPUBackend{I,F}}) where {I,F} = false
@inline isxpu(::Type{GPUBackend{I,F}}) where {I,F} = false

Base.show(io::IO, backend::GPUBackend{IntType,RealType}) where {IntType,RealType} = print(io,
                                                                                          "GPUBackend: {$IntType, $RealType}")
Base.show(io::IO, backend::CPUBackend{IntType,RealType}) where {IntType,RealType} = print(io,
                                                                                          "CPUBackend: {$IntType, $RealType}")
Base.show(io::IO, backend::xPUBackend{IntType,RealType}) where {IntType,RealType} = print(io,
                                                                                          "xPUBackend: {$IntType, $RealType}")

export CPUBackend, GPUBackend, xPUBackend, iscpu, isgpu, isxpu