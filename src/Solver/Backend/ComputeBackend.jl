#The daunting part here is what a backend should be. My thoughts about this is that a backend esentially provide some types
#and functions that are used by some other part of the code. I can't armonize this easily and I think that the best way to do this
#is esentially defining a backend in general like CPU, GPU, etc. The types and functions that are used by the backend must be defined
#by the caller of the backend.

#The sane expectations to me are simply ask for integer and real/floating point types. The method defined and the related structures
#must use these types to create the necessary structures. In my mind, reusing code for backend was a good idea, but I think that
#this is not the case in general.

"Abstract type for our available backends"
abstract type AbstractBackend{RealType,IntegerType} end

"This represents a CPU based backend"
struct CPUBackend{RealType,IntegerType} <: AbstractBackend{RealType,IntegerType}
end

"This represents a GPU based backend"
struct GPUBackend{RealType,IntegerType} <: AbstractBackend{RealType,IntegerType}
end

"This represents a heterogeneous based backend (Can use a mixture of CPU and GPU)"
struct xPUBackend{RealType,IntegerType} <: AbstractBackend{RealType,IntegerType}
end

"Function to check if a backend is a CPU based backend"
function iscpu end

"Function to check if a backend is a GPU based backend"
function isgpu end

"Function to check if a backend is a heterogeneous based backend"
function isxpu end

@inline iscpu(backend::CPUBackend) = true
@inline iscpu(backend::GPUBackend) = false
@inline iscpu(backend::xPUBackend) = false


@inline isgpu(backend::GPUBackend) = true
@inline isgpu(backend::CPUBackend) = false
@inline isgpu(backend::xPUBackend) = false

@inline isxpu(backend::xPUBackend) = true
@inline isxpu(backend::CPUBackend) = false
@inline isxpu(backend::GPUBackend) = false

Base.show(io::IO, backend::CPUBackend{RealType,IntegerType}) = print(io, "CPUBackend: {$RealType, $IntegerType}")
Base.show(io::IO, backend::GPUBackend{RealType,IntegerType}) = print(io, "GPUBackend: {$RealType, $IntegerType}")
Base.show(io::IO, backend::xPUBackend{RealType,IntegerType}) = print(io, "xPUBackend: {$RealType, $IntegerType}")