abstract type AbstractBackend{IntType,FloatType,ComplexType,VectorType} end

abstract type AbstractCPUBackend{IntType,FloatType,ComplexType,VectorType} <: AbstractBackend{IntType,FloatType,ComplexType,VectorType} end

abstract type AbstractGPUBackend{IntType,FloatType,ComplexType,VectorType} <: AbstractBackend{IntType,FloatType,ComplexType,VectorType} end

abstract type AbstractxPUBackend{IntType,FloatType,ComplexType,VectorType} <: AbstractBackend{IntType,FloatType,ComplexType,VectorType} end

struct GPUBackend{IntType,FloatType,ComplexType,VectorType} <: AbstractGPUBackend{IntType,FloatType,ComplexType,VectorType}
    device::IntType
end