This library allows you to use CPU, GPU and xPU backends for your computations. You need to specify what backend you want at setup time.

## Available backends

* `CPU`: `CPUBackend{IntType,FloatType}`
* `GPU`: `GPUBackend{IntType,FloatType}`
* `xPU`: `xPUBackend{IntType,FloatType}`

where `IntType` and `FloatType` are the integer and floating point types you want to use. For example, `CPUBackend{Int32,Float32}` will use 32-bit integers and 32-bit floating point numbers on the CPU.

## Usage

```julia
    using SchrodingerToolkit

    # Use the CPU backend
    backend = CPUBackend{Int32,Float32} # Note that we are using the type and not the instance.
    ...
    SolverParameters(backend, 2, (:ord4, :ord4), :tord2_1_1) # Simply pass the backend to the SolverParameters constructor.
    ...
```

If for some reason you're using a method that isn't supported by the backend you've chosen, you'll get an error message telling you so.
 