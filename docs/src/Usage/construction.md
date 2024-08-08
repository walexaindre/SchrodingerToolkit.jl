# How to define a problem

To construct your system, you simply need to define the following functions and structures:

## `SchrodingerPDEComponent`: Define independently all the components of the Schrodinger equation.
---
```@docs; canonical=false
SchrodingerPDEComponent
```

## `SchrodingerPDE`: Define the Schrodinger equation as a whole. 

### If your potential is polynomial then you need to use `SchrodingerPDEPolynomial`.

```@docs; canonical=false
SchrodingerPDEPolynomial
```

### If your potential is not polynomial then you need to use `SchrodingerPDENonPolynomial`.

```@docs; canonical=false
SchrodingerPDENonPolynomial
```

!!! note
    You must pass the set of components as a tuple to the constructor.

!!! warning "Important"
    You must define your functions in vectorized form where your components or variables are represented by columns. This means that your functions must be able to accept and return arrays.

## Example

```@example

using SchrodingerToolkit

const α₁ = 1.0 
const α₂ = 1.0 
const σ₁₁ = 1.0 
const σ₂₂ = 1.0     
const σ₁₂ = 2.0/3.0
const σ₂₁ = σ₁₂

const T = 40.0

const Ω = (-8.0,8.0) #Domain at each dimension

function f₁(x) # DF/Dx₁
    @. @views -σ₁₁*x[:,1] - σ₁₂*x[:,2]
end

function f₂(x) # DF/Dx₂ 
    @. @views -σ₂₁*x[:,1] - σ₂₂*x[:,2]
end

function F(x) # F
    @. @views -0.5*(σ₁₁*x[:,1]^2 + σ₂₂*x[:,2]^2) - σ₁₂*x[:,1]*x[:,2]
end

function ψ₁(x) # ψ₁ initial condition for component 1
    @. @views 1/sqrt(pi)*exp(-(x[:,1]-2.0)^2-(x[:,2]-2.0)^2)
end

function ψ₂(x) # ψ₂ initial condition for component 2
    @. @views 1/sqrt(pi)*exp(-(x[:,1]+2.0)^2-(x[:,2]+2.0)^2)
end

# Polynomial potential helper function to elude catastrophic cancellation.
#x: 
#-0.5r*σ₁₁ - 0.5x*σ₁₁ - y*σ₁₂
#y: 
#-0.5r*σ₂₂ - x*σ₁₂ - 0.5y*σ₂₂
function N(prev,next,idx)
    if idx == 1
        return @. @views -0.5*σ₁₁*next-0.5*σ₁₁*prev[:,1]-σ₁₂*prev[:,2]
    else
        return @. @views -0.5*σ₂₂*next-σ₁₂*prev[:,1]-0.5*σ₂₂*prev[:,2]
    end
end

# Define the components

C1 = SchrodingerPDEComponent(f₁, ψ₁, α₁)
C2 = SchrodingerPDEComponent(f₂, ψ₂, α₂)

# Define the Schrodinger equation

PDE = SchrodingerPDEPolynomial((Ω,Ω),(C1,C2), F, N, T)
```