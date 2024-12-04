using SchrodingerToolkit

#iuₜ+βΔu+(α₁|u|²+(α₁+2α₂)|v|²)u +γu+Γv= 0
#ivₜ+βΔv+(α₁|v|²+(α₁+2α₂)|u|²)v +γv+Γu= 0

#u(x,y,z,0)=cos(x+y+z), v(x,y,z,0)=cos(x-y-z)
#α₁=α₂= β=γ=1, Γ=1/5 or Γ=10

const α₁ = 1.0
const α₂ = 1.0
const β = 1.0
const γ = 1.0
const Γ = 1/5

const T=1.0

const Ω = (0,2*π) #Domain

function f₁(x)
    @. @views α₁*x[:,1] + 2*α₂*x[:,2]
end

function f₂(x)
    @. @views α₁*x[:,2] + 2*α₂*x[:,1]
end

function V₁(x)
    @. @views γ*x[:,1]
end

function V₂(x)
    @. @views γ*x[:,2]
end

#0.5*α₁(x²+y²)+x*y(α₁+2α₂)
function F(x)
    @. @views 0.5*α₁*(x[:,1]^2+x[:,2]^2)+x[:,1]*x[:,2]*(α₁+2*α₂)
end

function ψ₁(x)
    @. @views cos(x[:,1]+x[:,2]+x[:,3])
end

function ψ₂(x)
    @. @views cos(x[:,1]-x[:,2]-x[:,3])
end

C1 = SchrodingerPDEComponent(β, f₁, ψ₁, V₁, Γ)
C2 = SchrodingerPDEComponent(β, f₂, ψ₂, V₂, Γ)

PDE = SchrodingerPDEPolynomial((Ω, Ω, Ω), (C1, C2), F, T)