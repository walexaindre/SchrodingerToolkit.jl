using SchrodingerToolkit

const α₁ = 1.0
const α₂ = 1.0
const σ₁₁ = 1.0
const σ₂₂ = 1.0
const σ₁₂ = 2.0 / 3.0
const σ₂₁ = σ₁₂

const T = 40.0

const Ω = (-8.0, 8.0) #Domain

f₁(x) = @. @views σ₁₁ * x[:, 1] + σ₁₂ * x[:, 2]

f₂(x) = @. @views σ₂₁ * x[:, 1] + σ₂₂ * x[:, 2]

F(x) = @. @views 0.5 * (σ₁₁ * x[:, 1]^2 + σ₂₂ * x[:, 2]^2) + σ₁₂ * x[:, 1] * x[:, 2]

ψ₁(x) = @. @views 1 / sqrt(pi) * exp(-(x[:, 1] - 2.0)^2 - (x[:, 2] - 2.0)^2)

ψ₂(x) = @. @views 1 / sqrt(pi) * exp(-(x[:, 1] + 2.0)^2 - (x[:, 2] + 2.0)^2)

#x: 
#-0.5r*σ₁₁ - 0.5x*σ₁₁ - y*σ₁₂
#y: 
#-0.5r*σ₂₂ - x*σ₁₂ - 0.5y*σ₂₂
function N(prev, next, idx)
    if idx == 1
        return @. @views -0.5 * σ₁₁ * next - 0.5 * σ₁₁ * prev[:, 1] - σ₁₂ * prev[:, 2]
    else
        return @. @views -0.5 * σ₂₂ * next - σ₁₂ * prev[:, 1] - 0.5 * σ₂₂ * prev[:, 2]
    end
end

@inline V1(x) = @. @views 0 * (x[:, 1])

@inline V2(x) = V1(x)

Γ = 0

C1 = SchrodingerPDEComponent(α₁, f₁, ψ₁, V1, Γ)
C2 = SchrodingerPDEComponent(α₂, f₂, ψ₂, V2, Γ)

PDE = SchrodingerPDEPolynomic((Ω, Ω), (C1, C2), F, N, T)