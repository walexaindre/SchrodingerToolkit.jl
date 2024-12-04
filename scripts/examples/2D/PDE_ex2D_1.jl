using SchrodingerToolkit

#iuₜ +Δu+(|u|²+β|v|²)u - f₁= 0
#ivₜ +Δv+(|v|²+β|u|²)v - f₂= 0
#Domain Ω = [0,1]x[0,1]
#Initial conditions: u₀(x,y) = u(x,y,0) = cos(2πx)sin(πx), v₀(x,y) = v(x,y,0) = sin(πx)sin(πy)
#Solution: u(x,y,t) = eⁱᵗ cos(2πx)sin(πx)sin(πy), v(x,y,t) = eⁱᵗ sin(πx)sin(πy)

function F(x)
    
end