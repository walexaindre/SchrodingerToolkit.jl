abstract type SchrodingerPDE{N,RealType} end

"PDE component"
struct SchrodingerPDEComponent{Tv,Fn,Hn}
    σ::Tv #Dispersion coefficient
    f::Fn 
    ψ::Hn #Initial condition
end