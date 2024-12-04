using LinearAlgebra
using Symbolics

function PsiStr(idx, n = 5,nidx=0)

    np1 = "n+$(nidx+1)"
    np2 = nidx == 0 ? "n" : "n+$(nidx)"
    ldiff = length(np1) - length(np2)
    spaces = " "^ldiff

    uidx = n-idx

    out = "Ψ($idx) "
    for i in 1:(n - uidx)
        out *= "[$(i){$np1}] "
    end

    for i in (n - uidx + 1):n
        out *= "[$(i){$np2}$spaces] "
    end

    return out
end

function ApplyN(idx,n=5,nidx=0)
    np1 = "n+$(nidx+1)"
    np2 = nidx == 0 ? "n" : "n+$(nidx)"
    ldiff = length(np1) - length(np2)
    spaces = " "^ldiff



    out = "F(idx) = "
    out*= "F<$(PsiStr(idx,n,nidx))> - F<$(PsiStr(idx-1,n,nidx))> / Ψ[$idx{$np1}] - Ψ[$idx{$np2}$spaces]" 
end



function FwD(n=5)
    for i in 1:n
        println("   i=$i = > ",ApplyN(i,n))
    end

    
end

function BwD(n=5)
    for i in n:-1:1
        println("   i=$i = > ",ApplyN(i,n,1))
    end
    
end

