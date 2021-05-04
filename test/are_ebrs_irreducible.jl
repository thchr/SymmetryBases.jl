# PURPOSE:
# Check whether the EBR basis is minimal - i.e., whether each element is irreducible in the 
# sense that it can not be written as an integer cononical combination of other EBRs. 
# Since minimality is a necessary condition for it a basis to be a Hilbert basis, this is a
# quick way to check whether a given EBR basis is a Hilbert basis.
# In general, we find that the answer is: "EBRs are frequently a Hilbert basis, but not 
# generally".

using Crystalline
using SymmetryBases: has_posint_expansion
using JuMP, GLPK

has_tr = true
for sgnum in 1:230
    println("SG$sgnum")
    brs = bandreps(sgnum, 3, timereversal=has_tr, allpaths=true)
    B = matrix(brs, true)
    Nᴱᴮᴿ = size(B, 2)
    for (idx, b) in enumerate(collect(eachcol(B)))
        idxs = vcat(1:idx-1, idx+1:Nᴱᴮᴿ)
        isempty(idxs) && continue
        B′ = B[:, idxs]
        m  = has_posint_expansion(b, B′)
        if termination_status(m) == MOI.OPTIMAL # ← EBR is not irreducible in EBR basis
            println(idx)
            print("      EBR $idx decomposes via EBRs ")
            cs′ = value.(m[:c])
            cs′idxs = findall(≠(0), cs′)
            mults = cs′[cs′idxs]
            cs′idxs .+= ifelse.(cs′idxs .≥ idx, 1, 0) # correct for removed index

            println(join(cs′idxs, " + "), "   ", all(isone, mults) ? "" : mults)
        end
    end
end