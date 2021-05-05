# PURPOSE: compare computation of non-topological Hilbert basis (i.e. {AI+F}) via Smith
# normal decomposition or directly via EBR basis matrix; check equivalence and compare speed
# OUTCOME: directly using the EBR basis is actually faster in aggregate

using Test
using SymmetryBases, Crystalline

has_tr = true

t_sum = t′_sum = 0.0
for sgnum in 1:230
    !has_tr && sgnum ∈ (83, 174, 175, 176) && continue # computationally unfeasible
    print("SG ", sgnum)
    brs = bandreps(sgnum, 3, timereversal=has_tr, spinful=true)
    A = matrix(brs, true)

    # compute Hilbert basis for {AI+F} with conditions imposed either directly via EBR basis
    # or via a Smith normal decomposition
    t′ = @elapsed begin                             # directly via EBR basis
        C′ = SymmetryBases.PyNormaliz.Cone(inequalities = A)
        C′.Compute("HilbertBasis", "DualMode")
        c′ = C′.HilbertBasis()'
        nontopo_sb′ = collect(eachcol(A*c′))
    end
    print(" (", round(t′, digits=1), " s vs. ")

    t = @elapsed begin                              # via Smith normal form
        F = smith(A)
        dᵇˢ = count(!iszero, F.SNF)
        S = @view F.S[:,1:dᵇˢ]
        Λ = @view F.SNF[1:dᵇˢ]
        SΛ = S .* Λ'
        
        C = SymmetryBases.PyNormaliz.Cone(inequalities = SΛ)
        C.Compute("HilbertBasis", "DualMode")
        c = C.HilbertBasis()'

        nontopo_sb = collect(eachcol(SΛ*c))
    end
    print(round(t, digits=1), " s)")
    # only care about a possible speedup if calculation took some meaningful amount of time
    if max(t′, t) > 0.1
        print(": ", round(Int, (t′/t - 1)*100), "% speedup")
    end
    println()

    if sort(nontopo_sb) ≠ sort(nontopo_sb′)
        error("Approaches did not agree!")
    end

    global t′_sum += t′
    global t_sum  += t
end

# Print aggregate times
println("\nAggregate times:")
println("   w/o Smith: ", round(t′_sum, digits=1), " s")
println("   w/ Smith:  ", round(t_sum, digits=1), " s")