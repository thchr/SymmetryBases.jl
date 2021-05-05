using Crystalline
using SymmetryBases
using Test
include("consistency_check.jl")

test_nontopo = true
spinful = false
timereversal = true
algorithm = "DualMode" # DualMode or PrimalMode
for sgnum in 1:MAX_SGNUM[3]
    BRS = bandreps(sgnum, spinful=spinful, timereversal=timereversal)
    
    B = matrix(BRS, true)        # Matrix with columns of EBRs.
    
    F   = Crystalline.smith(B)   # Smith normal decomposition of B
    dᵇˢ = count(!iszero, F.SNF)  # "Dimensionality" of band structure
    Nⁱʳʳ, Nᴱᴮᴿ = size(B)

    # Print some simple stuff early on, to indicate that a calculation is running
    println("\nSG", sgnum, ": ", classification(BRS), 
            " (", dᵇˢ, " \"band structure dimensions\"; ", Nⁱʳʳ, " inequalities)")

    # Compatibility Hilbert basis  
    sb = compatibility_basis(F, BRS, algorithm=algorithm)
    nsᴴ = matrix(sb) 
    Nᴴ = length(sb) # Number of Hilbert basis elements

    # Nontopological Hilbert basis 
    sb_nontopo  = nontopological_basis(F, BRS, algorithm=algorithm)
    nsᴴ_nontopo = matrix(sb_nontopo)
    Nᴴ_nontopo  = length(sb_nontopo)

    # Splitting into trivial and fragile Hilbert basis elements
    trivial_idxs, fragile_idxs = split_fragiletrivial(sb_nontopo, B)

    # Write some stats about the obtained Hilbert basis
    println("   ", Nᴱᴮᴿ,       " EBRs")
    println("   ", Nᴴ,         " compatibility elements")
    println("   ", Nᴴ_nontopo, " nontopological elements")
    println("      ", size(trivial_idxs, 2), " trivial elements")
    println("      ", size(fragile_idxs, 2), " fragile elements")

    # Test consistency of basis, if requested
    if test_nontopo
        _test_hilbert_basis_consistency(BRS, F, nsᴴ, nsᴴ_nontopo)
    end
end