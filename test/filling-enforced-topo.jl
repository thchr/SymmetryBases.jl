using Crystalline, SymmetryBases
using Nemo: MatrixSpace, ZZ
using ProgressMeter

# ---------------------------------------------------------------------------------------- #
function calc_topo(n, topo_check, Bℤ, nontopo_M=nothing, trivial_M=nothing)

    if (topo_check == NONTRIVIAL || # ⇐ check only strict trivial/nontrivial distinction
        trivial_M === nothing) # ⇐ no fragile phases possible: use fast linear algebra instead of (slow) optimization
        # check "Z₂" TQC-type topology of Hilbert basis
        return calc_topology(n, Bℤ)

    elseif topo_check == FRAGILE
        # check both "Z₂" TQC-type topology and fragile topology of Hilbert basis
        return calc_detailed_topology(n, nontopo_M, trivial_M)
    end
end
# ---------------------------------------------------------------------------------------- #


sgnums       = 1:230 # 47, 83, 123
timereversal = false
topo_check   = NONTRIVIAL
verbose      = true
workhard     = false

for sgnum in sgnums
    print("SG ", sgnum)
    # prep-work to get Hilbert bases etc
    BRS  = bandreps(sgnum, spinful=false, timereversal=timereversal)
    B    = Crystalline.matrix(BRS, true) # Matrix with columns of EBRs
    isℤ₁ = classification(BRS) == "Z₁"
    if topo_check == NONTRIVIAL && isℤ₁
        println("\n── trivial symmetry indicator group ⇒ skipping\n") 
        continue
    end
    if (!workhard && topo_check == FRAGILE) && ((timereversal && sgnum ∈ (47, 83, 123)) ||
                                                (!timereversal && sgnum ∈ (83, 174, 175, 176)))
        println("\n── skipping computational quagmire\n")
        continue
    end
    println(); flush(stdout)

    F  = Crystalline.smith(B) # Smith normal decomposition of B
    sb = compatibility_basis(F, BRS)
    Bℤ = MatrixSpace(ZZ, size(B)...)(B)
    verbose && println("── computed sb ($(length(sb)) vectors)"); flush(stdout) 

    # compute non-topological Hilbert basis
    if topo_check == FRAGILE
        nontopo_sb = nontopological_basis(F, BRS)
        verbose && println("── computed nontopo_sb ($(length(nontopo_sb)) vectors)")
        flush(stdout)
        nontopo_M = Crystalline.matrix(nontopo_sb)

        trivial_idxs, fragile_idxs = split_fragiletrivial(nontopo_sb, B)
        verbose && println("── split fragile and trivial elements ($(length(trivial_idxs))|$(length(fragile_idxs)))"); flush(stdout)
        can_be_fragile = !isempty(fragile_idxs)
        trivial_M = can_be_fragile ? (@view nontopo_M[:, trivial_idxs]) : nothing
    else
        nontopo_M = trivial_M = nothing
    end

    # check if any fillings are solely topological
    μs = fillings(sb)
    unique_μs = sort(unique(μs))
    for μ in unique_μs
        idxs = findall(==(μ), μs)
        print("   μ = ", μ, 
              " "^(ndigits(maximum(μs))-ndigits(μ)), " (", length(idxs), " vectors)")
        fe_nontrivial = true
        fe_fragile = fe_mixed = (topo_check == FRAGILE)
        for i in idxs
            n = sb[i]
            if topo_check == FRAGILE && n ∈ (@view nontopo_sb[fragile_idxs])
                topo = FRAGILE # fast-path check, not requiring optimization
            else
                topo = calc_topo(n, topo_check, Bℤ, nontopo_M, trivial_M)
            end
            
            if topo == TRIVIAL
                fe_fragile = fe_mixed = fe_nontrivial = false
                break
            end
            topo == FRAGILE    && (fe_nontrivial = false)
            topo == NONTRIVIAL && (fe_fragile = false)
        end

        space = " "^(ndigits(length(sb)) - (ndigits(length(idxs))))
        if fe_nontrivial
            println(":", space, " nontrivial filling-enforced topology")
        elseif fe_fragile
            println(":", space, " fragile filling-enforced topology")
        elseif fe_mixed
            println(":", space, " mixed filling-enforced topology")
        else
            println()
        end
    end
    println()
end