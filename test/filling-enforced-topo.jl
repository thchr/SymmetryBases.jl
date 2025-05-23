using Crystalline, SymmetryBases
using ProgressMeter

# ---------------------------------------------------------------------------------------- #


sgnums       = 1:230 # 47, 83, 123
timereversal = false
check_fragile_split = false
verbose      = true

for sgnum in sgnums
    println("SG ", sgnum)
    # prep-work to get Hilbert bases etc
    brs  = bandreps(sgnum, spinful=false, timereversal=timereversal)
    B    = stack(brs) # matrix with columns of EBRs
    isℤ₁ = indicator_group_as_string(brs) == "Z₁"
    if isℤ₁
        print("── trivial symmetry indicator group") 
    end 
    flush(stdout)

    F  = Crystalline.smith(B) # Smith normal decomposition of B
    sb = compatibility_basis(F, brs)
    verbose && println("── computed sb ($(length(sb)) vectors)"); flush(stdout) 

    # compute non-topological Hilbert basis
    if check_fragile_split
        if !timereversal && sgnum ∈ (83, 174, 175, 176)
            println("── skipped trivial/fragile basis elements split: computational quagmire")
        else
            nontopo_sb = nontopological_basis(F, brs)
            verbose && println("── computed nontopo_sb ($(length(nontopo_sb)) vectors)")
            flush(stdout)

            trivial_idxs, fragile_idxs = split_fragiletrivial(nontopo_sb, B)
            verbose && println("── split fragile and trivial elements ($(length(trivial_idxs))|$(length(fragile_idxs)))"); flush(stdout)
        end
    end

    # check if any fillings are solely topological
    μs = fillings(sb)
    unique_μs = sort(unique(μs))
    for μ in unique_μs
        idxs = findall(==(μ), μs)
        print("   μ = ", μ, 
              " "^(ndigits(maximum(μs))-ndigits(μ)), " (", length(idxs), " vectors)")
        fe_nontrivial = true
        fe_fragile = fe_mixed = true
        for i in idxs
            n = sb[i]
            # check both "Z₂" (trivial/nontrivial) TQC-type topology and fragile topology
            topo = calc_detailed_topology(n, B, F)
            
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