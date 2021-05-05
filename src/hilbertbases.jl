
# All band structures can be written as 𝐧 = B𝐩 with pᵢ∈ℚ and nᵢ∈𝐍, and B a matrix whose 
# columns are EBRs. We can decompose B to a Smith normal form B = SΛT, such that all 
# allowable band structures can be written as 𝐧 = S𝐳. Here, S is an integer matrix with 
# elements Sᵢⱼ∈ℤ. To make nᵢ integer, we thus require zᵢ∈ℤ.

"""
    compatibility_bases(F::Smith, BRS::BandRepSet; algorithm, verbose)
    compatibility_bases(sgnum::Integer, D::Integer=3; kwargs...)

Computes the Hilbert basis associated with a Smith normal form `F` of the EBR matrix or from
a space group number `sgnum`, which respects all compatibility relations, returning a
`SymBasis` structure. The returned basis is a non-negative integer coefficient basis
for all possible band structures {BS}.

If the method is called with `sgnum::Integer`, the underlying `BandRepSet` is also returned.

## Keyword arguments

- `algorithm::String`: controls the algorithm used by Normaliz to compute the Hilbert
basis. Choices are `"DualMode"` (default) and `"PrimalMode"`
- `spinful::Bool`: Use single- (`false`, default) or double-valued (`true`) irreps.
- `timereversal::Bool`: Assume presence (`true`, default) or absence (`false`) of
time-reversal symmetry.
- `verbose::Bool`: whether to print progress info during the Normaliz computation
(`false`, default).
"""
function compatibility_bases(F::Smith, BRS::BandRepSet; 
                             algorithm::String="DualMode", verbose::Bool=false)
    # To restrict nᵢ to only positive integers, i.e. ℕ, the values of zᵢ must be such that 
    # ∑ⱼ Sᵢⱼzⱼ ≥ 0. This defines a set of inequalities, which in turn defines a polyhedral
    # integer cone. This is where (Py)Normaliz comes in.
    dᵇˢ = count(!iszero, F.SNF)           # "Dimensionality" of band structure
    S = @view F.S[:,OneTo(dᵇˢ)]           # All the nontrivial conditions on zⱼ

    C = PyNormaliz.Cone(inequalities = S) # Construct cone consistent with Sᵢⱼzⱼ ≥ 0
    verbose && C.setVerbose(true)         # Whether to print progress info
    C.Compute("HilbertBasis", algorithm)  # Compute¹ the Hilbert basis
    zsᴴ  = transpose(C.HilbertBasis())    # Columns are Hilbert basis vectors in 𝐳-space

    nsᴴ  = S*zsᴴ                          # Columns are Hilbert basis vectors in 𝐧-space

    return SymBasis(nsᴴ, BRS, true)       # Bases of all valid symmetry vectors in 𝐧-space
end

"""
    nontopological_bases(F::Smith, BRS::BandRepSet; algorithm, verbose)
    nontopological_bases(sgnum::Integer, D::Integer=3; kwargs...)

Computes the "non-topological" Hilbert basis associated with a Smith normal form `F` of the
EBR matrix or from a space group number `sgnum`, returning a `SymBasis` structure. 
The returned basis is a non-negative integer coefficient basis for all non-topological band
structures (i.e. both trivial and fragile-topological).

If the method is called with `sgnum::Integer`, the underlying `BandRepSet` is also returned.

For possible keyword arguments, see `compatibility_bases(..)`.
"""
function nontopological_bases(F::Smith, BRS::BandRepSet;
                              algorithm::String="DualMode", verbose::Bool=false)
    # To find _all_ nontopological bases we build a cone subject to the inequalities 
    # (SΛy)ᵢ ≥ 0 with yᵢ ∈ ℤ, which automatically excludes topological cases (since they
    # correspond to rational yᵢ)
    dᵇˢ = count(!iszero, F.SNF)  # "Dimensionality" of band structure
    S = @view F.S[:,OneTo(dᵇˢ)]  # All the nontrivial conditions on zⱼ
    Λ = @view F.SNF[OneTo(dᵇˢ)]  # Nonzero invariant factors of Smith normal decomposition
    SΛ = S .* Λ' # Equivalent to S*diagm(dᵇˢ, dᵇˢ, Λ); all the nontrivial conditions on yᵢ
    
    C_nontopo = PyNormaliz.Cone(inequalities = SΛ)    # Cone consistent with SᵢⱼΛⱼyⱼ ≥ 0
    verbose && C_nontopo.setVerbose(true)             # Whether to print progress info
    C_nontopo.Compute("HilbertBasis", algorithm)      # Compute¹ the Hilbert basis
    ysᴴ_nontopo = transpose(C_nontopo.HilbertBasis()) # Hilbert basis vectors in 𝐲-space

    nsᴴ_nontopo = SΛ*ysᴴ_nontopo                      # Hilbert basis vectors in 𝐧-space

    return SymBasis(nsᴴ_nontopo, BRS, false) # Bases of nontopological states in 𝐧-space
end

# Convenience accessors from a space group number and dimensionality alone
for f in (:compatibility_bases, :nontopological_bases)
    @eval begin
        function $f(sgnum::Integer, D::Integer=3; 
                    algorithm::String="DualMode", verbose::Bool=false,
                    spinful::Bool=false, timereversal::Bool=true, allpaths::Bool=false)
            BRS = bandreps(sgnum, D; allpaths=allpaths, spinful=spinful, timereversal=timereversal)
            B   = matrix(BRS, true) # Matrix with columns of EBRs.
            F   = smith(B)          # Smith normal decomposition of B

            return $f(F, BRS, algorithm=algorithm, verbose=verbose), BRS
        end
    end
end

"""
$(TYPEDSIGNATURES)

Compute the trivial and fragile indices of a _nontopological_ `SymBasis`, `sb_nontopo`, by 
determining whether or not each has a positive-coefficient expansion in elementary band
representations (EBRs).
The EBRs are given either through a `BRS::BandRepSet` or through its matrix representation
`B = matrix(BRS, includedim)`.

Note that both `sb_nontopo` and `B` must reference the same output space: in other words, if
the `sb_nontopo` includes a filling element, `includedim` must set to `true` for `B`; 
otherwise `false`.

Returns trivial indices, `trivial_idxs`, and fragile indices, `fragile_idxs`, into the basis
vectors in `sb_nontopo`.
"""
function split_fragiletrivial_bases(sb_nontopo::SymBasis, B::AbstractMatrix)
    if sb_nontopo.compatbasis
        throw(DomainError(sb_nontopo, "provided `SymBasis` must have `compatbasis=false`"))
    end
    # Every vector of nsᴴ_nontopo that has a non-negative integer coefficient expansion in
    # EBRs, i.e. in B, represents a trivial basis element. All other elements represent 
    # fragile basis elements. We can just test whether such a solution exists through
    # constrained optimization, and then partition into trivial and fragile categories
    Nᴱᴮᴿ = size(B, 2)
    trivial_idxs = Int[]; fragile_idxs = Int[]
    for (j, nᴴ_nontopoʲ) in enumerate(sb_nontopo)
        m = Model(GLPK.Optimizer)
        @variable(m, c[1:Nᴱᴮᴿ] >= 0, Int)
        @constraint(m, B*c .== nᴴ_nontopoʲ)
        optimize!(m)

        # Check to see what the termination status (i.e. result) of the optimization was 
        status = termination_status(m)
        if status == MOI.OPTIMAL         # A feasible solution was found ⇒ trivial!
            push!(trivial_idxs, j)
        elseif status == MOI.INFEASIBLE  # No feasible solution exists   ⇒ fragile!
            push!(fragile_idxs, j)
        else
            throw("Unexpected termination status $status")
        end
    end

    return trivial_idxs, fragile_idxs
end
function split_fragiletrivial_bases(sb_nontopo::SymBasis, BRS::BandRepSet)
    Nirr, Nrows = length(irreplabels(sb_nontopo)), length(first(sb_nontopo))
    includedim = (Nrows == Nirr+1) ? true : false

    @assert (includedim || Nrows == Nirr)
    @assert Nirr == length(first(BRS)) "Non-equal number of irreps in SymBasis and BandRepSet"    

    split_fragiletrivial_bases(sb_nontopo, matrix(BRS, includedim))
end

# Footnotes:
# ¹⁾ We want to control the algorithm used to calculate the Hilbert basis in Normaliz: I've
# found that Normaliz usually picks the "PrimalMode" algorithm (not always though), despite
# the documentation stating that "DualMode" usually is better for cones defined by
# inequalities. Through tests, I've found that there is usually a speedup of order 2-10,
# if the "DualMode" algorithm is picked instead; occassionally, up to ×100s. The speedup 
# seems to be especially large for the systems that are (very) hard to solve (e.g. 131). 
# To force "DualMode", we use the method `C.Compute(<quantity-to-compute>, <algorithm>)`, 
# which then calculates `<quantity-to-compute>` and stores it in `C` (here, "HilbertBasis");
# it can then subsequently be retrieved from `C` with `C.HilbertBasis()`.
# Some situations remain effectively out of reach still: for bosons without time-reversal
# symmetry, there are a handful of SGs for which nontopological_bases(..) does not seem to
# terminate in a reasonable amount of time (15h+ tested).