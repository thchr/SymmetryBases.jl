
# All band structures can be written as 𝐧 = B𝐩 with pᵢ∈ℚ and nᵢ∈𝐍, and B a matrix whose 
# columns are EBRs. We can decompose B to a Smith normal form B = SΛT, such that all 
# allowable band structures can be written as 𝐧 = S𝐳. Here, S is an integer matrix with 
# elements Sᵢⱼ∈ℤ. To make nᵢ integer, we thus require zᵢ∈ℤ.

"""
    compatibility_basis([F::Smith,] brs::Collection{<:BandRep}; kws...)
    compatibility_basis(sgnum::Integer, D::Integer=3; kwargs...)

Computes the Hilbert basis associated with a Smith normal form `F` of the EBR matrix or from
a space group number `sgnum`, which respects all compatibility relations, returning a
`SymBasis` structure. The returned basis is a non-negative integer coefficient basis
for all possible band structures {BS}.

If the method is called with `sgnum::Integer`, the underlying band representations are
also returned.

## Keyword arguments

- `algorithm::String`: controls the algorithm used by Normaliz to compute the Hilbert
basis. Choices are `"DualMode"` (default) and `"PrimalMode"`
- `spinful`: Use single- (`Val(false)`, default) or double-valued (`Val(true)`) irreps.
A plain `Bool` is also accepted, but is not type-stable.
- `timereversal::Bool`: Assume presence (`true`, default) or absence (`false`) of
time-reversal symmetry.
- `verbose::Bool`: whether to print progress info during the Normaliz computation
(`false`, default).
"""
function compatibility_basis(F::Smith, brs::Collection{<:BandRep};
                             algorithm::String="DualMode", verbose::Bool=false)
    # To restrict nᵢ to only positive integers, i.e. ℕ, the values of zᵢ must be such that 
    # ∑ⱼ Sᵢⱼzⱼ ≥ 0. This defines a set of inequalities, which in turn defines a polyhedral
    # integer cone. This is where (Py)Normaliz comes in.
    dᵇˢ = count(!iszero, F.SNF)           # "Dimensionality" of band structure
    S = @view F.S[:, 1:dᵇˢ]               # All the nontrivial conditions on zⱼ

    C = PyNormaliz.Cone(inequalities = eachrow(S)) # Cone consistent w/ Sᵢⱼzⱼ ≥ 0
    verbose && C.SetVerbose(true)         # Whether to print progress info
    C.Compute("HilbertBasis", algorithm)  # Compute¹ the Hilbert basis
    zsᴴ_py = C.HilbertBasis()             # Python list of lists (each element a Hilbert
                                          # basis vectors in 𝐳-space)
    zsᴴ_vs = pyconvert(Vector{Vector{Int}}, zsᴴ_py) # Convert to Julia vector of vectors
    zsᴴ = stack(zsᴴ_vs)                   # Columns are now Hilbert basis vectors (𝐳-space)
    nsᴴ  = S*zsᴴ                          # Columns are Hilbert basis vectors in 𝐧-space

    return SymBasis(nsᴴ, brs, true)       # Bases of all valid symmetry vectors in 𝐧-space
end

"""
    nontopological_basis([F::Smith,] brs::Collection{<:BandRep}; kws...)
    nontopological_basis(sgnum::Integer, D::Integer=3; kwargs...)

Computes the "nontopological" Hilbert basis associated with a Smith normal form `F` of the
EBR matrix or from a space group number `sgnum`, returning a `SymBasis` structure. 
The returned basis is a non-negative integer coefficient basis for all non-topological band
structures {AI+F} (i.e. both trivial and fragile-topological).

If the method is called with `sgnum::Integer`, the underlying band representations are
also returned.

For possible keyword arguments, see `compatibility_basis(..)`.
"""
function nontopological_basis(F::Smith, brs::Collection{<:BandRep};
                              algorithm::String="DualMode", verbose::Bool=false)
    # To find the nontopological basis we build a cone subject to the inequalities 
    # (SΛy)ᵢ ≥ 0 with yᵢ ∈ ℤ, which automatically excludes topological cases (since they
    # correspond to rational yᵢ)
    dᵇˢ = count(!iszero, F.SNF)  # "Dimensionality" of band structure
    S = @view F.S[:,OneTo(dᵇˢ)]  # All the nontrivial conditions on zⱼ
    Λ = @view F.SNF[OneTo(dᵇˢ)]  # Nonzero invariant factors of Smith normal decomposition
    SΛ = S .* Λ' # Equiv. to S*diagm(dᵇˢ, dᵇˢ, Λ); all nontrivial conditions on yᵢ
    
    C_nontopo = PyNormaliz.Cone(inequalities = eachrow(SΛ)) # Cone consistent w/ SᵢⱼΛⱼyⱼ ≥ 0
    verbose && C_nontopo.SetVerbose(true)          # Whether to print progress info
    C_nontopo.Compute("HilbertBasis", algorithm)   # Compute¹ the Hilbert basis
    ysᴴ_nontopo_py = C_nontopo.HilbertBasis()      # Hilbert basis vectors in 𝐲-space
    ysᴴ_nontopo_vs = pyconvert(Vector{Vector{Int}}, ysᴴ_nontopo_py) # Convert to Julia
    ysᴴ_nontopo = stack(ysᴴ_nontopo_vs)            # Hilbert basis vectors as cols (𝐲-space)

    nsᴴ_nontopo = SΛ*ysᴴ_nontopo                   # Hilbert basis vectors in 𝐧-space

    return SymBasis(nsᴴ_nontopo, brs, false) # Bases of nontopological states in 𝐧-space
end

# Convenience accessors from a space group number and dimensionality alone
for f in (:compatibility_basis, :nontopological_basis)
    @eval begin
        function $f(sgnum::Integer, D::Integer=3;
                    algorithm::String="DualMode", verbose::Bool=false,
                    spinful = Val(false), timereversal::Bool=true, allpaths::Bool=false)
            brs = bandreps(sgnum, D; allpaths, spinful, timereversal)
            return $f(brs; algorithm, verbose), brs
        end
        $f(brs::Collection{<:BandRep}; kws...) = $f(smith(stack(brs)), brs; kws...)
    end
end

"""
$(TYPEDSIGNATURES)

Compute the trivial and fragile indices of a _nontopological_ `SymBasis`, `sb_nontopo`, by 
determining whether or not each has a positive-coefficient expansion in elementary band
representations (EBRs).
The EBRs are given either through a `brs::Collection{<:BandRep}` or through its matrix
representation `B = stack(brs)`.

Returns trivial indices `trivial` and fragile indices `fragile` (indexing into the basis
vectors in `sb_nontopo`) as a named tuple with the corresponding field names.
"""
function split_fragiletrivial(ns::AbstractVector{<:AbstractVector{<:Integer}},
                              B::AbstractMatrix)
    # Every (not-nontrivial) vector of `ns` that has a non-negative integer coefficient
    # expansion in EBRs, i.e. in `B`, represents a trivial basis element. All other elements
    # represent fragile basis elements. We can just test whether such a solution exists
    # through constrained optimization, and then partition into trivial and fragile 
    # categories
    # NOTE: This method explicitly assumes that it is _never_ provided any nontrivial
    #       symmetry vectors. If it is, it will misclassify them as fragile.
    # TODO: Maybe just split into all variants, including nontrivial?
    Nᴱᴮᴿ = size(B, 2)
    trivial_idxs = Int[]; fragile_idxs = Int[]
    for (j, nʲ) in enumerate(ns)
        m = Model(GLPK.Optimizer)
        @variable(m, c[1:Nᴱᴮᴿ] >= 0, Int)
        @constraint(m, B*c .== nʲ)
        optimize!(m)

        # Check to see what the termination status (i.e. result) of the optimization was 
        status = termination_status(m)
        if status == MOI.OPTIMAL         # A feasible solution was found ⇒ trivial!
            push!(trivial_idxs, j)
        elseif status == MOI.INFEASIBLE  # No feasible solution exists   ⇒ fragile!
            push!(fragile_idxs, j)
        else
            error("Unexpected termination status $status")
        end
    end

    return (trivial=trivial_idxs, fragile=fragile_idxs)
end
function split_fragiletrivial(sb_nontopo::SymBasis, B::AbstractMatrix)
    if sb_nontopo.compatbasis
        throw(DomainError(sb_nontopo, "provided `SymBasis` must have `compatbasis=false`"))
    end
    return split_fragiletrivial(parent(sb_nontopo), B)
end
function split_fragiletrivial(sb_nontopo::SymBasis, brs::Collection{<:BandRep})
    return split_fragiletrivial(sb_nontopo, stack(brs))
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
# symmetry, there are a handful of SGs for which nontopological_basis(..) does not seem to
# terminate in a reasonable amount of time (15h+ tested).