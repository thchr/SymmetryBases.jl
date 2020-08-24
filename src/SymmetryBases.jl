module SymmetryBases

using Crystalline
using SmithNormalForm
using PyCall
using PrettyTables
using JuMP, GLPK
using Nemo # for `calc_topology`
using DocStringExtensions

import Base: OneTo, show, size, getindex, firstindex, lastindex, IndexStyle, length
import Crystalline: matrix, vecs, num, irreplabels, klabels, isspinful, istimeinvar

export SymBasis
export compatibility_bases, nontopological_bases, split_fragiletrivial_bases,
       has_posint_expansion, calc_detailed_topology, calc_topology, fillings,
       TopologyKind, trivial, nontrivial, fragile

# import the PyNormaliz library
const PyNormaliz = PyNULL()
function __init__() # https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
    copy!(PyNormaliz, pyimport("PyNormaliz"))
end

# -----------------------------------------------------------------------------------------

"""
    $(TYPEDEF)
    
with fields: $(TYPEDFIELDS)
"""
struct SymBasis <: AbstractVector{Vector{Int}}
    symvecs::Vector{Vector{Int}}
    irlabs::Vector{String}
    klabs::Vector{String}
    kvs::Vector{KVec}
    kv2ir_idxs::Vector{UnitRange{Int}} # pick k-point; find assoc. ir indices
    sgnum::Int
    spinful::Bool
    timeinvar::Bool
    allpaths::Bool
    compatbasis::Bool
end
function SymBasis(nsᴴ::AbstractMatrix{Int}, BRS::BandRepSet, compatbasis::Bool=true)
    kv2ir_idxs = [(f = irlab -> klabel(irlab)==klab; 
                   findfirst(f, BRS.irlabs):findlast(f, BRS.irlabs)) for klab in BRS.klabs]
    return SymBasis(collect(eachcol(nsᴴ)),
                    BRS.irlabs, BRS.klabs, BRS.kvs, kv2ir_idxs, 
                    BRS.sgnum, BRS.spinful, BRS.timeinvar, BRS.allpaths, compatbasis)
end

# accessors
matrix(sb::SymBasis) = hcat(sb.symvecs...)
vecs(sb::SymBasis)   = sb.symvecs
num(sb::SymBasis)    = sb.sgnum
irreplabels(sb::SymBasis) = sb.irlabs
klabels(sb::SymBasis)     = sb.klabs
isspinful(sb::SymBasis)   = sb.spinful
istimeinvar(sb::SymBasis) = sb.timeinvar
hasnonmax(sb::SymBasis)   = sb.allpaths
iscompatbasis(sb::SymBasis) = sb.compatbasis
fillings(sb::SymBasis)    = [nᴴ[end] for nᴴ in sb.symvecs]
length(sb::SymBasis)      = length(vecs(sb))

# define the AbstractArray interface for SymBasis
size(sb::SymBasis) = (length(sb),)
getindex(sb::SymBasis, keys...) = vecs(sb)[keys...]
firstindex(::SymBasis) = 1
lastindex(sb::SymBasis) = length(vecs(sb))
IndexStyle(::SymBasis) = IndexLinear()

# show method
function show(io::IO, ::MIME"text/plain", sb::SymBasis)
    Nⁱʳʳ = length(sb[1]) - 1

    # print a "title" line and the irrep labels
    println(io, iscompatbasis(sb) ? "Compatibility" : "Nontopological",
                " SymBasis (#", num(sb), "): ",
                length(sb), " Hilbert bases, sampling ",
                Nⁱʳʳ, " LGIrreps ",
                "(spin-", isspinful(sb) ? "½" : "1", " ",
                istimeinvar(sb) ? "w/" : "w/o", " TR)")

    k_idx = (i) -> findfirst(==(klabel(irreplabels(sb)[i])), klabels(sb)) # highlighters
    h_odd = Highlighter((data,i,j) -> i≤Nⁱʳʳ && isodd(k_idx(i)), crayon"light_blue")
    h_ν   = Highlighter((data,i,j) -> i==Nⁱʳʳ+1,                 crayon"light_yellow")

    pretty_table(io, 
        # table contents
        matrix(sb),
        # header
        eachindex(sb),
        # row names
        row_names = vcat(sb.irlabs, "ν"),
        # options/formatting/styling
        formatters = (v,i,j) -> iszero(v) ? "·" : string(v),
        vlines = [1,], hlines = [:begin, 1, Nⁱʳʳ+1, :end],
        row_name_alignment = :l,
        alignment = :c, 
        highlighters = (h_odd, h_ν), 
        header_crayon = crayon"bold"
        )

    # print k-vec labels
    print(io, "  KVecs (", hasnonmax(sb) ? "incl. non-maximal" : "maximal only", "): ")
    join(io, klabels(sb), ", ")
end

# -----------------------------------------------------------------------------------------

# All band structures can be written as 𝐧 = B𝐩 with pᵢ∈ℚ and nᵢ∈𝐍, and B a matrix whose 
# columns are EBRs. We can decompose B to a Smith normal form B = SΛT, such that all 
# allowable band structures can be written as 𝐧 = S𝐳. Here, S is an integer matrix with 
# elements Sᵢⱼ∈ℤ. To make nᵢ integer, we thus require zᵢ∈ℤ.

"""
    compatibility_bases(F::SmithNormalForm.Smith, BRS::BandRepSet; algorithm)
    compatibility_bases(sgnum::Integer; kwargs...)

Computes the Hilbert bases associated with a Smith normal form `F` of the EBR matrix or from
a space group number `sgnum`, which respects all compatibility relations. The resulting 
bases form a non-negative integer coefficient basis for all possible band structures.

If the method is called with `sgnum::Integer`, the associated `BandRepSet` is also returned.

Several keyword arguments `kwargs` are possible:

    - `algorithm::String`: controls the algorithm used by Normaliz to compute the Hilbert
    basis. Choices are `"DualMode"` (default) and `"PrimalMode"`
    - `spinful::Bool`: Use single- (`false`, default) or double-valued (`true`) irreps.
    - `timereversal::Bool`: Assume presence (`true`, default) or absence (`false`) of
    time-reversal symmetry.
    - `verbose::Bool`: whether to print progress info during the Normaliz computation
    (`false`, default).
"""
function compatibility_bases(F::SmithNormalForm.Smith, BRS::BandRepSet; 
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

    return SymBasis(nsᴴ, BRS, true), zsᴴ  # Bases of all valid symmetry vectors in 𝐧- and 𝐲-space
end

"""
    nontopological_bases(F::SmithNormalForm.Smith, BRS::BandRepSet; algorithm)
    nontopological_bases(sgnum::Integer; kwargs...)

Computes the "non-topological" Hilbert bases associated with a Smith normal form `F` of the
EBR matrix or from a space group number `sgnum`, which forms a non-negative, integer span of
all non-topological band structures (i.e. both trivial and fragile-topological).

If the method is called with `sgnum::Integer`, the associated `BandRepSet` is also returned.

For possible keyword arguments `kwargs`, see `compatibility_bases(..)`.
"""
function nontopological_bases(F::SmithNormalForm.Smith, BRS::BandRepSet;
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

    return SymBasis(nsᴴ_nontopo, BRS, false), ysᴴ_nontopo # Bases of nontopological states
end

# Convenience accessors from a space group number alone
for f in (:compatibility_bases, :nontopological_bases)
    @eval begin
        function $f(sgnum::Integer; algorithm::String="DualMode", spinful::Bool=false,
                                    timereversal::Bool=true)
            BRS = bandreps(sgnum; allpaths=false, spinful=spinful, timereversal=timereversal)
            B   = matrix(BRS, true)    # Matrix with columns of EBRs.
            F   = Crystalline.smith(B) # Smith normal decomposition of B

            return $f(F, BRS, algorithm=algorithm)..., BRS
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
        throw(DomainError(sb_nontopo, "Specified SymBasis must have compatbasis=false"))
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


# ---------------------------------------------------------------------------------------- #
#                  Methods to test the topology of a given symmetry vector n               #
# ---------------------------------------------------------------------------------------- #

"""
    $(TYPEDSIGNATURES)

Check whether a positive-integer coefficient expansion exists for `n` in the basis of the 
columns of `M`, i.e. whether there exists ``cᵢ∈ℕ`` (``=0,1,2,...``)  such that ``Mc=n``.
"""
function has_posint_expansion(n::AbstractVector{<:Integer}, M::AbstractMatrix{<:Integer})
    N = size(M, 2)
    # feasibility problem with positive-integer variables, subject to the condition Mc = n
    m = Model(GLPK.Optimizer)
    @variable(m, c[1:N] >= 0, Int)
    @constraint(m, M*c .== n)
    # try to solve the model
    optimize!(m)

    return m
end

@enum TopologyKind trivial=0 nontrivial=1 fragile=2

"""
    calc_detailed_topology(n, nontopo_M, trivial_M, M=nothing) -> ::TopologyKind

Evaluate whether a given (valid, i.e. regular) symmetry vector represents a band-combination
that is trivial, nontrivial, or fragile from a symmetry perspective.
This is done by comparing against the nontopological, trivial, and full symmetry bases
`nontopo_M`, `trivial_M`, and `M`, respectively, provided as matrices with columns of
symmetry basis vectors (i.e. checks whether a valid expansion exists in each).

If `trivial_M` is given as `nothing`, it is taken to imply that it is equal to `nontopo_M`,
meaning that there are no fragile phases. 
In this case, however, [`calc_topology`](@ref) will generally be far more efficient.

If `M` is _not_ `nothing` (i.e. a matrix representing the full symmetry basis), an 
additional sanity/safety check is carried out: otherwise not. Otherwise not necessary.

Returns a member value of the `TopologyKind` Enum (`trivial`, `nontrivial`, or 
`fragile`).
"""
function calc_detailed_topology(n::AbstractVector{<:Integer},
            nontopo_M::AbstractMatrix{<:Integer},
            trivial_M::Union{Nothing, AbstractMatrix{<:Integer}},
            M::Union{Nothing, <:AbstractMatrix{<:Integer}}=nothing)
    # check whether expansion exists in nontopological basis
    nontopo_m = has_posint_expansion(n, nontopo_M)
    # TODO: Maybe this should be done with `calc_topology` instead (since `cansolve` is
    #       significantly faster than JuMP's feasibility optimization)? Subsequent assesment
    #       of fragility should still be done with JuMP and `has_posint_expansion`. Would
    #       not necessarily require EBR as input; we could just use the `nontopo_M` as well
    #       (though, generally, the EBR basis would be more compact and hence better).

    # check to see what the termination status (i.e. result) of the optimization was 
    status′ = termination_status(nontopo_m)
    if status′ == MOI.OPTIMAL           # feasible solution found     ⇒ trivial/fragile
        # determine whether trivial or fragile
        if !isnothing(trivial_M)
            trivial_m = has_posint_expansion(n, trivial_M)
            if termination_status(trivial_m) ≠ MOI.OPTIMAL
                # expansion in trivial-only basis elements impossible ⇒ fragile
                return fragile
            end
        end
                # expansion in trivial-only elements feasible         ⇒ trivial
        return trivial
        
    elseif status′ == MOI.INFEASIBLE    # infeasible problem          ⇒ nontrivial
        if !isnothing(M) # do a sanity-check to verify that expansion exists in full basis
            m = has_posint_expansion(n, M)
            if termination_status(m) ≠ MOI.OPTIMAL
                throw("It must be possible to find an expansion in the full basis")
            end
        end

        return nontrivial

    else
        throw(DomainError(termination_status(nontopo_m), 
            "Unexpected termination status of nontopo_m: expected OPTIMAL or INFEASIBLE"))
    end

    return 
end

function calc_detailed_topology(n::AbstractVector{<:Integer}, 
            nontopo_sb::SymBasis, BRS::BandRepSet, sb::Union{Nothing, SymBasis}=nothing)
    
    nontopo_M = matrix(nontopo_sb)
    
    trivial_idxs, fragile_idxs = split_fragiletrivial_bases(nontopo_sb, matrix(BRS, true))
    can_be_fragile = !isempty(fragile_idxs)
    trivial_M = can_be_fragile ? (@view nontopo_M[:, trivial_idxs]) : nothing
    
    M = sb === nothing ? nothing : matrix(sb)

    return calc_detailed_topology(n, nontopo_M, trivial_M, M)
end

function calc_detailed_topology(n::AbstractVector{<:Integer}, sgnum::Integer; 
            spinful::Bool=false, timereversal::Bool=true)
    nontopo_sb, _, BRS = nontopological_bases(sgnum; spinful=spinful, timereversal=timereversal)
    sb, _, _           = compatibility_bases(sgnum; spinful=spinful, timereversal=timereversal)

    return calc_detailed_topology(n, nontopo_sb, BRS, sb)
end


# -----------------------------------------------------------------------------------------
# Trivial/nontrivial solution topology via BandRepSet and Nemo

"""
    $(TYPEDSIGNATURES)

Evaluate whether a given (valid, i.e. regular) symmetry vector `n` represents a
band-combination that is trivial or nontrivial from a symmetry perspective, i.e. whether it
has an integer-coefficient expansion in the elementary band representation (EBR) basis or
not (i.e. a rational-coefficient expansion).

No distinction is made between fragile and trivial symmetry vectors (see
[`calc_detailed_topology`](@ref)).

The EBR basis can be provided as `::BandRepSet`, `::Matrix{<:Integer}`, or `::fmpz_mat` (a
Nemo.jl-specific type).
The length of `n` must equal the EBR basis' number of irreps or the number of irreps plus 1
(i.e. include the band filling).
Evaluation of whether an integer-coefficient expansion exists is performed via Nemo.jl's
`cansolve`.

Returns a member value of the `TopologyKind` Enum (`trivial` or `nontrivial`).
"""
function calc_topology(n::AbstractVector{<:Integer}, Bℤ::fmpz_mat)
    nℤ = MatrixSpace(ZZ, length(n), 1)(n)
    # test whether an integer coefficient expansion exists for `nℤ` in the EBR basis `Bℤ`
    solvable, _ = cansolve(Bℤ, nℤ)

    return solvable ? trivial : nontrivial
end

function calc_topology(n::AbstractVector{<:Integer}, B::Matrix{<:Integer})
    length(n) == size(B, 2) || throw(DimensionMismatch())
    Bℤ = MatrixSpace(ZZ, size(B)...)(B)

    return calc_topology(n, Bℤ)
end

function calc_topology(n::AbstractVector{<:Integer}, BRS::BandRepSet)
    Nirr, Nn = length(irreplabels(BRS)), length(n)
    includedim = (Nn == Nirr+1) ? true : false
    @assert (includedim || Nn == Nirr)

    calc_topology(n, matrix(BRS, includedim))
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
end # module