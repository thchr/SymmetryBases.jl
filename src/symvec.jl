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

Returns a member value of the `TopologyKind` Enum (`TRIVIAL`, `NONTRIVIAL`, or 
`FRAGILE`).
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
                return FRAGILE
            end
        end
        # expansion in trivial-only elements feasible                ⇒ trivial
        return TRIVIAL
        
    elseif status′ == MOI.INFEASIBLE    # infeasible problem          ⇒ nontrivial
        if !isnothing(M) # do a sanity-check to verify that expansion exists in full basis
            m = has_posint_expansion(n, M)
            if termination_status(m) ≠ MOI.OPTIMAL
                error("It must be possible to find an expansion in the full basis")
            end
        end

        return NONTRIVIAL

    else
        throw(DomainError(status′, 
            "unexpected optimization termination status: expected OPTIMAL or INFEASIBLE"))
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

function calc_detailed_topology(n::AbstractVector{<:Integer}, sgnum::Integer, D::Integer=3;
            spinful::Bool=false, timereversal::Bool=true, allpaths::Bool=false)
    BRS = bandreps(sgnum, D, spinful=spinful, timereversal=timereversal, allpaths=allpaths)
    F   = smith(matrix(BRS, true))   
    nontopo_sb = nontopological_bases(F, BRS)
    sb         = compatibility_bases(F, BRS)

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

    return solvable ? TRIVIAL : NONTRIVIAL
end

function calc_topology(n::AbstractVector{<:Integer}, B::Matrix{<:Integer})
    length(n) == size(B, 1) || throw(DimensionMismatch("sizes of n and B are inconsistent"))
    Bℤ = MatrixSpace(ZZ, size(B)...)(B)

    return calc_topology(n, Bℤ)
end

function calc_topology(n::AbstractVector{<:Integer}, BRS::BandRepSet)
    Nirr, Nn = length(irreplabels(BRS)), length(n)
    includedim = (Nn == Nirr+1) ? true : false
    @assert (includedim || Nn == Nirr)

    calc_topology(n, matrix(BRS, includedim))
end

# -----------------------------------------------------------------------------------------
# test whether a band grouping respects compatibility relations, i.e. are in {BS}
# (another way would be to use a `SymBasis` and `has_posint_expansion` to test whether and
# integer conical combination exists; but that is _much_ slower (~30-150× at least))
"""
$(TYPEDSIGNATURES)

Test whether a symmetry vector `n` is a valid band grouping, i.e. whether it fulfils all
compatibility relations in the Brillouin zone and is non-negative. That is, test whether
`n` belong to the set of physical band structures {BS}.

## Implementation

Belongingness to {BS} is tested by comparing to a set of elementary band representations
(EBRs), provided either as a `BandRepSet`, a `Matrix{<:Integer}`, or a `Smith`
decomposition.
A symmetry vector ``\\mathbf{n}`` is in {BS} if

``
\\tilde{\\mathbf{S}}\\tilde{\\mathbf{S}}^{-1}\\mathbf{n} = \\mathbf{n}
``

where ``\\tilde{\\mathbf{S}}`` (``\\tilde{\\mathbf{S}}^{-1}``) denotes the nonsingular
columns (rows) of ``\\mathbf{S}`` (``\\mathbf{S}^{-1}``) in the Smith normal decomposition
of the EBR matrix ``\\mathbf{A} = \\mathbf{S}\\boldsymbol{\\Lambda}\\mathbf{T}``.

## Example
```julia-repl
julia> sb, brs = compatibility_bases(22, 3); # from Crystalline.jl
julia> n = sb[1];

# test a valid symmetry vector
julia> isbandstruct(n, brs)
true

# test an invalid symmetry vector
julia> n′ = copy(n);
julia> n′[1] += 1;              # invalid modification
julia> isbandstruct(n′, brs)
false
```
"""
function isbandstruct(n::AbstractVector{<:Integer}, F::Smith)
    dᵇˢ = count(!iszero, F.SNF)
    S   = @view F.S[:,OneTo(dᵇˢ)]     # relevant columns of S only
    S⁻¹ = @view F.Sinv[OneTo(dᵇˢ), :] # relevant rows of S⁻¹ only
    LHS = (S*S⁻¹ - I)*n

    return all(==(0), LHS)
end
isbandstruct(n::AbstractVector{<:Integer}, B::Matrix{<:Integer}) = isbandstruct(n, smith(B))
isbandstruct(n::AbstractVector{<:Integer}, BRS::BandRepSet) = isbandstruct(n, matrix(BRS, true))