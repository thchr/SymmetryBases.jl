# ---------------------------------------------------------------------------------------- #
#                  Methods to test the topology of a given symmetry vector n               #
# ---------------------------------------------------------------------------------------- #

function _throw_incompatible_or_nonphysical(n)
    error(DomainError(n, "`n` is not a physically realizable band grouping"))
end

"""
    $(TYPEDSIGNATURES)

Return whether `n` includes the connectivity as an element by comparing with size of `BRS`.
"""
function includes_connectivity(n::AbstractVector{<:Integer}, BRS::BandRepSet)
    Nirr, Nn = length(irreplabels(BRS)), length(n)
    if Nn == Nirr+1
        return true
    elseif N == Nirr
        return false
    else 
        error(DimensionMismatch("incompatible dimensions of `n` and `BRS`"))
    end
end

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
    calc_detailed_topology(n, B::Matrix{<:Integer}, [F::Smith=smith(B)]; allow_nonphysical)
    calc_detailed_topology(n, BRS::BandRepSet; allow_nonphysical)
    calc_detailed_topology(n, sgnum::Integer, [D::Integer=3]; kwargs...) 

Return whether a integer symmetry vector `n` is topologically trivial, nontrivial, or
fragile (to the extent that this classification is symmetry-identifiable). Returns a
value from the Enum [`TopologyKind`](@ref) (`TRIVIAL`, `NONTRIVIAL`, or `FRAGILE`).

## Implementation
The Smith normal decomposition `F` of the EBR matrix `B` (or, equivalently, a provided
`BandRepSet`) is used to first test whether `n` is nontrivial or not stably nontrivial
(i.e. trivial or fragile) using [`calc_topology`](@ref).
In the latter case, we resolve triviality vs. fragility by subsequently checking whether
`n` has a non-negative expansion in the EBR basis using [`has_posint_expansion`](@ref).

This approach is equivalent to checking whether `n` has a non-negative expansion only in
the compatibility basis (⇒ nontrivial), in the preceding *and* in the nontopological basis
(⇒ fragile), or in the preceding *and* in the EBR basis (⇒ trivial) - but is much faster
(due to not solving multiple non-negative feasibility questions).

## Keyword arguments
If `n` is not a compatible band structure (i.e., if `isbandstruct(n, BRS) = false`), an
error is thrown. This behavior can be controlled by two boolean keyword arguments:

- `allow_incompatible` (`false`): if `true`, disables the compatibility check entirely.
- `allow_nonphysical` (`false`): if `true`, allows negative symmetry content, but maintain
  requirement that `n` respects the compatibilty relations in an algebraic sense.
"""
function calc_detailed_topology(
            n::AbstractVector{<:Integer},
            B::AbstractMatrix{<:Integer},
            F::Smith=smith(B);
            allow_incompatible::Bool=false,
            allow_nonphysical::Bool=false)

    coarse_topo = calc_topology(n, F; allow_incompatible, allow_nonphysical) # checks `isbandstruct` as well
    if coarse_topo == TRIVIAL       # ⇒ trivial/fragile
        trivial_m = has_posint_expansion(n, B)
        if termination_status(trivial_m) ≠ MOI.OPTIMAL
            # expansion in trivial-only basis elements impossible ⇒ fragile
            return FRAGILE          # ⇒ fragile
        else
            return TRIVIAL          # ⇒ trivial
        end
    else                            # ⇒ nontrivial
        return coarse_topo
    end
end

function calc_detailed_topology(n::AbstractVector{<:Integer}, BRS::BandRepSet; kws...)
    B = matrix(BRS; includedim=includes_connectivity(n, BRS))
    return calc_detailed_topology(n, B; kws...)
end

function calc_detailed_topology(
            n::AbstractVector{<:Integer},
            sgnum::Integer,
            D::Integer=3;
            spinful::Bool=false,
            timereversal::Bool=true,
            allpaths::Bool=false,
            kws...)

    BRS = bandreps(sgnum, D; spinful, timereversal, allpaths)
    return calc_detailed_topology(n, BRS; kws...)
end


# To-be-removed-variants below:
"""
    calc_detailed_topology(n, nontopo_M, trivial_M, M=nothing) -> ::TopologyKind

!!! warning
    This method signature is deprecated and will be removed in a future version: use 
    [`calc_detailed_topology(::AbstractVector, ::AbstractMatrix, ::Smith)`](@ref) instead.

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
    # TODO: Remove this method: it is far slower than the `(n, ::Matrix, Smith)` variant
    #       and no better (but quite a bit more complicated)

    # check whether expansion exists in nontopological basis
    nontopo_m = has_posint_expansion(n, nontopo_M)

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
end

# TODO: Remove method
function calc_detailed_topology(
            n::AbstractVector{<:Integer}, 
            nontopo_sb::SymBasis,
            BRS::BandRepSet,
            sb::Union{Nothing, SymBasis}=nothing)

    nontopo_M = matrix(nontopo_sb)
    
    trivial_idxs, fragile_idxs = split_fragiletrivial(nontopo_sb, matrix(BRS; includedim=true))
    can_be_fragile = !isempty(fragile_idxs)
    trivial_M = can_be_fragile ? (@view nontopo_M[:, trivial_idxs]) : nothing
    
    M = sb === nothing ? nothing : matrix(sb)

    return calc_detailed_topology(n, nontopo_M, trivial_M, M)
end


# -----------------------------------------------------------------------------------------
# Trivial/nontrivial solution topology via Smith/BandRepSet

@doc """
$(TYPEDSIGNATURES)

Return whether a symmetry vector `n` is a band-combination that is trivial or nontrivial
from a symmetry perspective, i.e. whether it has an integer-coefficient expansion in the
elementary band representation (EBR) basis or not (i.e. a rational-coefficient expansion).
Returns a value from the Enum [`TopologyKind`](@ref) (`TRIVIAL` or `NONTRIVIAL`).

No distinction is made between fragile and trivial symmetry vectors: i.e., a `TRIVIAL`
return value may in fact be a `FRAGILE` state on more careful inspection (see
[`calc_detailed_topology`](@ref)).

## Input

The EBR basis can be provided as `::BandRepSet`, `::Matrix{<:Integer}`, or a `Smith`
decomposition.
The length of `n` must equal the EBR basis' number of irreps or the number of irreps plus 1
(i.e. include the band connectivity).

## Implementation

We check whether an integer-coefficient expansion exists via the Smith normal decomposition
of the EBR matrix ``\\mathbf{B} = \\mathbf{S}\\boldsymbol{\\Lambda}\\mathbf{T}``. If

```math
    (\\mathbf{S}^{-1}\\mathbf{n})_j = 0 \\mod \\lambda_j
```

for all ``j = 1, \\ldots, d^{\\text{bs}}`` (``d^{\\text{bs}}`` is the number of non-zero
diagonal elements of ``\\boldsymbol{\\Lambda}``, i.e. the invariant factors of
``\\mathbf{B}``), there exists a solution to ``\\mathbf{B}\\mathbf{c} = \\mathbf{n}`` with
integer coefficients ``c_j \\in \\mathbb{Z}``.

## Keyword arguments
If `n` is not a compatible band structure (i.e., if `isbandstruct(n, [...]) = false`), an
error is thrown. This behavior can be controlled by two boolean keyword arguments:

- `allow_incompatible` (`false`): if `true`, disables the compatibility check entirely.
- `allow_nonphysical` (`false`): if `true`, allows negative symmetry content, but maintain
  requirement that `n` respects the compatibilty relations in an algebraic sense.
"""
function calc_topology(
            n::AbstractVector{<:Integer},
            F::Smith;
            allow_incompatible::Bool=false,
            allow_nonphysical::Bool=false)

    if !allow_incompatible && !isbandstruct(n, F; allow_nonphysical)
        _throw_incompatible_or_nonphysical(n)
    end
    
    S⁻¹, Λ = F.Sinv, F.SNF # Λ = [λ₁, …, λ_{dᵇˢ}, 0, …, 0]
    dᵇˢ = count(!iszero, Λ)

    # n is trivial if (S⁻¹n)ⱼ = 0 mod λⱼ for j = 1, …, dᵇˢ. This is equivalent to checking
    # whether there exists an integer coefficient expansion for `n` in the EBR basis that
    # `F` represents (i.e., whether `Nemo.cansolve(B, n) == true`) but faster.
    # We do the matrix-vector product row-wise to check `mod(S⁻¹[1:dᵇˢ]*n)[i], Λ[i]) = 0`
    # for `i ∈ 1:dᵇˢ` without allocating unnecessarily
    is_trivial = all(1:dᵇˢ) do i
        Λᵢ = Λ[i]
        Λᵢ == 1 && return true # fast path: `mod(x, 1) = 0` for all integer `x`.
        S⁻¹ᵢ = @view S⁻¹[i,:]
        mod(dot(S⁻¹ᵢ, n), Λᵢ) == 0
    end
    return is_trivial ? TRIVIAL : NONTRIVIAL
end

function calc_topology(n::AbstractVector{<:Integer}, B::AbstractMatrix{<:Integer}; kws...)
    length(n) == size(B, 1) || throw(DimensionMismatch("incompatible dimensions of `n` and `B`"))
    return calc_topology(n, smith(B); kws...)
end

function calc_topology(n::AbstractVector{<:Integer}, BRS::BandRepSet; kws...)
    B = matrix(BRS; includedim=includes_connectivity(n, BRS))
    return calc_topology(n, B; kws...)
end

# -----------------------------------------------------------------------------------------
# test whether a band grouping respects compatibility relations, i.e. are in {BS}
# (another way would be to use a `SymBasis` and `has_posint_expansion` to test whether and
# integer conical combination exists; but that is _much_ slower (~30-150× at least))
@doc """
$(TYPEDSIGNATURES)

Test whether a symmetry vector `n` is a valid band grouping, i.e. whether it fulfils all
compatibility relations in the Brillouin zone and is non-negative. That is, test whether
`n` belong to the set of physical band structures {BS}.

## Keyword arguments

- `allow_nonphysical` (`false`): if `true`, allows negative symmetry content.
  This can be relevant if `n` contains negative content that may nevertheless respect the
  compatibility relations in a strictly algebraic sense.

## Implementation

Belongingness to {BS} is tested by comparing to a set of elementary band representations
(EBRs), provided either as a `BandRepSet`, a `Matrix{<:Integer}`, or a `Smith`
decomposition.
A symmetry vector ``\\mathbf{n}`` is in {BS} if

```math
    \\tilde{\\mathbf{S}}\\tilde{\\mathbf{S}}^{-1}\\mathbf{n} = \\mathbf{n}
```

where ``\\tilde{\\mathbf{S}}`` (``\\tilde{\\mathbf{S}}^{-1}``) denotes the nonsingular
columns (rows) of ``\\mathbf{S}`` (``\\mathbf{S}^{-1}``) in the Smith normal decomposition
of the EBR matrix ``\\mathbf{A} = \\mathbf{S}\\boldsymbol{\\Lambda}\\mathbf{T}``.

## Examples

```julia-repl
julia> sb, brs = compatibility_basis(22, 3); # from Crystalline.jl
julia> n = sb[1];

# test a valid symmetry vector
julia> isbandstruct(n, brs)
true

# test an invalid symmetry vector
julia> n′ = copy(n);
julia> n′[1] += 1;                   # invalid modification
julia> isbandstruct(n′, brs)
false

# test a symmetry vector with negative content
julia> n′′ = sb[1] + sb[2] - sb[3];  # contains negative elements
julia> isbandstruct(n′′, brs)
false
julia> isbandstruct(n′′, brs; allow_nonphysical=true)
true
```
"""
function isbandstruct(
            n::AbstractVector{<:Integer},
            F::Smith;
            allow_nonphysical::Bool=false)

    # check non-negativity
    (allow_nonphysical || all(≥(0), n)) || return false

    # check compatibility relations
    dᵇˢ = count(!iszero, F.SNF)
    S̃   = @view F.S[:,OneTo(dᵇˢ)]     # relevant columns of S only
    S̃⁻¹ = @view F.Sinv[OneTo(dᵇˢ), :] # relevant rows of S⁻¹ only
    return S̃*(S̃⁻¹*n) == n
end
isbandstruct(n::AbstractVector{<:Integer}, B::Matrix{<:Integer}; kws...) = isbandstruct(n, smith(B); kws...)
isbandstruct(n::AbstractVector{<:Integer}, BRS::BandRepSet; kws...) = isbandstruct(n, matrix(BRS; includedim=true); kws...)

# -----------------------------------------------------------------------------------------
# Stable topological indices

@doc """
$(TYPEDSIGNATURES)

Return the symmetry indicator indices of a symmetry vector `n` as well as its nontrivial
elementary factors, in the context of a set of elementary band representations (EBRs)
provided either as a `BandRepSet`, an integer matrix, or a `Smith` decomposition.

In detail, the method returns the nontrivial indices ``[\\nu_1, \\ldots, \\nu_n]`` and the
associated nontrivial elementary factors ``[\\lambda_1, \\ldots, \\lambda_n]`` of the EBR basis.
The indices ``\\nu_i`` are elements of a cyclic group of order ``\\lambda_i``, i.e. 
``\\nu_i ∈ \\mathbb{Z}_{\\lambda_i} = \\{0, 1, \\ldots, \\lambda_i-1\\}``.

## Implementation

The indices are computed using the Smith normal decomposition ``\\mathbf{B} = \\mathbf{S}
\\boldsymbol{\\Lambda}\\mathbf{T}`` of the EBR matrix ``\\mathbf{B}``. 
Specifically, denoting by ``\\mathbf{s}_i^{-1}`` the ``i``th nontrivial row of
``\\mathbf{S}^{-1}``, the symmetry indicator topological indices of a symmetry vector
``\\mathbf{n}`` are computed as ``\\nu_i = \\mathbf{s}_i^{-1}\\mathbf{n}``.[^HCP]

[^HCP]: [H.C. Po, J. Phys. Cond. Matter **32**, 263001 (2020)](https://doi.org/10.1088/1361-648X/ab7adb).

## Keyword arguments
If `n` is not a compatible band structure (i.e., if `isbandstruct(n, BRS) = false`), an
error is thrown. This behavior can be controlled by two boolean keyword arguments:

- `allow_incompatible` (`false`): if `true`, disables the compatibility check entirely.
- `allow_nonphysical` (`false`): if `true`, allows negative symmetry content, but maintain
  requirement that `n` respects the compatibilty relations in an algebraic sense.
"""
function indicators(
            n::AbstractVector{<:Integer},
            F::Smith;
            allow_incompatible::Bool=false,
            allow_nonphysical::Bool=false)

    if !allow_incompatible || !isbandstruct(n, F; allow_nonphysical)
        _throw_incompatible_or_nonphysical(n)
    end

    idxs = findall(x -> x≠0 && x≠1, F.SNF) # find nontrivial factor groups
    Λ    = @view F.SNF[idxs]               # nontrivial invariant factors
    S̃⁻¹  = @view F.Sinv[idxs, :]           # nontrivial rows of S⁻¹

    return mod.(S̃⁻¹*n, Λ), Λ
end
function indicators(n::AbstractVector{<:Integer}, B::AbstractMatrix{<:Integer}; kws...)
    length(n) == size(B, 1) || throw(DimensionMismatch("incompatible dimensions of `n` and `B`"))
    return indicators(n, smith(B); kws...)
end
function indicators(n::AbstractVector{<:Integer}, BRS::BandRepSet; kws...)
    B = matrix(BRS; includedim=includes_connectivity(n, BRS))
    return indicators(n, B; kws...)
end

# -----------------------------------------------------------------------------------------
# Decomposition in EBRs

@doc """
$(TYPEDSIGNATURES)

Return a decomposition of `n` in the columns of `B` as expansion coefficients `c`. I.e.,
denoting by ``\\mathbf{b}_i`` the columns of `B` and by ``c_i`` the elements of `c`, we find
`c` s.t.:

``\\mathbf{n} = \\sum_i c_i\\mathbf{b}_i``

Depending on the topology of `n`, the coefficients of `c` have the following attributes:
- `TRIVIAL`: `c` is a vector of non-negative integers,
- `FRAGILE`: `c` is a vector of integers,
- `NONTRIVIAL`: `c` is a vector of rationals.

The returned decomposition coefficients `c` are chosen to be the least norm coefficients
(under the constraint that the coefficients belong to the noted integer/rational domains).
In all cases, `c` is returned as a `Vector{Float64}`.

## Keyword arguments
If `n` is not a compatible band structure (i.e., if `isbandstruct(n, BRS) = false`), an
error is thrown. This behavior can be controlled by two boolean keyword arguments:

- `allow_incompatible` (`false`): if `true`, disables the compatibility check entirely.
- `allow_nonphysical` (`false`): if `true`, allows negative symmetry content, but maintain
  requirement that `n` respects the compatibilty relations in an algebraic sense.

"""
function decompose(
            n::AbstractVector{<:Integer},
            B::AbstractMatrix{<:Integer},
            F::Smith=smith(B);
            allow_incompatible::Bool=false,
            allow_nonphysical::Bool=false)
    
    topo = calc_topology(n, F; allow_incompatible, allow_nonphysical) # checks `isbandstruct` as well
    if topo === TRIVIAL
        # then an integer expansion must exist; find out if trivial or fragile
        m = has_posint_expansion(n, B)
        if termination_status(m) == MOI.OPTIMAL # ⇒ trivial
            return value.(m[:c]::Vector{JuMP.VariableRef}) # coefficients
        end
    end

    # `n` is either fragile or nontrivial: in either case, the best we can do is return
    # a possible decomposition using the generalized inverse of B
    Nⁱʳʳ, Nᴱᴮᴿ = size(F.S, 1), size(F.T, 1)
    Λᵍ = diagm(Nⁱʳʳ, Nᴱᴮᴿ, map(λ -> iszero(λ) ? 0.0 : inv(λ), F.SNF))' # TODO: avoid constructing full matrix
    Bᵍ = F.Tinv * Λᵍ * F.Sinv # generalized inverse w/ rational coefficients
    # the general solution is Bᵍn + (I-BᵍB)y w/ y∈ℤᵈ cf. Ben & Israel book, p. 98. Since
    # Bᵍ generates the column space while (I-BᵍB) generates the null space, the two
    # contributions are orthogonal; so the norm of the solution is ‖Aᵍn‖ + ‖(I-AᵍA)y‖.
    # Thus, to get the least norm solution, we simply drop the null space term - this is
    # basically following the usual reasoning for _real_ least norm least squares solutions
    # (see e.g. https://math.stackexchange.com/a/2253614/)
    c = Bᵍ*n 
    
    # check that the coefficients indeed solve the system exactly, not just as a least
    # squares solution
    res = norm(B*c - n)
    res < Crystalline.DEFAULT_ATOL || error("nonzero residual $res of solution c=$c")

    return c
end
function decompose(n::AbstractVector{<:Integer}, BRS::BandRepSet; kws...)
    B = matrix(BRS; includedim=includes_connectivity(n, BRS))
    return decompose(n, B; kws...)
end