# ---------------------------------------------------------------------------------------- #
#                  Methods to test the topology of a given symmetry vector n               #
# ---------------------------------------------------------------------------------------- #

using Crystalline: includes_connectivity

"""
    $(TYPEDSIGNATURES)

Check whether a positive-integer coefficient expansion exists for `n` in the basis of the 
columns of `M`, i.e. whether there exists ``cᵢ∈ℕ`` (``=0,1,2,...``)  such that ``Mc=n``.
"""
function has_posint_expansion(n::AbstractVector{<:Integer}, M::AbstractMatrix{<:Integer})
    N = size(M, 2)
    # feasibility problem with positive-integer variables, subject to the condition Mc = n
    m = Model(GLPK.Optimizer) # TODO: can swap for HiGHS?
    @variable(m, c[1:N] >= 0, Int)
    @constraint(m, M*c .== n)
    # try to solve the model
    optimize!(m)

    return m
end

"""
    calc_detailed_topology(n, B::Matrix{<:Integer}, [F::Smith=smith(B)]; allow_negative)
    calc_detailed_topology(n, brs::BandRepSet; allow_negative)
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
the compatibility basis `sb` (⇒ nontrivial), in `sb` *and* in the nontopological basis
but *not* in the EBR basis (⇒ fragile), or in `sb` *and* in the EBR basis (⇒ trivial) -
but is much faster (due to not solving multiple non-negative feasibility problems).

## Keyword arguments
If `n` is not a compatible band structure (i.e., if `iscompatible(n, brs) = false`), an
error is thrown. This behavior can be controlled by two boolean keyword arguments:

- `allow_incompatible` (`false`): if `true`, disables the compatibility check entirely.
- `allow_negative` (`false`): if `true`, allows negative symmetry content, but maintain
  requirement that `n` respects the compatibilty relations in an algebraic sense.
"""
function calc_detailed_topology(
            n::AbstractVector{<:Integer},
            B::AbstractMatrix{<:Integer},
            F::Smith=smith(B);
            allow_incompatible::Bool=false,
            allow_negative::Bool=false)

    coarse_topo = calc_topology(n, F; allow_incompatible, allow_negative) # checks `iscompatible` as well
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

function calc_detailed_topology(
    n::AbstractVector{<:Integer},
    brs::Union{Collection{<:NewBandRep}, BandRepSet};
    kws...
)
    B = stack(brs)
    if !includes_connectivity(n, brs)
        return calc_detailed_topology(n, (@view B[1:end-1, :]); kws...)
    end

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

    brs = bandreps(sgnum, D; spinful, timereversal, allpaths)
    return calc_detailed_topology(n, brs; kws...)
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

`c` is returned as a `Vector{Rational{Int}}`.

## Keyword arguments
If `n` is not a compatible band structure (i.e., if `iscompatible(n, brs) = false`), an
error is thrown. This behavior can be controlled by two boolean keyword arguments:

- `allow_incompatible` (`false`): if `true`, disables the compatibility check entirely.
- `allow_negative` (`false`): if `true`, allows negative symmetry content, but maintain
  requirement that `n` respects the compatibilty relations in an algebraic sense.
- `seek_minimal_norm` (`true`): if `true`, additionally solves an integer quadratic problem
  to explicitly minimize the norm of the decomposition vector. Useful for ensuring a
  maximally simple expansion. Typically, however, this is achieved even with
  `seek_minimal_norm = false`.
  Setting to `false` will improve performance, usually substantially.
"""
function decompose(
            n::AbstractVector{<:Integer},
            B::AbstractMatrix{<:Integer},
            F::Smith=smith(B);
            allow_incompatible::Bool=false,
            allow_negative::Bool=false,
            seek_minimal_norm::Bool=true)
    
    c::Union{Vector{Float64}, Nothing} = nothing # placeholder, for branching
    topo = calc_topology(n, F; allow_incompatible, allow_negative) # checks `iscompatible` as well
    if topo === TRIVIAL
        # then an integer expansion must exist; find out if trivial or fragile
        m = has_posint_expansion(n, B)
        if termination_status(m) == MOI.OPTIMAL # ⇒ trivial
            c = value.(m[:c]::Vector{JuMP.VariableRef})::Vector{Float64} # coefficients
        end
    end

    if seek_minimal_norm || isnothing(c)
        # `n` is either fragile or nontrivial: in either case, the best we can do is return
        # a possible decomposition using the generalized inverse of B
        Nⁱʳʳ, Nᴱᴮᴿ = size(F.S, 1), size(F.T, 1)
        Λᵍ = diagm(Nⁱʳʳ, Nᴱᴮᴿ, map(λ -> iszero(λ) ? 0.0 : inv(λ), F.SNF))' # TODO: avoid constructing full matrix
        Bᵍ = F.Tinv * Λᵍ * F.Sinv # generalized inverse w/ rational coefficients
        if isnothing(c)
            # situation: there is no positive, integer-coefficient solution; must have
            # negative or rational coefficients - we then use the following general solution
            # procedure: namely, c = Bᵍn + (I-BᵍB)y w/ y∈ℤᵈ cf. Ben & Israel book, p. 98.
            # Here, Bᵍ generates the column space while (I-BᵍB) generates the null space,
            # but unlike the real-coefficient solution case (https://math.stackexchange.com/a/2253614/),
            # this apparently does not necessitate that the two spaces are orthogonal;
            # instead, a smaller-norm solution might be found by varying y.
            c = Bᵍ*n
        end
        if seek_minimal_norm
            c = _minimize_decomposition_norm(c::Vector{Float64}, B, Bᵍ)
        end
    end

    # cast to a rational vector
    c′ = rationalize.(c; tol=1e-4)

    # check that the coefficients indeed solve the system exactly, not just as a least
    # squares solution
    res = norm(B*c′ - n)
    res < Crystalline.DEFAULT_ATOL || error("nonzero residual $res of solution c′=$c′")

    return c′
end
function decompose(
    n::AbstractVector{<:Integer},
    brs::Union{Collection{<:NewBandRep}, BandRepSet};
    kws...
)
    B = stack(brs)
    if !includes_connectivity(n, brs)
        return decompose(n, (@view B[1:end-1, :]); kws...)
    end
    return decompose(n, B; kws...)
end

function _minimize_decomposition_norm(c, B, Bᵍ)
    # all solutions have the form `c + N*y` where y is an integer vector
    N = round.(Int, I-Bᵍ*B) # null-space matrix
    N ≈ I-Bᵍ*B || error("failed to compute integer-valued null-space matrix")

    idxs = findall(!iszero, eachcol(N)) 
    N′ = N[:,idxs] # nonzero columns of null-space

    # set up & solve norm-minimization problem: min_y ‖c + Ny‖ for yᵢ∈ℤ
    m = Model( 
        # copied whole-sale from from Pajarito's README.md
        optimizer_with_attributes(
            Pajarito.Optimizer,
            "oa_solver" => optimizer_with_attributes(
                HiGHS.Optimizer,
                MOI.Silent() => true,
                "mip_feasibility_tolerance" => 1e-8,
                "mip_rel_gap" => 1e-6),
            "conic_solver" => optimizer_with_attributes(
                Hypatia.Optimizer,
                MOI.Silent() => true),
            )
    )
    set_attribute(m, "verbose", false)
    @variable(m, y[1:length(idxs)], Int)
    @objective(m, Min, sum(abs2, c+N′*y))
    optimize!(m)

    # obtain new coefficient solution as `c+N*y`
    y_optim = round.(Int, value.(m[:y]::Vector{JuMP.VariableRef}))
    Δc = N′*y_optim

    return c + Δc
end