# Regression tests for the Hilbert bases, and a check that spinful input works.
#
# The reference values below were obtained from Bilbao's tabulated EBRs, i.e. from what
# `bandreps` returned before Crystalline v0.7. An exhaustive A/B over all 230 space groups
# found the compatibility basis to be *identical* to the one obtained from the EBRs that
# Crystalline now computes itself, once the irrep rows are matched by label — as it must
# be, since the basis depends only on the cone spanned by the EBRs and not on which EBRs
# span it. These spot-checks guard that agreement against future regression.
using Crystalline, SymmetryBases, Test

# sgnum => (number of irreps, number of Hilbert basis vectors, Σ fillings, extremal filling)
const COMPATIBILITY_BASIS_REFS = Dict(
    1   => (8,     1,   1, (1, 1)),
    2   => (16,  256, 256, (1, 1)),
    13  => (20,   36,  72, (2, 2)),
    68  => (22,   18,  64, (2, 4)),
    81  => (16,   89, 146, (1, 2)),
    147 => (16,  169, 322, (1, 2)),
    174 => (36,  192, 612, (1, 4)),
)

@testset "`compatibility_basis` against pre-v0.7 reference values" begin
    for (sgnum, (Nⁱʳʳ, Nᴴ, Σμ, μ_extrema)) in sort!(collect(COMPATIBILITY_BASIS_REFS))
        sb, brs = compatibility_basis(sgnum, 3)
        @test length(irreplabels(sb)) == Nⁱʳʳ
        @test length(sb) == Nᴴ
        @test sum(fillings(sb)) == Σμ
        @test extrema(fillings(sb)) == μ_extrema

        # every Hilbert basis vector must be a non-negative, compatible symmetry vector
        @test all(nᴴ -> iscompatible(nᴴ, brs), sb)
    end
end

@testset "spinful input" begin
    sgnum = 13
    sb, brs = compatibility_basis(sgnum, 3; spinful = Val(true))
    @test isspinful(sb)
    @test isspinful(first(brs))
    @test length(sb) == 81
    @test all(nᴴ -> iscompatible(nᴴ, brs), sb)

    # the spinless basis of the same space group differs, and is the default
    sb_spinless, _ = compatibility_basis(sgnum, 3)
    @test !isspinful(sb_spinless)
    @test length(sb_spinless) ≠ length(sb)

    # a plain `Bool` is accepted as well, though it is not type-stable
    sb′, _ = compatibility_basis(sgnum, 3; spinful = true)
    @test collect(sb′) == collect(sb)
end
