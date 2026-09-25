using Test

@testset "SymmetryBases" begin
    include("regression.jl")
end

# NB: `full_hilbert_scan.jl` sweeps both Hilbert bases across every space group; it is not
# included here, since the nontopological basis is intractable for some space groups.
