module SymmetryBases

using Crystalline
using PyCall
using PrettyTables
using JuMP, GLPK
using Nemo # for `calc_topology`
using DocStringExtensions
using LinearAlgebra

import Base: OneTo, show, size, getindex, firstindex, lastindex, IndexStyle, length
import Crystalline: matrix, vecs, num, irreplabels, klabels, isspinful, istimeinvar

# ---------------------------------------------------------------------------------------- #

const PyNormaliz = PyNULL()
function __init__() 
    # import the PyNormaliz library
    # https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
    copy!(PyNormaliz, pyimport("PyNormaliz"))

    # try to stop Nemo spamming its welcome message on package load, following hint from 
    # https://github.com/Nemocas/Nemo.jl/issues/817#issuecomment-613991144, unfortunately
    # this doesn't really seem to work due to the variable not really being available at the
    # right time (... precompilation?), see discussion in 
    # https://discourse.julialang.org/t/redirecting-stdout-to-avoid-banners-while-import-ing/37633/11
    ENV["NEMO_PRINT_BANNER"] = "false"
end

# ---------------------------------------------------------------------------------------- #

include("types.jl")
export SymBasis, fillings
export TopologyKind, TRIVIAL, NONTRIVIAL, FRAGILE

include("show.jl")

include("hilbertbases.jl")
export compatibility_basis, nontopological_basis, split_fragiletrivial

include("symvec.jl")
export has_posint_expansion, calc_detailed_topology, calc_topology,
       isbandstruct
       
# ---------------------------------------------------------------------------------------- #

end # module