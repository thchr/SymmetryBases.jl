module SymmetryBases

using Crystalline
using PyCall
using PrettyTables
using JuMP, GLPK
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