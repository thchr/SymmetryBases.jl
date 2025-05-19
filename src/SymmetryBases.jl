module SymmetryBases

using Crystalline
using PythonCall: pynew, pyimport, pycopy!, pyconvert
using PrettyTables
using JuMP
using GLPK     # for `has_posint_expansion`
using HiGHS    # ┌ 
using Pajarito # │ for `decompose` with `ensure_min_norm = true`
using Hypatia  # └
using DocStringExtensions
using LinearAlgebra

import Base: OneTo, show, size, getindex, firstindex, lastindex, IndexStyle, length, parent
import Crystalline: num, irreplabels, klabels, isspinful

# ---------------------------------------------------------------------------------------- #

const PyNormaliz = pynew()
function __init__() 
    # import the PyNormaliz library
    try
        pycopy!(PyNormaliz, pyimport("PyNormaliz"))
    catch
        @warn "PyNormaliz could not be imported: related functionality of SymmetryBases.jl is nonfunctional"
    end
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
       isbandstruct, indicators, decompose
       
# ---------------------------------------------------------------------------------------- #

end # module