"""
    $(TYPEDEF)
    
with fields: $(TYPEDFIELDS)
"""
struct SymBasis <: AbstractVector{Vector{Int}}
    symvecs::Vector{Vector{Int}}
    irlabs::Vector{String}
    klabs::Vector{String}
    kvs::Vector{<:KVec} # TODO: make parameteric
    kv2ir_idxs::Vector{UnitRange{Int}} # pick k-point; find assoc. ir indices
    sgnum::Int
    spinful::Bool
    timeinvar::Bool
    compatbasis::Bool
end
function SymBasis(nsᴴ::AbstractMatrix{Int}, BRS::BandRepSet, compatbasis::Bool=true)
    kv2ir_idxs = [(f = irlab -> klabel(irlab)==klab; 
                   findfirst(f, BRS.irlabs):findlast(f, BRS.irlabs)) for klab in BRS.klabs]
    return SymBasis(collect(eachcol(nsᴴ)),
                    BRS.irlabs, BRS.klabs, BRS.kvs, kv2ir_idxs, 
                    BRS.sgnum, BRS.spinful, BRS.timeinvar, compatbasis)
end

# accessors
matrix(sb::SymBasis) = hcat(sb.symvecs...)
vecs(sb::SymBasis)   = sb.symvecs
num(sb::SymBasis)    = sb.sgnum
irreplabels(sb::SymBasis) = sb.irlabs
klabels(sb::SymBasis)     = sb.klabs
isspinful(sb::SymBasis)   = sb.spinful
istimeinvar(sb::SymBasis) = sb.timeinvar
fillings(sb::SymBasis)    = [nᴴ[end] for nᴴ in sb.symvecs]

# define the AbstractArray interface for SymBasis
size(sb::SymBasis) = (length(vecs(sb)),)
getindex(sb::SymBasis, keys...) = vecs(sb)[keys...]
IndexStyle(::SymBasis) = IndexLinear()

# ---------------------------------------------------------------------------------------- #

@enum TopologyKind begin
    TRIVIAL    = 0
    NONTRIVIAL = 1
    FRAGILE    = 2
end