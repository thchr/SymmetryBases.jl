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
    timereversal::Bool
    compatbasis::Bool
end
function SymBasis(nsᴴ::AbstractMatrix{Int}, brs::BandRepSet, compatbasis::Bool=true)
    kv2ir_idxs = [(f = irlab -> klabel(irlab)==klab; 
                   findfirst(f, brs.irlabs):findlast(f, brs.irlabs)) for klab in brs.klabs]
    return SymBasis(collect(eachcol(nsᴴ)),
                    brs.irlabs, brs.klabs, brs.kvs, kv2ir_idxs, 
                    brs.sgnum, brs.spinful, brs.timereversal, compatbasis)
end

# accessors
parent(sb::SymBasis) = sb.symvecs
num(sb::SymBasis)    = sb.sgnum
irreplabels(sb::SymBasis) = sb.irlabs
klabels(sb::SymBasis)     = sb.klabs
isspinful(sb::SymBasis)   = sb.spinful
fillings(sb::SymBasis)    = [nᴴ[end] for nᴴ in sb.symvecs]

# define the AbstractArray interface for SymBasis
size(sb::SymBasis) = (length(parent(sb)),)
getindex(sb::SymBasis, keys...) = parent(sb)[keys...]
IndexStyle(::SymBasis) = IndexLinear()
