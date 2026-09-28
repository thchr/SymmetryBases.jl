"""
    $(TYPEDEF)
    
with fields: $(TYPEDFIELDS)
"""
struct SymBasis{D} <: AbstractVector{Vector{Int}}
    symvecs::Vector{Vector{Int}}
    irlabs::Vector{String}
    klabs::Vector{String}
    kvs::Vector{KVec{D}}
    kv2ir_idxs::Vector{UnitRange{Int}} # pick k-point; find assoc. ir indices
    sgnum::Int
    spinful::Bool
    timereversal::Bool
    compatbasis::Bool
end
function SymBasis(
    nsᴴ::AbstractMatrix{Int},
    brs::Collection{<:BandRep{D}},
    compatbasis::Bool=true
) where {D}
    irlabs, klabs = irreplabels(brs), klabels(brs)
    kv2ir_idxs = [(f = irlab -> klabel(irlab)==klab;
                   findfirst(f, irlabs):findlast(f, irlabs)) for klab in klabs]
    br = first(brs)
    # NB: materialize the columns, rather than `collect(eachcol(nsᴴ))`: the latter gives a
    #     vector of views, which does not match the field type and keeps `nsᴴ` alive
    symvecs = [Vector{Int}(nᴴ) for nᴴ in eachcol(nsᴴ)]
    sb = SymBasis{D}(
        symvecs, irlabs, klabs, position.(littlegroups(brs)), kv2ir_idxs,
        num(br), isspinful(br), br.timereversal, compatbasis
    )
    return sb
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
