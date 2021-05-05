# PURPOSE: log how long it takes to calculate the Hilbert bases for {BS} and {AI+F}

using SymmetryBases

has_tr = true

t_sb   = Vector{Float64}(undef, 230)
t_ntsb = Vector{Float64}(undef, 230)
N_sb   = Vector{Int}(undef, 230)
N_ntsb = Vector{Int}(undef, 230)
for sgnum in 1:230
    println(sgnum)
    t = @elapsed begin
        sb, _ = compatibility_basis(sgnum, 3, timereversal=has_tr)    # basis for {BS}
    end
    t_sb[sgnum] = t
    N_sb[sgnum] = length(sb)
    println("   ", "sb:   ", round(t, digits=2), " s")
    
    t′ = @elapsed begin
        ntsb, _ = nontopological_basis(sgnum, 3, timereversal=has_tr) # basis for {AI+F}
    end
    t_ntsb[sgnum] = t′
    N_ntsb[sgnum] = length(ntsb)
    println("   ", "ntsb: ", round(t′, digits=2), " s")
end

#=
using PlotlyJS

p_sb   = plot(x = 1:230, y = t_sb,   type="bar")
p_ntsb = plot(x = 1:230, y = t_ntsb, type="bar")
=#