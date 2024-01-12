function show(io::IO, ::MIME"text/plain", sb::SymBasis)
    Nⁱʳʳ = length(sb[1]) - 1

    # print a "title" line and the irrep labels
    println(io, "SymBasis (#", num(sb), "): ",
                length(sb), " Hilbert vectors, sampling ",
                Nⁱʳʳ, " LGIrreps ",
                "(spin-", isspinful(sb) ? "½" : "1", " ",
                sb.timereversal ? "w/" : "w/o", " TR)")

    k_idx = (i) -> findfirst(==(klabel(irreplabels(sb)[i])), klabels(sb)) # highlighters
    h_odd = Highlighter((data,i,j) -> i≤Nⁱʳʳ && isodd(k_idx(i)), crayon"light_blue")
    h_ν   = Highlighter((data,i,j) -> i==Nⁱʳʳ+1,                 crayon"light_yellow")

    pretty_table(io, 
        # table contents
        matrix(sb);
        # column/row names
        header = eachindex(sb),
        row_labels = vcat(sb.irlabs, "μ"),
        # options/formatting/styling
        formatters = (v,i,j) -> iszero(v) ? "·" : string(v),
        vlines = [1,],
        hlines = [:begin, 1, Nⁱʳʳ+1, :end],
        row_label_alignment = :l,
        alignment = :c, 
        highlighters = (h_odd, h_ν), 
        header_crayon = crayon"bold"
        )

    # print k-vec labels
    print(io, "  KVecs: ")
    join(io, klabels(sb), ", ")
end