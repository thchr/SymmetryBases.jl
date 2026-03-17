function show(io::IO, ::MIME"text/plain", sb::SymBasis)
    Nⁱʳʳ = length(sb[1]) - 1

    # print a "title" line and the irrep labels
    println(io, "SymBasis (#", num(sb), "): ",
                length(sb), " Hilbert vectors, sampling ",
                Nⁱʳʳ, " LGIrreps ",
                "(spin-", isspinful(sb) ? "½" : "1", " ",
                sb.timereversal ? "w/" : "w/o", " TR):")

    k_idx = (i) -> findfirst(==(klabel(irreplabels(sb)[i])), klabels(sb)) # highlighters
    h_odd = TextHighlighter((data,i,j) -> i≤Nⁱʳʳ && isodd(k_idx(i)), crayon"light_blue")
    h_ν   = TextHighlighter((data,i,j) -> i==Nⁱʳʳ+1,                 crayon"light_yellow")

    pretty_table(io, 
        # table contents
        stack(sb);
        # column/row names
        column_labels = eachindex(sb),
        row_labels = vcat(sb.irlabs, "μ"),
        # options/formatting/styling
        formatters = [(v,i,j) -> iszero(v) ? "·" : string(v)],
        table_format = TextTableFormat(;
            horizontal_line_at_beginning = true,
            horizontal_line_after_column_labels = true,
            horizontal_line_after_data_rows = true,
            horizontal_lines_at_data_rows = [Nⁱʳʳ],
            @text__no_vertical_lines(),
            vertical_line_after_row_label_column = true
        ),
        row_label_column_alignment = :l,
        alignment = :c, 
        highlighters = [h_odd, h_ν], 
        style = TextTableStyle(column_label = crayon"bold"),
        new_line_at_end = false
    )
end