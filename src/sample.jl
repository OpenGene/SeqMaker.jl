function sample(assembly, chr, pos, temp_len)
    start = max(1, pos - div(temp_len, 2))
    if !haskey(assembly, chr)
        error("The chosen assembly doesn't have $chr")
    end
    start = min(start, length(assembly[chr]) - temp_len)
    seq = assembly[chr].sequence[start:start+temp_len-1]
    return seq, start
end