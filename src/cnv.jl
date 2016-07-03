function simulate_copy_number(chr, start, seqlen, cnv_list)
    for cnv in cnv_list
        if chr == cnv["chrom"]
            # check if it is intersect
            overlap = intersect(cnv["start"]:cnv["end"], start:start+seqlen)
            if length(overlap) > seqlen / 3
                return cnv["copy"]
            end
        end
    end
    return 1.0
end