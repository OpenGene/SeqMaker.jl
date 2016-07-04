function simulate_copy_number(chr, pos_in_chr, cnv_list)
    for cnv in cnv_list
        if chr == cnv["chrom"]
            # check if it is intersect
            if cnv["start"] < pos_in_chr < cnv["end"]
                return cnv["copy"]
            end
        end
    end
    return 1.0
end