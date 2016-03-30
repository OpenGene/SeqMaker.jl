function simulate_fusion(seq, chr, start, profile_fusion)
    seqstr = seq.seq
    len = length(seq)
    changed = false
    for fusion in profile_fusion
        if fusion_pos_in_left(fusion["left"], chr, start, len)
            if rand() > fusion["rate"]
                return make_fusion_seq(seq, start, fusion["left"], fusion["right"])
            end
        elseif fusion_pos_in_right(fusion["right"], chr, start, len)
            if rand() > fusion["rate"]
                return make_fusion_seq(seq, start, reverse(fusion["right"]), reverse(fusion["left"]))
            end
        end
    end
    return seq
end

function make_fusion_seq(seq, start, original, mate)
    len = length(seq)
    l = 0
    if original["strand"]=="+"
        l = start - original["pos"]
    else
        l = original["pos"] - start
    end
    l = min(len-1, max(1, l))
    return seq
end

function fusion_pos_in_left(fusion, chr, start, len)
    if fusion["strand"] == "+"
        return fusion["chrom"] == chr && fusion["pos"]-len<start && fusion["pos"]<start
    else
        return fusion["chrom"] == chr && fusion["pos"]>start && fusion["pos"]<start+len
    end
end

function fusion_pos_in_right(fusion, chr, start, len)
    if fusion["strand"] == "+"
        return fusion["chrom"] == chr && fusion["pos"]>start && fusion["pos"]<start+len
    else
        return fusion["chrom"] == chr && fusion["pos"]-len<start && fusion["pos"]<start
    end
end

function reverse(fusion)
    strand = "+"
    if fusion["strand"] == "+"
        strand = "-"
    end
    return Dict("chrom"=>fusion["chrom"], "pos"=>fusion["pos"], "strand"=>fusion["strand"])
end