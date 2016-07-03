function simulate_fusion(seq, chr, start, profile_fusion, assembly)
    seqstr = seq.seq
    len = length(seq)
    changed = false
    l=0
    r=0
    for fusion in profile_fusion
        on_fusion = false
        if is_in_fusion(fusion["left"], chr, start, len)
            on_fusion = true
            if fusion["left"]["strand"]=="+"
                l = fusion["left"]["pos"] - start
            else
                l = start - fusion["left"]["pos"]
            end
            l = min(len-1, max(1, l))
            r = len - l
        elseif is_in_fusion(fusion["right"], chr, start, len)
            on_fusion = true
            if fusion["right"]["strand"]=="+"
                r = (start+len) - fusion["right"]["pos"]
            else
                r = fusion["right"]["pos"] - (start-len)
            end
            r = min(len-1, max(1, r))
            l = len - r
        end
        if on_fusion && rand() < fusion["rate"]
            return make_fusion_seq(l, r, fusion["left"], fusion["right"], assembly)
        end
    end
    return seq
end

function make_fusion_seq(l, r, fusion_left, fusion_right, assembly)
    leftseq = nothing
    rightseq = nothing
    if fusion_left["strand"]=="+"
        leftseq = assembly[fusion_left["chrom"]].sequence[fusion_left["pos"]-l-1:fusion_left["pos"]]
    else
        leftseq = assembly[fusion_left["chrom"]].sequence[fusion_left["pos"]:fusion_left["pos"]+l-1]
        leftseq = ~leftseq
    end
    if fusion_right["strand"]=="+"
        rightseq = assembly[fusion_right["chrom"]].sequence[fusion_right["pos"]:fusion_right["pos"]+r-1]
    else
        rightseq = assembly[fusion_right["chrom"]].sequence[fusion_right["pos"]-r-1:fusion_right["pos"]]
        rightseq = ~rightseq
    end
    # connect left/right seq to a fusion seq
    fusion_seq = dna(uppercase(leftseq.seq * rightseq.seq))
    # println(">fusion_left:", l, "_right:", r)
    # println(fusion_seq.seq)
    return fusion_seq
end

function is_in_fusion(fusion, chr, start, len)
    if fusion["strand"] == "+"
        return fusion["chrom"] == chr && start < fusion["pos"] && fusion["pos"] < start + len
    else
        return fusion["chrom"] == chr && start-len < fusion["pos"] && fusion["pos"] < start
    end
end


function reverse(fusion)
    strand = "+"
    if fusion["strand"] == "+"
        strand = "-"
    end
    return Dict("chrom"=>fusion["chrom"], "pos"=>fusion["pos"], "strand"=>strand)
end