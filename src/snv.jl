function simulate_snv(seq, chr, start, profile_snv)
    seqstr = seq.seq
    for snv in profile_snv
        if chr == snv["chrom"] && start <= snv["pos"] && start+length(seqstr) > snv["pos"]
            seqstr = mutate(seqstr, start, snv)
            # println(seqstr)
        end
    end
    return dna(seqstr)
end

function mutate(seqstr, start, snv)
    # simulate the mutation rate
    if rand() > snv["rate"]
        return seqstr
    end
    offset = snv["pos"] - start
    left = seqstr[1:offset]
    center = ASCIIString(snv["alt"])
    right = seqstr[offset+length(snv["ref"])+1:end]
    seqstr = left * center * right
    return seqstr
end