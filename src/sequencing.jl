
# simulate the sequencing process based on the config
# return a fastq pair
function pair_end_seq(dna_template, config)

end

function single_end_seq(dna_template, config)
    seq_arr = Char[ nt for nt in dna_template.seq]
    readlen = min(length(seq), config["readlen"])
    rand_err = rand(readlen)
    has_error = false
    for e in rand_err
        if e < config["seq_error_rate"]
            has_error = true
            break
        end
    end
    rand_qual = rand(readlen)
    min_qual = config["normal_base_qual"]["min"]
    max_qual = config["normal_base_qual"]["max"]
    quals = Int(round(rand_qual * min_qual + (1.0 - rand_qual) * max_qual))
    qual_arr = [Char(q+33) for q in quals]
    if !has_error
        sequence = Sequence(ASCIIString(seq_str))
        quality = Quality(ASCIIString(qual_str))
        return sequence, quality
    end
end