const SNV_SUBSTITUTION = 0.9
const SNV_INSERTION = 0.05
const SNV_DELETION = 0.05

# simulate the sequencing process based on the config
# return a fastq pair
function pair_end_seq(dna_template, config)
    r1seq, r1qual = sequence_simulation(dna_template, config)
    r2seq, r2qual = sequence_simulation(~dna_template, config)
    random_index1 = false
    random_index2 = false
    if haskey(config, "random_index1")
        random_index1 = config["random_index1"]
    end
    if haskey(config, "random_index2")
        random_index2 = config["random_index2"]
    end
    r1name, r2name = name_simulation(true, random_index1, random_index2)
    read1 = FastqRead(r1name, r1seq, "+", r1qual)
    read2 = FastqRead(r2name, r2seq, "+", r2qual)
    pair = FastqPair(read1, read2)
    return pair
end

function random_index(index_len = 8)
    atcg = ['A', 'T', 'C', 'G']
    rnd = round(rand(index_len) * 3) + 1
    arr = [atcg[Int(rnd[i])] for i in 1:index_len]
    return join(arr)
end

function name_simulation(pairend = true, random_index1 = false, random_index2 = false)
    const device = "SEQMAKER"
    const run = 1
    const chip = SEQMAKER_VERSION
    const fixed_index = "ATCGATCG"
    index1 = fixed_index
    index2 = fixed_index
    if random_index1
        index1 = random_index()
    end
    if random_index2
        index2 = random_index()
    end
    lane = Int(round(rand() * 3)) + 1
    tile = Int(round(rand() * 99)) + 1
    x = Int(round(rand() * 9999)) + 1
    y = Int(round(rand() * 9999)) + 1

    firstpart = "@$device:$run:$chip:$lane:$tile:$x:$y"
    if pairend
        return "$firstpart 1:N:0:$index1", "$firstpart 2:N:0:$index2"
    else
        return "$firstpart 1:N:0:$index1"
    end
end

function sequence_simulation(dna_template, config)
    seq_arr = Char[ nt for nt in dna_template.seq]
    readlen = min(length(dna_template), config["readlen"])
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
    quals = round(rand_qual * min_qual + (1.0 - rand_qual) * max_qual)
    qual_arr = [Char(Int(q)+33) for q in quals]
    if !has_error
        sequence = Sequence(ASCIIString(seq_arr[1:readlen]))
        quality = Quality(ASCIIString(qual_arr))
        return sequence, quality
    else
        seqlen = 0
        dnapos = 0
        qualpos = 0
        for seq_err in rand_err
            # no error for this base
            if seq_err >= config["seq_error_rate"]
                seqlen += 1
                dnapos += 1
                qualpos += 1
                seq_arr[seqlen] = dna_template.seq[dnapos]
            else
                err = rand()
                # substitution
                if err < SNV_SUBSTITUTION
                    seqlen += 1
                    dnapos += 1
                    qualpos += 1
                    seq_arr[seqlen] = get_other_base(dna_template.seq[dnapos])
                    qual_arr[seqlen] = simulate_error_qual(config["seq_error_qual"]["min"], config["seq_error_qual"]["max"])
                # insertion
                elseif err < SNV_SUBSTITUTION + SNV_INSERTION
                    seqlen += 1
                    seq_arr[seqlen] = get_rand_base()
                    qual_arr[seqlen] = simulate_error_qual(config["seq_error_qual"]["min"], config["seq_error_qual"]["max"])
                # deletion
                elseif err < SNV_SUBSTITUTION + SNV_INSERTION + SNV_DELETION
                    dnapos += 1
                end
            end
            if seqlen>=readlen || dnapos>=length(dna_template)
                break
            end
        end
        sequence = Sequence(ASCIIString(seq_arr[1:seqlen]))
        quality = Quality(ASCIIString(qual_arr[1:seqlen]))
        return sequence, quality
    end
end

function simulate_error_qual(minqual, maxqual)
    rnd = rand()
    q = Int(round(minqual*rnd + maxqual * (1.0 - rnd)))
    return Char(q+33)
end

function get_other_base(base)
    const others = Dict(
        'A'=>['T', 'C', 'G'],
        'T'=>['A', 'C', 'G'],
        'C'=>['T', 'A', 'G'],
        'G'=>['T', 'C', 'A']
        )
    if !haskey(others, uppercase(base))
        return 'N'
    end
    rnd = Int(round(rand() * 2)) + 1
    return others[uppercase(base)][rnd]
end

function get_rand_base()
    const bases = ['A', 'T', 'C', 'G']
    rnd = Int(round(rand() * 3)) + 1
    return bases[rnd]
end