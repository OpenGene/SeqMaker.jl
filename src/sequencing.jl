const SNV_SUBSTITUTION = 0.98
const SNV_INSERTION = 0.01
const SNV_DELETION = 0.01
const per_base_qual_adjust = rand(5, 1000)

# simulate the sequencing process based on the config
# return a fastq pair
function pair_end_seq(dna_template, config)
    const fixed_index = "ATCGATCG"
    index1 = fixed_index
    index2 = fixed_index
    if config["random_index1"]
        index1 = random_index()
    end
    if config["random_index2"]
        index2 = random_index()
    end
    # simulate the duplication rate
    read_num = Int(round((rand() * config["duplication_rate"])/0.5))
    reads = []
    for i in 1:read_num
        r1seq, r1qual = sequence_simulation(dna(dna_template.seq * config["read1_adapter"]), config)
        r2seq, r2qual = sequence_simulation(dna((~dna_template).seq * config["read2_adapter"]), config)
        r1name, r2name = name_simulation(true, index1, index2)
        read1 = FastqRead(r1name, r1seq, "+", r1qual)
        read2 = FastqRead(r2name, r2seq, "+", r2qual)
        pair = FastqPair(read1, read2)
        push!(reads, pair)
    end
    return reads
end

function random_index(index_len = 8)
    atcg = ['A', 'T', 'C', 'G']
    rnd = round(rand(index_len) * 3) + 1
    arr = [atcg[Int(rnd[i])] for i in 1:index_len]
    return join(arr)
end

function name_simulation(pairend, index1, index2)
    const device = "SEQMAKER"
    const run = 1
    const chip = SEQMAKER_VERSION
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

    # simulate quality raise (by coordination) and quality drop
    const qual_curve = [0.95, 1.06, 1.14, 1.16, 1.168,
        1.176, 1.173, 1.172, 1.172, 1.170,
        1.168, 1.165, 1.163,
        1.159, 1.158, 1.153,
        1.148, 1.150, 1.143,
        1.137, 1.133,
        1.127, 1.130, 1.120,
        1.116, 1.112,
        1.105,
        1.096, 1.097, 1.093,
        1.087, 1.083, 1.081,
        1.079, 1.075, 1.073,
        1.066, 1.063,
        1.054, 1.051,
        1.046,
        1.038, 1.03, 1.021, 1.01, 1.02, 1.0]
    const base_index = Dict('A'=>1, 'T'=>2, 'C'=>3, 'G'=>4, 'N'=>5)
    # simulate last base with low quality
    const qual_last_base = 0.9
    qual_adjust_ratio = [0.0 for i in 1:readlen]
    # interpolation with qual_curve
    for i in 1:readlen
        if i == readlen
            qual_adjust_ratio[i] = qual_last_base
            continue
        end
        curve_len = length(qual_curve)
        pos = curve_len*i/readlen
        p = pos - floor(pos)
        left = min(curve_len, Int(floor(pos)) + 1)
        right = min(curve_len, left + 1)
        # linear interpolation
        qual_adjust_ratio[i] = qual_curve[left] * (1-p) + qual_curve[right] * p
    end

    rand_err = rand(readlen)
    has_error = false
    for i in 1:readlen
        if rand_err[i] < config["seq_error_rate"] / qual_adjust_ratio[i]
            has_error = true
            break
        end
    end

    rand_qual = rand(readlen)
    min_qual = config["normal_base_qual"]["min"]
    max_qual = config["normal_base_qual"]["max"]
    quals = rand_qual * min_qual + (1.0 - rand_qual) * max_qual
    # adjust the quality and make phred qual string
    qual_arr = [Char(Int( round(quals[i] * qual_adjust_ratio[i] * (per_base_qual_adjust[base_index[seq_arr[i]], i]*0.02+0.99) ) )+33) for i in 1:readlen]

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
