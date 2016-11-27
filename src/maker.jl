function make_seq(panel_file::AbstractString, profile_file::AbstractString, output_folder::AbstractString, depth)
    profile = JSON.parsefile(profile_file)
    config = profile["config"]
    fill_default_config!(config)

    # load assembly
    assembly = load_assembly(config["assembly"])
    panel, panel_size = load_bed(panel_file, assembly)

    if !isdir(output_folder)
        mkpath(output_folder)
    end

    if depth<=0
        depth = config["depth"]
    end

    io = nothing
    filename = joinpath(output_folder, "SEQMAKER-$SEQMAKER_VERSION")
    filename *= "-" * basename(panel_file)
    filename *= "-" * basename(profile_file)

    if config["pair-end"]
        io = fastq_open_pair("$filename.R1.fq", "$filename.R2.fq", "w")
    else
        io = fastq_open("$filename.R1.fq", "w")
    end

    temp_max = config["template_len"]["max"]
    temp_min = config["template_len"]["min"]
    temp_len_mean = ( temp_min + temp_max )/2
    read_num = (depth * panel_size) / (temp_len_mean * config["duplication_rate"])
    read_num = Int(floor(read_num))
    rand_pos = rand(read_num)
    rand_temp_len = rand(read_num)
    for i in 1:length(rand_pos)
        panel_pos = Int(floor(rand_pos[i] * panel_size))
        # linear interpolation
        temp_len = Int(floor(rand_temp_len[i] * temp_min + (1.0-rand_temp_len[i]) * temp_max))
        chr, pos_in_chr = seek_in_panel(panel, panel_pos)
        cycles = 1
        if haskey(profile, "cnv")
            copy_num = simulate_copy_number(chr, pos_in_chr, profile["cnv"])
            cycles = Int(round((rand() * copy_num)/0.5))
        end
        for c in 1:cycles
            # vibration
            this_pos_in_chr = pos_in_chr + Int(round((rand() - 0.5) * 0.4 * temp_len))
            this_temp_len = temp_len + Int(round((rand() - 0.5) * 0.1 * temp_len))
            # sampling
            original_seq, start = sample(assembly, chr, this_pos_in_chr , this_temp_len)
            seq = simulate_snv(original_seq, chr, start, profile["snv"])
            seq = simulate_fusion(seq, chr, start, profile["fusion"], assembly)
            # simulate watson/crick strand
            if rand()> 0.5
                seq = ~seq
            end
            if config["pair-end"]
                reads = pair_end_seq(seq, config)
                for pair in reads
                    fastq_write_pair(io, pair)
                end
            end
        end
    end
    close(io)
end

function seek_in_panel(panel, panel_pos)
    # binary search to find the contig
    left = 1
    right = length(panel)
    contig = 0
    while left < right
        center = div(left+right, 2)
        rec = panel[center]
        if rec["start_in_panel"]<=panel_pos && rec["end_in_panel"]>=panel_pos
            contig = center
            break
        elseif rec["start_in_panel"]>panel_pos
            right = center - 1
        elseif rec["end_in_panel"]<panel_pos
            left = center + 1
        end
    end
    if contig == 0
        contig = left
    end
    rec = panel[contig]
    offset = panel_pos - rec["start_in_panel"]
    pos_in_chr = rec["from"] + offset
    return rec["chr"], pos_in_chr
end

function fill_default_config!(config)
    if !haskey(config, "duplication_rate")
        config["duplication_rate"] = 1.0
    end
    if !haskey(config, "random_index1")
        config["random_index1"] = false
    end
    if !haskey(config, "random_index2")
        config["random_index2"] = false
    end
    if !haskey(config, "read1_adapter")
        config["read1_adapter"] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    end
    if !haskey(config, "read2_adapter")
        config["read2_adapter"] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    end
end