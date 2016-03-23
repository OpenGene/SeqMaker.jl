function make_seq(panel_file::AbstractString, profile_file::AbstractString, output_folder::AbstractString = ".")
    profile = JSON.parsefile(profile_file)
    panel, panel_size = load_bed(panel_file)
    if !isdir(output_folder)
        mkpath(output_folder)
    end

    config = profile["config"]

    # load assembly
    assembly = load_assembly(config["assembly"])

    temp_max = config["template_len"]["max"]
    temp_min = config["template_len"]["min"]
    temp_len_mean = ( temp_min + temp_max )/2
    read_num = (config["depth"] * panel_size) / temp_len_mean
    read_num = Int(floor(read_num))
    rand_pos = rand(read_num)
    rand_temp_len = rand(read_num)
    for i in 1:length(rand_pos)
        panel_pos = Int(floor(rand_pos[i] * panel_size))
        # linear interpolation
        temp_len = Int(floor(rand_temp_len[i] * temp_min + (1.0-rand_temp_len[i]) * temp_max))
        chr, pos_in_chr = seek_in_panel(panel, panel_pos)
        seq, start = sample(assembly, chr, pos_in_chr, temp_len)
        seq = simulate_snv(seq, chr, start, profile["snv"])
        # simulate watson/crick strand
        if rand()> 0.5
            seq = ~seq
        end
    end
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