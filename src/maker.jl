function make_seq(panel_file::AbstractString, profile_file::AbstractString, output_folder::AbstractString = ".")
    profile = JSON.parsefile(profile_file)
    panel, panel_size = load_bed(panel_file)
    if !isdir(output_folder)
        mkpath(output_folder)
    end
    config = profile["config"]
    temp_len_mean = (config["template_len"]["max"] + config["template_len"]["min"])/2
    read_num = (config["depth"] * panel_size) / temp_len_mean
    read_num = Int(floor(read_num))
    rand_pos = rand(read_num)
    for pos in rand_pos
        panel_pos = Int(floor(pos * panel_size))
        chr, pos_in_chr = seek_in_panel(panel, panel_pos)
    end
    #assembly = load_assembly(config["assembly"])
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