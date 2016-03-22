function make_seq(panel_file::AbstractString, profile_file::AbstractString, output_folder::AbstractString = ".")
    profile = JSON.parsefile(profile_file)
    panel, panel_size = load_bed(panel_file)
    if !isdir(output_folder)
        mkpath(output_folder)
    end
    config = profile["config"]
    temp_len_mean = (config["template_len"]["max"] + config["template_len"]["min"])/2
    read_num = (config["depth"] * panel_size) / temp_len_mean
    assembly = load_assembly(config["assembly"])
end