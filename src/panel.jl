function load_bed(bed_file::AbstractString)
    io = open(bed_file)
    bed_file = readall(io)
    lines = split(bed_file, '\n')
    panel = []
    start_in_panel = 1
    panel_size = 0
    for line in lines
        line = rstrip(line, '\n')
        cols = split(line)
        if length(cols)<4
            continue
        end
        chr = ASCIIString(cols[1])
        from = parse(Int64, ASCIIString(cols[2]))
        to = parse(Int64, ASCIIString(cols[3]))
        len = abs(to - from) + 1
        contig_name = ASCIIString(cols[4])
        panel_size = start_in_panel+len-1
        rec = Dict(
            "chr"=>chr,
            "name"=>contig_name,
            "from"=>from, "to"=>to,
            "start_in_panel"=>start_in_panel,
            "end_in_panel"=> panel_size
            )
        push!(panel, rec)
        start_in_panel += len
    end
    return panel, panel_size
end