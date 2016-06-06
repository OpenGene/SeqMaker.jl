function load_bed(bed_file::AbstractString, assembly)
    io = open(bed_file)
    bed_file = ""
    # TODO: work around for readall missing in master
    if isdefined(Base, :readstring)
        bed_file = readstring(io)
    else
        bed_file = readall(io)
    end
    lines = split(bed_file, '\n')
    panel = []
    start_in_panel = 1
    panel_size = 0
    for line in lines
        line = rstrip(line, '\n')
        cols = split(line)
        if length(cols)<3
            continue
        end
        chr = ASCIIString(cols[1])
        if !haskey(assembly, chr)
            warn("$chr is not in assembly")
            continue
        end
        chrlen = length(assembly[chr])
        from = parse(Int64, ASCIIString(cols[2]))
        to = parse(Int64, ASCIIString(cols[3]))

        # extend from/to by 50*2 bp
        from += sign(from - to) * 50
        to += sign(to - from) * 50

        # clamp from 1 ~ chromo length
        from = max(1, min(chrlen, from))
        to = max(1, min(chrlen, to))

        len = abs(to - from) + 1
        contig_name = "unknown"
        if length(cols)>=4
            contig_name = ASCIIString(cols[4])
        end
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