module SeqMaker

using JSON
using OpenGene
using OpenGene.Reference

const SEQMAKER_VERSION = "V1"

# package code goes here

include("compat.jl")
include("panel.jl")
include("sample.jl")
include("snv.jl")
include("cnv.jl")
include("fusion.jl")
include("sequencing.jl")
include("maker.jl")

# test
export ngs,
    alk,
    wgs

function ngs(outdir="output", panel_file="", profile_file=""; depth=0)
    if panel_file==""
        panel_file = joinpath(dirname(@__FILE__), "../data/panels/lung_cancer_hg19.bed")
    end
    if profile_file==""
        profile_file = joinpath(dirname(@__FILE__), "../data/profiles/example.json")
    end
    make_seq(panel_file, profile_file, outdir, depth)
end

function alk(profile_file=""; depth=50)
    if profile_file==""
        profile_file = joinpath(dirname(@__FILE__), "../data/profiles/fusion.json")
    end
    return ngs("alk", joinpath(dirname(@__FILE__), "../data/panels/alk.bed"), profile_file, depth=depth)
end

function wgs(outdir="wgs", profile_file=""; depth=30)
    return ngs(outdir, joinpath(dirname(@__FILE__), "../data/panels/hg19_whole_genome.bed"), profile_file, depth=depth)
end

end # module
