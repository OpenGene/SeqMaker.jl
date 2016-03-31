module SeqMaker

using JSON
using OpenGene
using OpenGene.Reference

const SEQMAKER_VERSION = "V1"

# package code goes here

include("panel.jl")
include("sample.jl")
include("snv.jl")
include("fusion.jl")
include("sequencing.jl")
include("maker.jl")

# test
export ngs,
    alk

function ngs(outdir="output", panel_file="", profile_file=""; depth=0)
    if panel_file==""
        panel_file = joinpath(Pkg.dir("SeqMaker"), "data/panels/lung_cancer_hg19.bed")
    end
    if profile_file==""
        profile_file = joinpath(Pkg.dir("SeqMaker"), "data/profiles/example.json")
    end
    make_seq(panel_file, profile_file, outdir, depth)
end

function alk(profile_file=""; d=50)
    return ngs("alk", joinpath(Pkg.dir("SeqMaker"), "data/panels/alk.bed"), profile_file, depth=d)
end

end # module
