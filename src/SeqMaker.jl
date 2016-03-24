module SeqMaker

using JSON
using OpenGene
using OpenGene.Reference

const SEQMAKER_VERSION = "V1"

# package code goes here

include("panel.jl")
include("sample.jl")
include("snv.jl")
include("sequencing.jl")
include("maker.jl")

# test
export seqtest
function seqtest(outdir="output")
    panel_file = joinpath(Pkg.dir("SeqMaker"), "data/panels/lung_cancer_hg19.bed")
    profile_file = joinpath(Pkg.dir("SeqMaker"), "data/profiles/fusion.json")
    make_seq(panel_file, profile_file, outdir)
end

end # module
