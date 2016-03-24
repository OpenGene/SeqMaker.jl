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
make_seq("data/panels/lung_cancer_hg19.bed", "data/profiles/fusion.json", "output")

end # module
