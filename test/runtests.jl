using SeqMaker
using OpenGene
using Base.Test

function test_ngs()
    ngs(depth = 0.01)
    return true
end

# write your own tests here
info("Test starts, this may take minutes to hours if it need to download the genome assembly")
@test test_ngs()
