# SeqMaker
Generate simulated sequencing data based on Human Genome Assembly (hg19/hg38)

# Features
* support sequencing error simulation
* support SNV simulation
* support dbSNP (in dev)
* support fusion (in dev)

# Usage
```julia
# clone SeqMaker first
Pkg.clone("https://github.com/sfchen/SeqMaker.jl.git")

using SeqMaker

# it will download human genome assembly data automatically
seqtest("foldername")
```
