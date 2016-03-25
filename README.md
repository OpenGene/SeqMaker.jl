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

# it will download human genome assembly data automatically, please make sure your system can access internet
# ngs(outdir, panel_file, profile_file)
# outdir is default to "output"
# panel_file is default to "data/panels/lung_cancer_hg19.bed"
# profile is default to "data/profiles/example.json"
ngs()
```
