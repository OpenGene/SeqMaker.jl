# SeqMaker
Generate simulated sequencing data based on Human Genome Assembly (hg19/hg38)

# Features
* support sequencing error simulation
* support SNV simulation
* support dbSNP (in dev)
* support fusion (in dev)

# Panel
* A file describes the regions of your target capturing
* In bed format, each line is a record
* If you don't do capturing, you can use a whole genome panel

# Profile
* A file describes the sequencing simulation config and the mutation to simulate
* In json format

# Usage
```julia
# clone SeqMaker first
Pkg.clone("https://github.com/sfchen/SeqMaker.jl.git")

using SeqMaker

# it will download human genome assembly data automatically, please make sure your system can access internet
# ngs(outdir, panel_file, profile_file; depth=0)
# outdir is default to "output"
# panel_file is default to "data/panels/lung_cancer_hg19.bed"
# profile is default to "data/profiles/example.json"

# use all default settings, in this case, depth must be set in the profile.json
ngs()

# set the depth to 100X
ngs(depth=100)

# set the output folder and depth
# depth can be float, like 0.1
ngs("myout", depth=100)

# set the output folder and panel
ngs("myout", "panelfile.bed")

# set the output folder and profile
ngs("myout", "", "profile.json")

# set the output folder, panel and profile
ngs("myout", "panelfile.bed", "profile.json")

# set all parameters
ngs("myout", "panelfile.bed", "profile.json", depth=100)
```
