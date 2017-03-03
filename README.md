# SeqMaker

Generate simulated sequencing data based on Human Genome Assembly (hg19/hg38), written in [Julia](http://julialang.org/) language  
This tool can be used to test or benchmark bioinformatics algorithms, software or pipelines

# Features
* support sequencing error simulation
* support SNV simulation
* support fusion simulation
* support duplication simulation
* support single molecule indexing
* support CNV

## Julia
Julia is a fresh programming language with `C/C++` like performance and `Python` like simple usage  
On Ubuntu, you can install Julia by `sudo apt-get install julia`, and type `julia` to open Julia interactive prompt

# Examples
```julia
# launch Julia

# clone SeqMaker first
Pkg.clone("https://github.com/OpenGene/SeqMaker.jl.git")

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

# simulate whole genome sequencing with a depth of 10x
wgs("myout",  depth=10)
```

# Panel
* A file describes the regions of your target capturing
* In bed format, each line is a record
* If you don't do capturing, you can use a whole genome panel

```tsv
chr9    133588266   133763062   ABL1
chr14   105235686   105262088   AKT1
chr19   40736224    40791443    AKT2
chr2    29415640    30144432    ALK
chrX    66764465    66950461    AR
chr11   108093211   108239829   ATM
chr3    142168077   142297668   ATR
chr2    111876955   111926024   BCL2L11
chr7    140419127   140624564   BRAF
chr17   41196312    41277500    BRCA1
chr13   32889611    32973805    BRCA2
chr11   69455855    69469242    CCND1
chr12   58141510    58149796    CDK4
chr7    92234235    92465908    CDK6
chr5    149432854   149492935   CSF1R
chr1    162601163   162757190   DDR2
```

# Profile
* A file describes the sequencing simulation config and the mutation to simulate
* In json format

```json
{
    "config":{
        "depth":300,
        "pair-end":true,
        "readlen":151,
        "assembly":"hg19",
        "seq_error_rate":0.001,
        "duplication_rate":3.0,
        "random_index1":false,
        "random_index2":true,
        "read1_adapter":"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "read2_adapter":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "template_len":{
            "min":130,
            "max":230
        },
        "normal_base_qual":{
            "min":30,
            "max":37
        },
        "seq_error_qual":{
            "min":8,
            "max":20
        }
    },
    "fusion":[
        {
            "name":"ALK-intron19-EML4-intron13",
            "left":
            {
                "chrom":"chr2",
                "strand":"-",
                "pos":29447873
            },
            "right":
            {
                "chrom":"chr2",
                "strand":"+",
                "pos":42526793
            },
            "rate":0.1
        }
    ],
    "snv":[
        {
            "name":"EGFR-L861Q",
            "chrom":"chr7",
            "pos":55259524,
            "ref":"T",
            "alt":"A",
            "rate":0.2
        },
        {
            "name":"KRAS-G12D",
            "chrom":"chr12",
            "pos":25398284,
            "ref":"G",
            "alt":"A",
            "rate":0.75
        }
    ],
    "cnv":[
        {
            "name":"MET-Amplification",
            "chrom":"chr7",
            "start":116312406,
            "end":116438440,
            "copy":3.0
        }
    ]
}
```

# Cite SeqMaker
If you use SeqMaker for your research, you can cite this tool as:
```
Chen, S., Han, Y., Guo, L., Hu, J., & Gu, J. (2016, December). SeqMaker: A next generation sequencing simulator with variations, sequencing errors and amplification bias integrated. In Bioinformatics and Biomedicine (BIBM), 2016 IEEE International Conference on (pp. 835-840). IEEE.
```
