# circfull: a tool to detect and quantify full-length circRNA isoforms from circFL-seq

## Introduction

circfull is a circRNA analysis pipeline to detect and quantity full-length circRNA for circFL-seq data

## Availability

circfull is a free software, which can be downloaded from https://github.com/yangence/circfull

## Softwares and packages dependencies
1) Softwares of minimap2 (>=2.1), bedtools (>=2.29.2), samtools (>=1.6) need to be installed by users.
2) porechop(0.2.4) can be downloaded from https://github.com/rrwick/Porechop.
We recommand replace adapters.py in porechop directory with circfull/scripts/bin/adapter.py which contain circFL-seq primers
3) TideHunter (1.0) and trf (4.09) are included in this program.
4) Python 3.x.x and corresponding versions of pysam, numpy, pyfasta, pandas, python-intervals, pyfasta, sklearn, interval, mappy, progressbar, docopt.


## Installation
Install latest release from pip
```
pip install circfull
```

Install latest release from source codes
```
git clone https://github.com/yangence/circfull.git
cd circfull/script
pip install -r requirements.txt
python setup.py install
```

## Required files:

Users can prepare the external files under the following instructions:

1) Indexed genome fasta file

```
samtools faidx $genome
```

2) Tabix indexed gene annotation GTF file

```
grep -v '#' $gtf |sort -k 1,1 -k 4,4n |bgzip >sort.gtf.gz
tabix sort.gtf.gz
```

3) Repetitive elements (optional)

Repetitive elements can be download from https://genome.ucsc.edu/cgi-bin/hgTables 

## Usage
```
Usage: circfull <command> [options]

Command:
    RG             Reference guide detection
    DNSC           De novo self-correction to get consensus sequence
    cRG            Use RG algorithem with sequences corrected by DNSC
    strand         Identify the strand origion of circFL-seq reads
    sRG            RG algorithem with stranded circFL-seq reads.
    mRG            Merge circFL results from RG and cRG/sRG.
    anno           Annotate full-length circRNA.
    geneExp        Quantify gene expression to evaluate linear RNA residual
    var            Calling variants in circRNA.
```


### RG
```
Usage: circfull RG -f fastq -g genome -a anno [-t threads] [-r] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -g genome                   Fasta file of genome.
    -a anno                     Tabix indexed gtf file of gene annotation.
    -t threads                  Number of threads [default: 20].
    -r                          Filter out false discovery.
    -o output                   Output dir [default: circFL_out].

```

### DNSC
```
Usage: circfull DNSC -f fastq  [-t threads] [-c] [-m mem] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -c                          Perform cluster of consensus sequence.
    -m mem                      The maximum memory usage, only worked when -c is designated. [default: 100G]
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
```

### cRG
```
Usage: circfull cRG -f DNSC -g genome -a anno [-t threads] [-u] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f DNSC                     DNSC directory.
    -g genome                   Fasta file of genome.
    -a anno                     Tabix indexed gtf file of gene annotation.
    -t threads                  Number of threads [default: 20].
    -u                          Stranded reads.
    -o output                   Output dir [default: circFL_out].
```

### strand
```
Usage: circfull strand -f fastq -F fastq2 -a anno [-r RG_dir] [-l hang_len] [-t threads] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    Raw circFL-seq fastq file with adapter remain.
    -F fastq2                   Trimmed circFL-seq fastq file.
    -a anno                     Tabix indexed gtf file of gene annotation.
    -r RG_dir                   RG output dir. [default: circFL_out]
    -l hang_len                 Length of sequence to find primers in both end. [default: 100]
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
```

### sRG
```
Usage: circfull sRG -g genome -a anno [-s strandDir] [-t threads] [-r] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -g genome                   Fasta file of genome.
    -a anno                     Tabix indexed gtf file of gene annotation.
    -t threads                  Number of threads [default: 20].
    -s dir                      Strand dir [default: circFL_out].
    -r                          Filter out false discovery.
    -o output                   Output dir [default: circFL_out].
```

### mRG
```
Usage: circfull mRG -f fastq -g genome [-r RG] [-c cRG] [-s sRG] [-t threads] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -g genome                   Fasta file of genome.
    -r RG                       RG directory [default: circFL_out].
    -c cRG                      cRG directory [default: circFL_out].
    -s sRG                      sRG directory [default: circFL_out].
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
```

### anno
```
Usage: circfull anno -b bed -a anno [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -b bed                      Full-length circRNA position information in BED format (12 columns).
    -a anno                     Tabix indexed gtf file of gene annotation.
    -o output                   Output file [default: circ_anno.txt].
```

### geneExp
```
Usage: circfull geneExp -f fastq -r ref [-t threads] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -r ref                      Gene reference
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
```

### var
```
Usage: circfull var -f fastq -g genome [-r RG_dir] [-t threads] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -g genome                   Fasta file of genome.
    -r RG_dir                   RG output dir [default: circFL_out].
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
```




## Output files
| File name         |  Details | 
|   :---            | ---        |
| circFL_Normal.txt       | normal circRNA isoforms detected by `RG`, `cRG`, `sRG` |
| circFL_Normal_pass.txt  | final circRNAs isoforms from `mRG` |
| circFL_Fusion.txt       | fusion circRNA isoforms detected by `RG` |
| novoCluster.txt         | circRNA isoforms detected by `DNSC` with `-c` |
| circ_Var.txt            | variants in circRNAs detected by `var` |
| geneCountDf.txt         | gene expression quantified by `geneExp` |
| strandedFq.fastq        | transcript-strand circFL-seq reads in FASTQ format calcuated by `strand` |
| circ_anno.txt           | full-length circRNA annotated by `anno` |

### circFL_Normal.txt and circFL_Normal_pass.txt
| No. | Column name     |  Details | 
|:---:|   :---          | ---        |
|  1  | circID          | circRNA ID |
|  2  | isoID           | isoform ID |
|  3  | chr             | chromosome |
|  4  | start           | start coordinate of circRNA |
|  5  | end             | end coordinate of circRNA |
|  6  | len             | length of circRNA |
|  7  | exonNum         | number of exons in circRNA |
|  8  | exon_start      | start coordinate of exon |
|  9  | exon_end        | end coordinate of exon |
|  10 | motif           | motif of BS |
|  11 | leftSeq         | motif of start junction |
|  12 | rightSeq        | motif of end junction |
|  13 | exon_leftSeq    | motif of exon start junction |
|  14 | exon_rightSeq   | motif of exon end junction |
|  15 | strand          | strand direction |
|  16 | geneName        | gene name |
|  16 | readCount       | Count of reads supporting this circRNA isoform |
|  17 | readID          | ID of reads supporting this circRNA isoform |

### circFL_Fusion.txt
| No. | Column name     |  Details | 
|:---:|   :---          | ---        |
|  1  | circID          | circRNA ID |
|  2  | isoID           | isoform ID |
|  3  | chr_first       | chromosome of one locus |
|  4  | start_first           | start coordinate of one locus |
|  5  | end_first             | end coordinate of one locus |
|  6  | len_first             | length of one locus |
|  7  | exonNum_first         | number of exons in one locus |
|  8  | exon_start_first      | start coordinate of exon in one locus |
|  9  | exon_end_first        | end coordinate of exon in one locus |
|  10 | exon_leftSeq_first    | motif of exon start junction of one locus |
|  11 | exon_rightSeq_first   | motif of exon end junction of one locus |
|  12 | strand_first          | strand direction of one locus |
|  13 | geneName_first        | gene name of one locus |
|  14 | chr_second            | chromosome of another locus |
|  15 | start_second          | start coordinate of another locus |
|  16 | end_second            | end coordinate of another locus |
|  17 | len_second            | length of another locus |
|  18 | exonNum_second        | number of exons in another locus |
|  19 | exon_start_second     | start coordinate of exon in another locus |
|  20 | exon_end_second       | end coordinate of exon in another locus |
|  21 | exon_leftSeq_second   | motif of exon start junction of another locus |
|  22 | exon_rightSeq_second  | motif of exon end junction of another locus |
|  23 | strand_second         | strand direction of another locus |
|  24 | geneName_second       | gene name of another locus |
|  25 | readCount             | read counts supporting this circRNA isoform |
|  26 | readID                | ID of reads supporting this circRNA isoform |


### novoCluster.txt
| No. | Column name     |  Details | 
|:---:|   :---          | ---        |
|  1  | ID              | read ID |
|  2  | cluster         | cluster ID |
|  3  | strand          | 1/0 denotes sense/antisense of transcript strand, respectively |
|  4  | readLen         | length of read |
|  5  | start           | start position of tandem repeats of read |
|  6  | end             | end position of tandem repeats of read |
|  7  | consLen         | length of consensus sequence |
|  8  | copyNum         | number of tandem repeats |
|  9  | usage           | proportion of tandem repeats |
|  10 | consensus       | consensus sequence |

### circ_Var.txt
| No. | Column name     |  Details | 
|:---:|   :---          | ---        |
|  1  | isoID           | isoform ID |
|  2  | chr             | chromosome |
|  3  | pos             | coordinate of variants |
|  4  | ref             | reference base |
|  5  | alt             | alternative base |
|  6  | refCount        | read counts supporting reference base |
|  7  | altCount        | read counts supporting alternative base |
|  8  | ratio           | read counts ratio of alternative base to total base of the variant |

### geneCountDf.txt
| No. | Column name     |  Details | 
|:---:|   :---          | ---        |
|  1  | geneID          | ensembl gene ID |
|  2  | geneName        | gene name |
|  3  | num             | read counts supporting this gene |
|  4  | ratio           | read counts ratio of this gene to all genes |

### circ_anno.txt
| No. | Column name     |  Details | 
|:---:|   :---          | ---        |
|  1  | chr             | chromosome |
|  2  | start           | start coordinate of circRNA |
|  4  | end             | end coordinate of circRNA |
|  5  | isoID           | isoform ID |
|  6  | strand          | strand direction |
|  7  | exon_start      | start coordinate of exon |
|  8  | exon_end        | end coordinate of exon |
|  9  | len             | length of circRNA |
|  10  | start_type           | annotation for start coordinate of exon |
|  11  | end_type             | annotation for end coordinate of exon |
|  12  | geneName             | gene name |
|  13  | detail_type          | annotation details for isoform |
|  14  | type                 | annotation for isoform |

## Examples
Details in [examples](https://github.com/yangence/circfull/example).
## Copyright and License Information

Copyright (C) 2021 Zelin Liu (zlliu@bjmu.edu.cn). See the [LICENSE](https://github.com/yangence/circfull/blob/master/LICENSE) file for license rights and limitations.
