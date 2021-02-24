## identify full-length circRNA

### required files
```
rawFq=M1.example.fq # included in example directory
gtfFile=gencode.v19.annotation.sort.gtf.gz # sorted and tabix indexed
genome=hg19.fa
```
### chop and split raw reads to get clean reads
```
cleanFq=M1.clean.fq
thread=80
porechop-runner.py -t $thread -i $rawFq -o $cleanFq --barcode_threshold 95 --check_reads 1000
```

### circfull with reference guide method to identify full-length circRNA from clean reads
```
outDir=circFL_example
circfull  RG  -f $cleanFq  -g $genome -a $gtfFile  -t $thread -o $outDir # output circFL_Normal.txt and circFL_Fusion.txt in $outDir/RG
```
Notices: More threads require more memory

### circfull with strand module to get transcript-strand reads
```
circfull strand -f $rawFq  -F  $cleanFq   -t $thread -a $gtfFile -l 100 -o $outDir -r $outDir # output strandedFq.fastq in $outDir/strand
```

### circfull with reference guide method to identify full-length circRNA from stranded reads
```
circfull sRG -g $genome -a $gtfFile  -t $thread -s $outDir  -o $outDir # output circFL_Normal.txt in $outDir/sRG
```

### circfull with de novo self correction method to get consensus sequence from clean reads
```
circfull DNSC -f $cleanFq -t $thread -o $outDir # get consensus sequnece
circfull DNSC -f $cleanFq -t $thread -o $outDir -c # output novoCluster.txt in $outDir/sRG
```

### circfull with reference guide method to identify full-length circRNA from consensus sequence
```
circfull  cRG -t $thread -g $genome -a $gtfFile -f $outDir -o  $outDir # output circFL_Normal.txt in $outDir/cRG
```

### merge results and filter out low-quality circRNA
```
rmsk=rmsk.bed.gz # optional
circfull  mRG  -m $rmsk -t $thread -g $genome -f $cleanFq -r $outDir -c  $outDir -s $outDir -o  $outDir # output circFL_Normal_pass.txt in $outDir/mRG
```

## annotate full-length circRNA
```
circfull anno -b ${outDir}/mRG/circFL_Normal_pass.bed  -a $gtfFile -o ${outDir}/mRG/circFL_Normal_pass_anno.txt
```

## evaluate gene expression
### required files
```
umgtfFile=gencode.v19.annotation.gtf # uncompressed
```
### build gene reference
```
gene_ref='ref_gene.fa'
createGeneRef $umgtfFile $genome $gene_ref
```

### evaluate gene expression with geneExp module
```
circfull geneExp -f $cleanFq -r $gene_ref -o $outDir -t $thread # output geneCountDf.txt in $outDir/geneExp
```

## calling variants in full-length circRNA
```
circfull var -f $cleanFq -g $genome -t $thread -r $outDir -o $outDir
```
