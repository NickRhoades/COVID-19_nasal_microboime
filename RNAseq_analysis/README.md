Basic RNAseq analysis

***Step 1 Trimming using TrimGalore
These trimming parameters were determined by examining a FASTQC report of raw sequencing reads
For all steps "${ID}" is equivelent to individual sample IDs
```
module unload miniconda2/4.4.10
module load fastqc
module load trim_galore
module load python/2.7.5




trim_galore --paired --clip_R1 13 --three_prime_clip_R1 2 --clip_R2 13 --three_prime_clip_R2 2 --gzip -o ./ -q 30 -length 50 -fastqc ${ID}READ1.fastq.gz ${ID}READ2.fastq.gz
```

*** Step 2 Alignment

```
module load hisat2
#module unload python
#module load python/3
module load samtools


mkdir ./${ID}hisat
hisat2 -p 40 -k 1 --min-intronlen 50 --max-intronlen 500000 -x /bigdata/messaoudilab/mzulu/References/Homo_Sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa -1 ./${ID}READ1_val_1.fq.gz  -2 ./${ID}READ2_val_2.fq.gz -S ./${ID}hisat/${ID}sam --summary-file ./${ID}hisat/${ID}_summary.txt



samtools view -bS ./${ID}hisat/${ID}sam > ./${ID}hisat/accepted_hits.bam
```

*** Moving into and environment and Counting reads

Start R in the main directory and Load the required packages again
```
R
library(systemPipeR)
library(GenomicFeatures)
library(BiocParallel)
library(ape)
library(DESeq2)
library(edgeR)
```

- Read in targets file. (This file can be found in this repository and contains information about each sample) 
```
targets <- read.delim("targets.txt", comment.char = "#")
targets
```
- Create the args object again
```
args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
file.exists(outpaths(args))
```

- The following performs read counting with summarizeOverlaps in parallel mode with multiple cores
```
library("GenomicFeatures"); library(BiocParallel)
txdb <- loadDb("./")
eByg <- exonsBy(txdb, by=c("gene"))
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=FALSE, inter.feature=TRUE, singleEnd=TRUE))
```
- Wait until read counting is done, then write countDFeByg into an excel file
```
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(outpaths(args))
write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
```
- Generate RPKM normalized expression values from the countDFeByg file
```
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
```






