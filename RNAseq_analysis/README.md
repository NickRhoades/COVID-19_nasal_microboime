Basic RNAseq analysis

***Step 1 Trimming using TrimGalore***
These trimming parameters were determined by examining a FASTQC report of raw sequencing reads
For all steps "${ID}" is equivelent to individual sample IDs
```
module unload miniconda2/4.4.10
module load fastqc
module load trim_galore
module load python/2.7.5




trim_galore --paired --clip_R1 13 --three_prime_clip_R1 2 --clip_R2 13 --three_prime_clip_R2 2 --gzip -o ./ -q 30 -length 50 -fastqc ${ID}READ1.fastq.gz ${ID}READ2.fastq.gz
```

***Step 2 Alignment***

```
module load hisat2
#module unload python
#module load python/3
module load samtools


mkdir ./${ID}hisat
hisat2 -p 40 -k 1 --min-intronlen 50 --max-intronlen 500000 -x /bigdata/messaoudilab/mzulu/References/Homo_Sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa -1 ./${ID}READ1_val_1.fq.gz  -2 ./${ID}READ2_val_2.fq.gz -S ./${ID}hisat/${ID}sam --summary-file ./${ID}hisat/${ID}_summary.txt



samtools view -bS ./${ID}hisat/${ID}sam > ./${ID}hisat/accepted_hits.bam
```

***Step 3: Moving into and environment and Counting reads***

- Start R in the main directory and Load the required packages again
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
args <- systemArgs(sysma="hisat2.param", mytargets="targets.txt")
file.exists(outpaths(args))
```

- The following performs read counting with summarizeOverlaps in parallel mode with multiple cores
```
txdb <- loadDb("./Homo_Sapiens/Homo_sapiens.sqlite")
eByg <- exonsBy(txdb, by=c("gene"))
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=FALSE, inter.feature=TRUE, singleEnd=FAlSE))
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

***Step 4 DEG analysis***
- DEG Analysis with EdgeR (target file edits can change comparison this should be done fore 2 group and 5 group) 
- Load libraries 
```
library(systemPipeR)
library(edgeR)
```
- Read in raw counts 
```
countDF <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
targets <- read.delim("targets.txt", comment="#")
cmp <- readComp(file="targets.txt", format="matrix", delim="-")
edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")

desc <- read.delim("/bigdata/messaoudilab/mzulu/References/Homo_Sapiens/Human_genes_GRCh38.91.p10.txt", row.names=1)

edgeDF <- cbind(edgeDF, desc[rownames(edgeDF),])
write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)
edgeDF <- read.delim("results/edgeRglm_allcomp_2group.xls", row.names=1, check.names=FALSE)

```

