*Load packages and modules.* 
```conda activate virtual_environment
export PATH="/rhome/username/.local/bin:$PATH"
pip install quast
pip install pRESTO
module load trim_galore/0.4.2
module load bwa/0.7.17
module load samtools/1.9
module load bedtools/2.29.2
module load BBTools/38.60
moduel load fastqc/
```
*Create an index of your reference genome*
```
bwa index -p SARS_CoV_2_Wuhan_NC_045512_2 -a is SARS_CoV_2_Wuhan_NC_045512_2.fasta
```
*Run fastqc on all raw files*
```
for file in *.gz; do fastqc $file; done
```
*Trim adapters from reads. Quality filter.*
```
module load trim_galore/0.4.2
for file in *gz; do  trim_galore -q 32 --illumina --length 120 $file; done
for file in *.txt.gz_trimmed.fq.gz; do cp $file "${file%.txt.gz_trimmed.fq.gz}"_trimmed.fq.gz; done
for file in *.txt.gz_trimmed.fq.gz; do rm $file; done

for file in *trimmed.fq.gz; do gunzip $file; done
```
*Extract number of reads from trimming reports*
```
for file in *.txt; do sed -n '4p;23p;25p' $file > "${file%.test}"_stats.txt; done 
cat *_stats.txt > all_trimming_stats.txt
awk '{ printf("%-10s ",$0); if(NR%3==0) printf("\n");}' all_trimming_stats.txt > all_trimming_stats_tabbed.xls
```
*Remove primers.*
Primer files must be separated by "forward" (Vprimer") and "reverse" (Cprimer") and must be in fasta format within the same directory. Use max error rate >0.40. 
```
for file in *R1_trimmed.fq; do
MaskPrimers.py align -s $file -p /bigdata/messaoudilab/apinski/Scripts/viral_genome_sequencing/forwardPrimers.fa \
    --maxlen 50 --maxerror 0.55 --mode cut --pf VPRIMER \
    --outname "${file%.fq}" --log "${file%.fq}"_Vprimer.log
ParseLog.py -l "${file%.fq}"_Vprimer.log -f ID PRIMER ERROR;
done

for file in *R2_trimmed.fq; do
MaskPrimers.py align -s $file -p /bigdata/messaoudilab/apinski/Scripts/viral_genome_sequencing/reversePrimers.fa \
    --maxlen 50 --maxerror 0.55 --mode cut --pf CPRIMER --revpr --skiprc \
    --outname "${file%.fq}" --log "${file%.fq}"_Cprimer.log
ParseLog.py -l "${file%.fq}"_Cprimer.log -f ID PRIMER ERROR; 
done
```
*Extract primer removal data*
```
grep -w 'SEQ_FILE\|SEQUENCES\|PASS\|FAIL' removePrimers.stdout | awk 'ORS=NR%6?FS:RS' | tr ' ' \\t | awk  '{print $2,$4,$6,$8}' | tr ' ' \\t > remove_primer_info.txt
```
*Subsample here if sequences > 1 M per* 
```
module load seqtk/281537b 
seqtk sample -s100 read1.fq 1000000 > sub1.fq
```
*Concatenate Read1 and Read2*
```
for file1 in *R1_trimmed.fq.gz; do file2=${file1/R1_trimmed.fq.gz/R2_trimmed.fq.z} cat $file1 $file2 > cat_$file1; done
```
*Align reads to reference genome*
```
for file in *fq; do
bwa mem -M -t 48 SARS_CoV_2_Wuhan_NC_045512_2 $file > "${file%.fq.gz}".sam;
done
for file in *sam; do samtools view -S -b $file > "${file%.sam}".bam; done
for file in *.bam; do samtools sort $file -o "${file%.bam}".sorted.bam; done
for file in *.sorted.bam; do  bedtools bamtofastq -i "$file" -fq "${file%.bam}".fq; done
for file in *sorted.bam; do samtools index $file "${file%.bam}".sorted.bam.bai; done
```
*Assess coverage: create the coverage.R script then run the following*
```
for file in *trimmed.sorted.bam; do
samtools depth $file > "${file%.sorted.bam}".coverage
bamCoverage -b $file -o "${file%.bam}".bw;
done
for file in *.coverage; do Rscript coverage.R $file "${file%.coverage}"; done 
```
######coverage.R
```
args <- commandArgs(TRUE)
library(reshape)
library(lattice, pos=10)
library(ggplot2)
args <- commandArgs(TRUE)
file1 <- args[1]
outprefix <- args[2]
f1<- read.table(file1, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
f1<- rename(f1, c(V1="chr", V2="position", V3="depth"))
#Find the stats on coverage
summary(f1$depth)
f1summarydf <-data.frame(unclass(summary(f1$depth)), check.names=FALSE, stringsAsFactors=FALSE)
write.table(f1summarydf, paste(outprefix,"_summary.txt", sep=""), sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
``` 
#Create histogram of coverage. Change Y axis as needed under scale_y_continous
```
pdf(paste(outprefix,"_coverage_histogram.pdf", sep=""))
plot_f1 <- ggplot(f1, aes(x=f1$depth))+geom_histogram(binwidth=1, position="identity")+xlab("depth")+ylab("count")+ggtitle(paste(outprefix, "SARS-CoV-2 Depth Distribution", sep=" "))+theme(plot.title=element_text(hjust=0.5))+theme_bw()+scale_x_continuous(breaks=c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000), limits=c(0, 8500))+scale_y_continuous(limits=c(0, 250))
plot_f1
dev.off()
```
#Density
```
pdf(paste(outprefix, "_coverage.pdf", sep=""))
xyplot(depth ~ position, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'), y=list(relation='same')), data=f1, main="SARS-CoV-2 Coverage by Position")
dev.off()
```
*Call variants*
```
for file in *sorted.bam; do samtools mpileup -uf /bigdata/messaoudilab/apinski/Scripts/viral_genome_sequencing/SARS_CoV_2_Wuhan_NC_045512_2.fasta  $file | bcftools call -c  | vcfutils.pl vcf2fq > "${file%.sorted.bam}"_consensus.fq; done
```
*Assess assembly and extract files*
```
for file in *.fasta; do quast.py $file -r /bigdata/messaoudilab/apinski/Scripts/viral_genome_sequencing/SARS_CoV_2_Wuhan_NC_045512_2.fasta  -o "${file%.fasta}"_quast_results;  done
for subdir in *; do cp $subdir/report.txt $subdir.txt; done
```
