# *QIIME2 Updated 4/21/20 (last edited 7/4/21)*
This pipeline is based on QIIME2 2019.10
And follows the Moving Pictures tutorial that can be found here. [Moving Picutres Tutorial](https://docs.qiime2.org/2019.10/tutorials/moving-pictures/)
QIIME2 has very good documentation and almost all questions can be answered here. [QIIME2](https://qiime2.org/)

#### **The scripts I have added here are to run QIIME2 analysis in the UCR cluster slurm framework.**
This allow the user to run scripts in the background without having to wait for resources.

## *Importing Data*
To run this batch script you will need to create a directory with all of you sequence files in a fastq.gz. 

Unless you edit the script the directory must be called **seqs**

To use this script file names must be in the Casava format ex. SampleID_L001_R1_001.fastq.gz

For importing other formats of data see. [QIIME2 importing data](https://docs.qiime2.org/2019.10/tutorials/importing/)

While QIIME2 version 2019.10 does not vary significantly for the most up to date version 
it is always suggested that the most up-to-date version of the program is installed/

## **Getting into a QIIME environment:**  
This is just a little trick to make sure that our cluster and QIIME2 are speaking the same language. 
Then entering the QIIME2 virtual environment. (Not required on non-UCR cluter systems.

```
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
conda init bash
source activate qiime2-2019.10
```
## **Generating diversity metrics:**
````
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'   --input-path seqs/   --input-format CasavaOneEightSingleLanePerSampleDirFmt   --output-path demux-paired-end.qza

qiime demux summarize   --i-data demux-paired-end.qza   --o-visualization demux.qzv
````

## **Sequence quality control, feature table construction**
This script combines multiple steps of the qiime pipeline.
In this step we use the program dada2 but there are other options. See QIIME2 documentation for other options and troubleshooting.
It will take a while so I suggest you use the batch script.
Trimming parameters were determined after examination of the demux.qzv file. 
These trimming parameters will vary based on experiment and should be selected based on when quality score beging to significantly fall off.
````
qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end.qza --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --p-trim-left-f 1 --p-trim-left-r 1 --p-trunc-len-f 235 --p-trunc-len-r 180 --p-n-threads n --o-denoising-stats stats-dada2.qza

qiime feature-table summarize --i-table table-dada2.qza --o-visualization table.qzv 
````
## **Assign taxonomy using the SILVA database to remove potential contaminate reads**
While the final taxonomy for this analysis was assigned using the eHOMD database be first used a broad microbial database (SILVA) to identify contaminating sequences.

````
qiime feature-classifier classify-sklearn --i-classifier silva-132-99-nb-classifier.qza --i-reads rep-seqs-dada2.qza --o-classification taxonomy-silva.qza --p-n-jobs 60
````
## **Remove poorly classified reads and contaminating taxa**
We removed any sequence that was not classified to the Phyla level. 
Additionaly we removed any read identified as mitochondria or chlorplast.
Finally we removed two taxa that were found in high abundance within our extraction controls:
1) k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Thermoanaerobacterium;s__saccharolyticum
2) k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Myxococcales;f__0319-6G20;g__;s__
````
qiime taxa filter-table \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy-silva.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast, k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Thermoanaerobacterium;s__saccharolyticum, k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Myxococcales;f__0319-6G20;g__;s__ \
  --o-filtered-table table-no-contamination.qza
````
## **Assign taxonomy using the eHOMD database**

````
qiime feature-classifier classify-sklearn --i-classifier eHOMD_flu_classifier.qza --i-reads rep-seqs-dada2.qza --o-classification taxonomy-eHOMD.qza --p-n-jobs n
````
## **Tree building**
````
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza --p-n-threads n
````
## **Generating diversity metrics:**
You may need to edit the sampling depth number based on your study.
```
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table-no-contamination.qza --p-sampling-depth 8000 --m-metadata-file Nasal_covid_map.txt --output-dir core-metrics-results-8k
```
## **Making taxa bar plots:**
```
qiime taxa barplot --i-table ./ccore-metrics-results-8k/rarified-table.qza --i-taxonomy taxonomy-eHOMD.qza --m-metadata-file Nasal_covid_map.txt --o-visualization taxa-bar-plots-eHOMD-8k.qzv
````
