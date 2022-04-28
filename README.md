A Quick RNA-seq Pipeline Tutorial
================

**Keren Li**, \<<li.keren.cn@gmail.com>\>


## What is RNA-seq?

RNA-Seq (named as an abbreviation of RNA sequencing) is the leading technique which uses next-generation sequencing to reveal the presence and quantity of RNA in a biological sample at a given moment, analyzing the continuously changing cellular transcriptome. ([Wikipedia](https://en.wikipedia.org/wiki/RNA-Seq))

Accurate alignment of high-throughput RNA-seq data is a challenging and yet unsolved problem because of the non-contiguous transcript structure, relatively short read lengths and constantly increasing throughput of the sequencing technologies.

This tutorial will guide you through the procedures from downloading RNA-seq reads to alignment on QUEST, as well as differential gene expression analysis. You can also do similar operations on your laptop/desktop.


## Get to know your QUEST allocation


Northwestern University's high performance computing (HPC) system is called **QUEST**. It is used for computationally intensive tasks and extremely large datasets. Students in STAT465 are invited to use QUEST alocation **`e31675`**. 

#### Logging in to QUEST

To login to Quest from a Linux or Mac computer, open a terminal and at the command line enter:

```bash
$ ssh -X <netid>@quest.northwestern.edu
```

Windows users will first need to download an ssh client, such as [_PuTTY_](https://www.chiark.greenend.org.uk/~sgtatham/putty/), which will allow you to interact with the remote Unix command line. Use the following information to set up your connection: <br />
**Hostname** : quest.northwestern.edu  <br />
**Username** : your Northwestern NetID  <br />
**Password** : your Northwestern NetID password  

If you'll be working with GUI (graphical window) applications, A FastX client is recommended.

#### Getting started with QUEST

To get started with QUEST, let's visit the allocation folder and create a new folder `RNAseqExample` along with three subfolder within.

```bash
$ cd /projects/e31675
$ mkdir RNAseqExample
$ cd RNAseqExample/
$ mkdir genome gtf scripts reads bams
```
The samples used in this tutorial are `SRR391535`, `SRR391536`, `SRR391537`, `SRR391538`, `SRR391539`, and `SRR391541`from [_Comparing the response to chronic ozone of five agriculturally important legume species_](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP009826) at NCBI, and much credit goes to [this RNA-seq tutorial](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/rna-seq-tutorial-with-reference-genome/).

Let's first download the reference genome and annotation file for Glycine max (soybean) from [EnsemblePlants](http://plants.ensembl.org/Glycine_max/Info/Index). 

```bash
$ cd gtf/
$ wget http://ftp.ensemblgenomes.org/pub/plants/release-53/fasta/glycine_max/dna/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz
$ gunzip Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz
$ wget http://ftp.ensemblgenomes.org/pub/plants/release-53/gff3/glycine_max/Glycine_max.Glycine_max_v2.1.53.gff3.gz
$ gunzip Glycine_max.Glycine_max_v2.1.53.gff3.gz
```

Let's inspect content of `Glycine_max.Glycine_max_v2.1.53.gff3`. `less` is a command line utility similar to `more`, but has more advanced features and allows you to navigate both forward and backward through the file.

```bash
less Glycine_max.Glycine_max_v2.1.53.gff3
```

The GFF3 file downloaded has a messed up `Parent` field. This is unusual situation. We need remove strings: `gene:`, `transcript:`, and `CDS:` and convert the GFF3 file to GTF file to make STAR work properly later.

```bash
$ sed -i 's/gene://g;s/CDS://g;s/transcript://g;' Glycine_max.Glycine_max_v2.1.53.gff3
$ module load cufflinks
$ gffread Glycine_max.Glycine_max_v2.1.53.gff3 -T -o Glycine_max.Glycine_max_v2.1.53.gtf
```
To call any package pre-installed on QUEST, load the corresponding module. In the above bash script, we load module `cufflinks` to call `gffread` app.

#### Downloading SRA files

`fastq-dump` is a tool for downloading sequencing reads from NCBI’s Sequence Read Archive (**SRA**). These sequence reads will be downloaded as FASTQ files. Load module `sratoolkit` to call `fastq-dump` to download `SRR391535` .SRA file and converting it to FASTQ.

```bash
$ cd reads/
$ module load sratoolkit
$ fastq-dump  SRR391535
```

`SRR391535` is a single-end RNA-seq data. To download paired-end RNA-seq data, an argument `--split-files` is required for `fastq-dump`, which splits the FASTQ reads into two files: one file for mate 1s (`...1`), and another for mate 2s (`...2`). 


#### Submit your scripts to QUEST

Although we can call apps directly in the command line, it is recommended to create scripts to be submitted to QUEST.


Now let’s create your first script `download.sh` in the path `/projects/e31675/RNAseqExample/scripts/`:

```bash
$ cd ../scripts/
$ nano download.sh
```

Within `nano` we add the shebang line and the SLURM directives, as well as `fastq-dump` command. **Remember** to replace `$YOUR_EMAIL` by your northwestern email.

```bash
#!/bin/sh
#SBATCH -A e31675
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mail-user=$YOUR_EMAIL
#SBATCH --job-name="download"
#SBATCH --mem=100000

module purge all

module load sratoolkit

cd /projects/e31675/RNAseqExample/reads/
fastq-dump SRR391535
fastq-dump SRR391536
fastq-dump SRR391537 
fastq-dump SRR391538 
fastq-dump SRR391539 
fastq-dump SRR391541 
```

To submit your first script to QUEST, run

```bash
$ sbatch download.sh
Submitted batch job #######
```
`#######` from the screen output is the QUEST job id assigned.

After a cup of coffee (or two cups!), come back to check the path `/projects/e31675/RNAseqExample/reads/`, or check your job status directly:

```bash
$ checkjob #######
```

For details of QUEST please see [**Quest Overview**](https://www.it.northwestern.edu/research/user-services/quest/overview.html).

## STAR: alignment tool for RNA-seq

**STAR** (Spliced Transcripts Alignment to a Reference) is a popular aligner designed to address spliced alignment problems for RNA-seq data mapping.

#### STAR alignment strategy

Credit of this section goes to [this tutorial](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html). 

The central idea of the STAR seed finding phase is the sequential search for a Maximal Mappable Prefix (MMP), the longest matching sequences.

![Fig1](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step1.png) 

The different parts of the read that are mapped separately are called ‘seeds’. So the first MMP that is mapped to the genome is called seed1.

STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP, which will be seed2.

![Fig2](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step2.png) 

This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping.

If STAR does not find an exact matching sequence for each part of the read due to mismatches or indels, the previous MMPs will be extended.

![Fig3](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step3.png) 


If extension does not give a good alignment, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.

![Fig4](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step4.png) 


The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of ‘anchor’ seeds, or seeds that are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.).

![Fig5](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step5.png) 

To align with STAR, there are two steps: generate a genome index and align reads to the genome.

#### Generating a genome index

The following script is used to generate a genome index for an RNA-seq alignment purpose.

```bash
STAR --runMode genomeGenerate \  
--runThreadN 10 \
--genomeDir $GenomePath \
--genomeFastaFiles $FastaFile \
--sjdbGTFfile $GTFFile \
--sjdbOverhang 100 
```

The options used to geneate a genome index are as follows: <br />
`--runMode`: choose `genomeGenerate` mode to generate a genome index. Its default value is `alignReads`. <br />
`--runThreadN`: 10 threads are used (set 2 or 4 threads for your own laptop/desktop).  <br />
`--genomeDir`: the genome index is stored in `GenomePath`. <br />
`--genomeFastaFiles`: `FastaFile` is the file name of the FASTA including the full path.  <br />
`--sjdbGTFfile`: `GTFFile` is the file name of the annotated transcript file including the full path. It is optional, though using annotations is highly recommended for splicing alignment by STAR. <br />
`--sjdbOverhang`: specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, it should be the maximum read length -1. 

If a GFF3 formatted annotations instead of GTF format is used, `--sjdbGTFtagExonParentTranscript Parent` is required.

Create a STAR script to geneate a genome index for yeast:

```bash
$ nano STAR_generate.sh
```

```bash
#!/bin/sh
#SBATCH -A e31675
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mail-user=$YOUR_EMAIL
#SBATCH --job-name="genome_build"
#SBATCH --mem=100000

module purge all
module load STAR

cd /projects/e31675/RNAseqExample/

STAR --runMode genomeGenerate \
--runThreadN 10 \
--genomeDir ./genome/ \
--genomeFastaFiles ./gtf/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa \
--outFileNamePrefix Glycine_max \
--sjdbGTFfile ./gtf/Glycine_max.Glycine_max_v2.1.53.gtf \
--sjdbOverhang 100
```

**Important**: be ware when you copy and paste those scripts. STAR cannot handle extra empty spaces very well.

Submit the above script to QUEST 

```bash
$ sbatch STAR_generate.sh
```

#### Aligning reads to the genome 

Once you build the genome index, RNA-seq reads can be aligned using the following script.

```bash
STAR --runThreadN 10 \
--genomeDir $GenomePath \
--sjdbGTFfile $GTFFile \
--readFilesIn $ReadFile \
--readFilesCommand zcat \
--readFilesPrefix $FilePrefix \
--limitBAMsortRAM 100000000000 \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts 
```

The new options in read alignment are  as follows: <br />
`--readFilesIn`: `ReadFile` contains names of the RNA-seq files, e.g. FATSA abd FASTQ files. For paired-end reads, both files have to be supplied. <br />
`--readFilesCommand`: `zcat` is the uncompression command if the input read files are in gzipped format (*.gz). <br />
`--readFilesPrefix`: `FilePrefix` will be added in front of the input file names. <br /> 
`--limitBAMsortRAM`: maximum available RAM (bytes) for sorting BAM. <br /> 
`--outFilterMultimapNmax`: max number of multiple alignments allowed for a read.  The read is considered unmapped if exceeded. We only take uniquely mapped reads by setting value 1. <br />
`--outSAMtype`: option `BAM SortedByCoordinate` outputs sorted by coordinate bam files. <br /> 
`--quantMode`: with option `GeneCounts`, STAR will also count number reads per gene while mapping. 


If you want to turn off the splicing alignments for purposes other than RNA-seq data, please re-generate the genome index without annotations argument `--sjdbGTFfile` in the previous step and use the following options to align. <br />
`--alignIntronMax`: maximum intron size. Use 1 to turn off splicing. <br />
`--alignMatesGapMax`: maximum genomic distance between mates.

Now we are ready for alignment through STAR.

```bash
$ nano STAR_mapping.sh
```

```bash
#!/bin/sh
#SBATCH -A e31675
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mail-user=$YOUR_EMAIL
#SBATCH --job-name="alignment"
#SBATCH --mem=100000

module purge all
module load STAR

cd /projects/e31675/RNAseqExample/reads/

for file in *.fastq
do
file1=${file/.fastq/}
STAR --runThreadN 10 \
--genomeDir ../genome/ \
--sjdbGTFfile ../gtf/Glycine_max.Glycine_max_v2.1.53.gtf \
--readFilesIn $file \
--outFileNamePrefix ../bam/$file1 \
--limitBAMsortRAM 100000000000 \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--alignIntronMin 10 \
--alignIntronMax 5000
done
```
To submit the alignment script on QUEST, run

```bash
$ sbatch STAR_mapping.sh
```

RNA-seq reads aligned on genes are now stored in the `*ReadsPerGene.out.tab` files.



## Exploring BAM files with `samtools`

A BAM file (.bam) is the binary version of a SAM file. A SAM file (.sam) is a tab-delimited text file that contains sequence alignment data. 

`samtools` is a suite of programs for interacting with high-throughput sequencing data. To view and convert SAM/BAM/CRAM files, use `samtools view`. With no options or regions specified, prints all alignments in the specified input alignment file (in SAM, BAM, or CRAM format) to standard output in SAM format (with no header).

```bash
$ module load samtools
$ cd /projects/e31675/RNAseqExample/bam
$ samtools view SRR391535Aligned.sortedByCoord.out.bam | less
```
Here `|` is the piping symbol, which is the process of redirecting the output of one command to the input of another command. 

First we create the index file (BAI) for the BAM file with the following command:
```bash
$ samtools index SRR391535Aligned.sortedByCoord.out.bam
```
It generates an index file `SRR391535Aligned.sortedByCoord.out.bam.bai`. Then we extract the data for reference (chromosome) `2`:
```bash
$ samtools view SRR391535Aligned.sortedByCoord.out.bam 2 | less
```
Now save aligned result in (chromosome) `2` to a BAMfile using `>` symbol with specified columns. `>` takes the output of a command and redirects it into a file.

```bash
$ samtools view -h SRR391535Aligned.sortedByCoord.out.bam 2 | awk '{print $3, "\t", $4, "\t", $8}' > SRR391535.chr2.sam
```
`-h` means include the header in the output.

#### SAM format

The SAM format consists of a header and an alignment section:

|Col|Field  |	Type	|Brief description                      |
|---| :---  | :---      | :---                                  |
|1	| QNAME	| String	| Query template NAME                   |
|2	| FLAG	| Int	    | bitwise FLAG                          |
|3	| RNAME	| String	| References sequence NAME              |
|4	| POS	| Int	    | 1- based leftmost mapping POSition    |
|5	| MAPQ	| Int	    | MAPping Quality                       |
|6	| CIGAR	| String	| CIGAR string                          |
|7	| RNEXT	| String	| Ref. name of the mate/next read       |
|8	| PNEXT	| Int	    | Position of the mate/next read        |
|9	| TLEN	| Int	    | observed Template LENgth              |
|10	| SEQ	| String	| segment SEQuence                      |
|11	| QUAL	| String	| ASCII of Phred-scaled base QUALity+33 |

#### CIGAR Format

|CIGAR  |Code   |	BAM Integer	Description                                 |Consumes query |Consumes reference |
|---    | :---  | :---                                                      | :---          |:---               |
|M	    |0      |	alignment match (can be a sequence match or mismatch)   |	yes         |	yes             |
|I	    |1      |	insertion to the reference                              |	yes         |	no              |
|D	    |2      |	deletion from the reference                             |	no          |	yes             |
|N	    |3      |	skipped region from the reference                       |	no          |	yes             |
|S	    |4      |	soft clipping (clipped sequences present in SEQ)        |	yes         |	no              |
|H	    |5      |	hard clipping (clipped sequences NOT present in SEQ)    |	no          |	no              |
|P	    |6      |	padding (silent deletion from padded reference)         |	no          |	no              |
|=	    |7      |	sequence match                                          |	yes         |	yes             |
|X	    |8      |	sequence mismatch                                       |	yes         |	yes             |


## Analysis of Counts with `DESeq2`:

`**DESeq2**` is a differential gene expression analysis R package based on the negative binomial distribution.



#### Installing `DESeq2`

It could be easier for you to work on your personal computer rather than the server. So you can download the `*ReadsPerGene.out.tab` files from STAR alignment results onto your computer. 

We first install `DESeq2` package on [Rstudio server](https://rstudio.questanalytics.northwestern.edu/):

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```
Load necessary packages and set working directory:

```R
library("DESeq2")
library("data.table")
setwd("/projects/e31675/RNAseqExample/bam/")
```

#### Loading gene counts

Read gene counts from the `*ReadsPerGene.out.tab` files:

```R
countFileList <- list.files(pattern = "ReadsPerGene.out.tab$")
replicateNames <- gsub("ReadsPerGene.out.tab", "", countFileList)

countGene <- NULL
for(countFile in countFileList){
  countGeneRep <- fread(countFile, header = FALSE, skip =4)
  countGene <- cbind(countGene, countGeneRep$V2)
}
colnames(countGene) <- replicateNames
dim(countGene)
# [1] 57147     6
```
First four lines of the `*ReadsPerGene.out.tab` files are some summaries, so we choose to skip the first four lines.

#### Constructing DESEQDataSet object

Design specifies how the counts from each gene depend on our variables in the metadata.
For this dataset the factor we care about is our treatment status (`condition`).

```R
sampleNames <- c("Leaf tissue 1", "Leaf tissue 2", "Leaf tissue 3",
                 "Leaf tissue 4", "Leaf tissue 5", "Leaf tissue 6")
sampleCondition <- c("ambient", "ambient", "elevated",
                     "elevated", "elevated", "ambient")
metaData <- data.frame(sampleName = sampleNames, 
                       fileName = replicateNames)
metaData$condition <- factor(sampleCondition)
metaData
```
```R
##      sampleName  fileName condition
## 1 Leaf tissue 1 SRR391535   ambient
## 2 Leaf tissue 2 SRR391536   ambient
## 3 Leaf tissue 3 SRR391537  elevated
## 4 Leaf tissue 4 SRR391538  elevated
## 5 Leaf tissue 5 SRR391539  elevated
## 6 Leaf tissue 6 SRR391541   ambient
```
```R
dds <- DESeqDataSetFromMatrix(countData = countGene,
                              colData = metaData,
                              design = ~ condition)
dds
```
```R
## class: DESeqDataSet 
## dim: 57147 6 
## metadata(1): version
## assays(1): counts
## rownames(57147): GLYMA_01G000100 GLYMA_01G000200 ... GLYMA_U009100 GLYMA_U022500
## rowData names(0):
## colnames(6): SRR391535 SRR391536 ... SRR391539 SRR391541
## colData names(3): sampleName fileName condition
```
Now we are ready to run `DESeq()` function:

```R
dds <- DESeq(dds)
```
```R
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
```
Assign differential expression results to `res`. Column `padj` is the adjust p-value. Check 
```R
res <- results(dds)
res[1:5]
```

```R
## log2 fold change (MLE): condition elevated vs ambient 
## Wald test p-value: condition elevated vs ambient 
## DataFrame with 5 rows and 6 columns
##                   baseMean log2FoldChange     lfcSE       stat    pvalue      padj
##                  <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
## GLYMA_01G000100   7.125903      1.4344115  0.792251  1.8105521 0.0702102  0.157418
## GLYMA_01G000200   8.316737     -1.3244503  0.806630 -1.6419550 0.1005993  0.207693
## GLYMA_01G000300   0.293603     -0.0522698  4.080473 -0.0128097 0.9897796        NA
## GLYMA_01G000400 232.363948     -0.0377731  0.158426 -0.2384267 0.8115501  0.889824
## GLYMA_01G000500   2.205610     -0.7472891  1.332581 -0.5607834 0.5749452  0.718275
```
Sort output by `padj` and save it to file "Gmax_DESeq2.txt".

```R
resNames <- rownames(res)
res <- as.data.table(res)
res[, gene:=resNames]
setcolorder(res, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
setorder(res, padj, na.last=TRUE)
res[padj<=0.05]
```

```R
##                   gene     baseMean log2FoldChange      lfcSE      stat        pvalue          padj
##     1: GLYMA_13G097900  1210.331930      4.3090027 0.15133852 28.472610 2.558510e-178 9.717478e-174
##     2: GLYMA_09G049100  4192.189573      5.0867612 0.18397169 27.649695 2.814424e-168 5.344732e-164
##     3: GLYMA_06G141600  1936.964864      3.4357875 0.12809214 26.822782 1.752779e-158 2.219077e-154
##     4: GLYMA_08G277000   853.382940      4.2198096 0.16140063 26.144938 1.125204e-150 1.068409e-146
##     5: GLYMA_01G127200 12227.418568      3.2534156 0.12598843 25.823129 4.876659e-147 3.704408e-143
##    ---                                                                                             
## 12741: GLYMA_01G053800   135.335931      0.6802200 0.28437770  2.391960  1.675868e-02  4.995774e-02
## 12742: GLYMA_02G061700   247.459076     -0.4839882 0.20234625 -2.391881  1.676226e-02  4.996449e-02
## 12743: GLYMA_17G149600   859.044978      0.2952993 0.12347136  2.391642  1.677319e-02  4.999313e-02
## 12744: GLYMA_10G180600  4769.466904     -0.2339578 0.09782527 -2.391588  1.677565e-02  4.999489e-02
## 12745: GLYMA_04G006600     4.919786     -2.7690815 1.15785012 -2.391572  1.677641e-02  4.999489e-02
```

```R
fwrite(res, file="Gmax_DESeq2.txt")
```

For visualization part, please see `/projects/e31675/RNAseqExample/scripts/DESeq2Ex.R` on QUEST.

**References and links**:

Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics (Oxford, England)*, *29*(1), 15–21. https://doi.org/10.1093/bioinformatics/bts635

Love M. I., Huber W., & Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. doi: 10.1186/s13059-014-0550-8.

Dobin, A. (2022).  [STAR manual 2.7.10a.](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) 

[Introduction to RNA-Seq using high-performance computing.](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html) 

[Quest Overview](https://www.it.northwestern.edu/research/user-services/quest/overview.html)

[RNA-seq Tutorial (with Reference Genome)](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/rna-seq-tutorial-with-reference-genome/)


</div>
