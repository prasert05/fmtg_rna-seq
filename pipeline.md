# **RNA-Seq analysis pipeline**
Script created & modified by Prasert Yodsawat and Jiratchaya Nuanpirom

We provided an analysis pipeline from the beginning. All tools are run on a specific environment and we divided the analysis into 8 parts.
* [Pre-installed software](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#pre-installed-software)
* [Creating working directory file](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#creating-working-directory-file)
* [Creating conda environments with a specific version of package](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#creating-conda-environments-with-a-specific-version-of-package)
* [Part 1: Quality check & Quality control](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#part1-quality-check--quality-control)
* [Part 2: In silico read normalization & *De novo* transcriptome assembly ](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#part2-in-silico-read-normalization--de-novo-transcriptome-assembly)
* [Part 3: Finding and removing vector contamination](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#part3-finding-and-removing-vector-contamination)
* [Part 4: Finding coding regions](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#part4-finding-coding-regions)
* [Part 5: Transcripts clustering](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#part5-transcripts-clustering)
* [Part 6: Foreign sequences removal](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#part6-foreign-sequences-removal)
* [Part 7: Assembly completeness evaluation](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#part7-assembly-completeness-evaluation)
* [Part 8: Assembly statistics](https://github.com/prasert05/fmtg_rnaseq/blob/main/pipeline.md#part8-assembly-statistics)

** **

## Pre-installed software
Bioconda, bioinformatics tools repository belongs to conda project, was preferred to use as it easy to manage the downloaded recipes and working environment. To install bioconda, click on the badge [![Bioconda](https://img.shields.io/static/v1?label=install&message=bioconda&color=brightgreen)](https://bioconda.github.io/)

## Creating working directory file
```bash
# For quality check
mkdir 01_quality_check 01_quality_check/before 01_quality_check/after

# For read quality and adapter sequence removal
mkdir 02_trim

# For in silico read normalization
mkdir 03_normalization

# For de novo assembly
mkdir 04_assembly_trinity

# For vector decontamination
mkdir 05_vecscreen

# For coding region prediction
mkdir 06_coding

# For transcripts clustering
mkdir 07_cluster

# For foreign sequences removal
mkdir 08_final

# For evaluating transcript completeness
mkdir 09_evaluate

# For assembly statistics
mkdir 10_stats
```

## Creating conda environments with a specific version of package
```bash
# Create environment named 'qc' containing FastQC (v0.11.8), Cutadapt (v2.7) and MultiQC (v1.8) packages
conda create --name qc fastqc=0.11.8 cutadapt=2.7 multiqc=1.8 -y

# Create environment named 'assembly' containing Trinity (v2.8.5) package
conda create --name assembly trinity=2.8.5 -y

# Create environment named 'vecscreen' containing blast (v2.10.1) and SeqKit (v0.14.0) packages
conda create --name vecscreen blast=2.10.1 seqkit=0.14.0 -y

# Create environment named 'coding' containing transdecoder (v5.5.0), r-ggplot2 (v3.3.2) and bioconductor-seqlogo (v1.56.0) packages
conda create --name coding transdecoder=5.5.0 r-ggplot2=3.3.2 bioconductor-seqlogo=1.56.0 -y

# Create environment named 'cluster' containing CD-HIT (v4.8.1) package
conda create --name cluster cd-hit=4.8.1 -y

# Create environment named 'final' containing SeqKit (v0.14.0) package
conda create --name final seqkit=0.14.0 -y
	
# Create environment named 'evaluate' containing BUSCO (v4.1.4) packages
conda create --name evaluate busco=4.1.4 -y

# Create environment named 'stats' containing SeqKit (v0.14.0) packages
conda create --name stats seqkit=0.14.0 -y
```

** **
## **[PART1] Quality check & Quality control**

RNA-Seq reads were checked before and after quality control (QC). 

### Activate conda environment: "qc"
```bash
conda activate qc
```

### Quality check before QC

Quality checking was performed by FastQC, and the QC results were combined with multiQC.

```bash
fastqc \
-t 32 \
-o 01_quality_check/before/ \
raw_data/*

cd 01_quality_check/before
multiqc --filename before_qc_multiqc_report.html .
cd ../../
```

### Quality control (QC)

RNA-Seq reads were trimmed by preserving paired reads contain quality score greater than or equal to 20, and retaining minimum length to 25 bp. Adapter sequences were used based on the library preparation kit e.g. Illumina TruSeq for this study. Please see [Cutadapt user guide](https://cutadapt.readthedocs.io/en/stable/guide.html#user-guide) and [Illumina Adapter Sequences Document](https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html) for more information.

```bash
for (( i = 0; i <= 3; i++ )); do
for (( j = 1; j <= 2; j++ )); do
cutadapt \
--cores 32 \
--quality-cutoff 20,20 --minimum-length 25 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o '02_trim/trimmed.TgS'$i'_rep'$j'_1.fq' \
-p '02_trim/trimmed.TgS'$i'_rep'$j'_2.fq' \
'raw_data/TgS'$i'_rep'$j'_1.fq' \
'raw_data/TgS'$i'_rep'$j'_2.fq' \
> '02_trim/TgS'$i'_rep'$j'.log'
done
done
```

### Quality check after QC
```bash
fastqc \
-t 32 \
-o 01_quality_check/after/ \
02_trim/trimmed.TgS*.fq

cd 01_quality_check/after
multiqc --filename after_qc_multiqc_report.html .
cd ../../
```

### Deactivate conda environments
```bash
conda deactivate
```

** **
## **[PART2] In silico read normalization & *De novo* transcriptome assembly**

Before assembly, the RNA-Seq reads can be reduced due to the exceeding size of total reads. This can save the computational time when assembling a large dataset. The *in silico* read normalization is set as the default option when assemble with Trinity. In our experience, with the limitation of the computational resource, the main Trinity assembly process is usually crashed when Jellyfish is performed, due to the RAM insufficiency. Therefore, the process *in silico* read normalization should be run separately not only customize the usage of RAM and CPUs but also maximize coverage of reads before going to the main assembly process. For more information, please see [Trinity's In silico Read Normalization](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Insilico-Normalization).

### Activate conda environment: "assembly"
```bash
conda activate assembly
```

### *in silico* read normalization before assembly

```bash
insilico_read_normalization.pl \
--CPU 30 --JM 62G \
--seqType fq --max_cov 30 --pairs_together \
--left \
02_trim/trimmed.TgS0_rep1_1.fq,\
02_trim/trimmed.TgS0_rep2_1.fq,\
02_trim/trimmed.TgS1_rep1_1.fq,\
02_trim/trimmed.TgS1_rep2_1.fq,\
02_trim/trimmed.TgS2_rep1_1.fq,\
02_trim/trimmed.TgS2_rep2_1.fq,\
02_trim/trimmed.TgS3_rep1_1.fq,\
02_trim/trimmed.TgS3_rep2_1.fq \
--right \
02_trim/trimmed.TgS0_rep1_2.fq,\
02_trim/trimmed.TgS0_rep2_2.fq,\
02_trim/trimmed.TgS1_rep1_2.fq,\
02_trim/trimmed.TgS1_rep2_2.fq,\
02_trim/trimmed.TgS2_rep1_2.fq,\
02_trim/trimmed.TgS2_rep2_2.fq,\
02_trim/trimmed.TgS3_rep1_2.fq,\
02_trim/trimmed.TgS3_rep2_2.fq \
--output 03_normalization/
```

### Transcriptome assembly
```bash
Trinity \
--CPU 30 --max_memory 60G \
--seqType fq --no_normalize_reads \
--left 03_normalization/left.norm.fq \
--right 03_normalization/right.norm.fq \
--output 04_assembly_trinity/
```

### Deactivate conda environments
```bash
conda deactivate
```

** **
## **[PART3] Finding and removing vector contamination**

These processes were performed to make sure that the data was not contaminated by any artificial sequences including the remaining adapter sequences not passed in the QC process and primer sequences used to amplify the cDNA library before sequencing. These contaminants were removed using standalone VecScreen tool and NCBI UniVec database.

### Activate conda environment: "vecscreen"
```bash
conda activate vecscreen
```

### Download the VecScreen standalone, the script to filter the VecScreen results and adapter database
```bash
cd 05_vecscreen/
wget https://ftp.ncbi.nlm.nih.gov/blast/demo/vecscreen
chmod +x vecscreen
wget ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/VSlistTo1HitPerLine.awk
chmod +x VSlistTo1HitPerLine.awk
wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
cd ..
```

### Make BLAST database for VecScreen standalone
```bash
makeblastdb \
-dbtype nucl \
-parse_seqids \
-blastdb_version 4 \
-in 05_vecscreen/UniVec \
-out 05_vecscreen/UniVec
```

### Screening vector contamination
```bash
05_vecscreen/vecscreen \
-f 3 \
-d 05_vecscreen/UniVec \
-i 04_assembly_trinity/Trinity.fasta \
-o 05_vecscreen/VecScreen_output.txt
```

### Disriminating hits to "Weak", "Suspect Origin", "no hits" and "errors" from result and filtering out
```bash
05_vecscreen/VSlistTo1HitPerLine.awk \
weak=0 suspect=0 no_hits=0 errors=0 \
05_vecscreen/VecScreen_output.txt \
> 05_vecscreen/filtered_VecScreen_output.txt
```

### Extract list of sequences with "Strong" and "Moderate" match adapters to files
```bash
awk \
'{print $2}' \
05_vecscreen/filtered_VecScreen_output.txt | sort | uniq \
> 05_vecscreen/seqlist_filtered_VecScreen_output.txt
```

### Filter out sequences with adapters
```bash
seqkit \
grep -v \
-f 05_vecscreen/seqlist_filtered_VecScreen_output.txt \
04_assembly_trinity/Trinity.fasta \
-o 05_vecscreen/vecscreen_Trinity.fasta
```

### Deactivate conda environments
```bash
conda deactivate
```

* **
## **[PART4] Finding coding regions**
### Activate conda environment: "coding"
```bash
conda activate coding
```

### Extract the long open reading frames
```bash
TransDecoder.LongOrfs \
-t 05_vecscreen/vecscreen_Trinity.fasta \
--output_dir 06_coding/
```

### Predict the likely coding regions
```bash
TransDecoder.Predict \
-t 05_vecscreen/vecscreen_Trinity.fasta \
--output_dir 06_coding/ \
--single_best_only
```

### Move result files
```bash
mv pipeliner* 06_coding/
mv *.transdecoder.* 06_coding/
```

### Deactivate conda environments
```bash
conda deactivate
```

** **
## **[PART5] Transcripts clustering**
### Activate conda environment: "cluster"
```bash
conda activate cluster
```

### Cluster redundant sequences
```bash
cd-hit-est \
-T 30 -M 60000 \
-c 0.95 -G 0 -aS 1.00 \
-d 0 \
-i 06_coding/vecscreen_Trinity.fasta.transdecoder.cds \
-o 07_cluster/cdhit_coding_vecscreen_Trinity.fasta
```

### Deactivate conda environments
```bash
conda deactivate
```

** **
## **[PART6] Foreign sequences removal**

To make sure that the transcript assembly was not composed of many foreign sequence reads contaminated in our RNA-Seq reads, the assembled transcripts were searched for the sequences belong to the species that might contaminated from the process during the sample collection and library preparation. The foreign sequences can be the sequences of bacteria or fungi on the surface, from human, mouse or rat, and plastids like chloroplast or mitochondria. These contaminants can be removed by searching the assembled transcripts against nucleotide databases and filter out. 

### Activate conda environment: "final"
```bash
conda activate final
```

### Remove FASTA description after first whitespace
```bash
seqkit replace \
--pattern "\s.+" \
--out-file 08_final/final_cdhit_coding_vecscreen_Trinity.fasta \
07_cluster/cdhit_coding_vecscreen_Trinity.fasta
```

### Remove sequence that flag as "Exclude", "mask/trim" and "duplicated" during TSA submission (remove_list.txt)
```bash
seqkit \
grep -v \
-f 08_final/remove_list.txt \
08_final/final_cdhit_coding_vecscreen_Trinity.fasta \
-o 08_final/fil_final_cdhit_coding_vecscreen_Trinity.fasta
```

### Deactivate conda environments
```bash
conda deactivate
```

** **
## **[PART7] Assembly completeness evaluation**
### Activate conda environment: "evaluate"
```bash
conda activate evaluate
```

### Transcript evaluation by BUSCO using "arthropoda_odb10" as lineage dataset
```bash
busco \
--cpu 30 \
--mode tran \
--in 08_final/fil_final_cdhit_coding_vecscreen_Trinity.fasta \
--lineage_dataset arthropoda_odb10 \
--out busco_output
```

### Move output folder
```bash
mv busco_* 09_evaluate/
```

### Deactivate conda environments
```bash
conda deactivate
```

** **
## **[PART8] Assembly statistics**
### Activate conda environment: "evaluate"
```bash
conda activate evaluate
```
### Generate FASTA statistics of raw and final transcripts assembly
```bash
seqkit stats --all --tabular \
04_assembly_trinity/Trinity.fasta > 10_stats/stats_all.txt

seqkit stats --all --tabular \
08_final/fil_final_cdhit_coding_vecscreen_Trinity.fasta | \
tail -n 1 \
>> 10_stats/stats_all.txt
```

### Deactivate conda environments
```bash
conda deactivate
```

** **
