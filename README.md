# mitoribo
Bioinformatics workflow for ribosome sequencing of yeast

## Installation
The following commands will install Anadama2 (a workflow management library) and this GitHub repository to run the workflow. The GitHub repository contains a `.yml` file to create a conda environment with all the necessary tools. Because TopHat2 requires Python 2 and Anadama2 requires Python 3, we will need to run the workflow with the same python path that contains Anadama2. Everything else will then use Python 2.
```
pip install anadama2
git clone https://github.com/WillNickols/mitoribo
cd mitoribo
conda env create -f mitoribo.yml
conda activate mitoribo
```

## Databases
Create a folder called `refseq` and download the yeast genome GFF into it (such as [this one](http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz)). Split this file so that the appended chromosomes are their own fasta file (e.g., `tail -n 151990 saccharomyces_cerevisiae.gff > genome.fna`). Also, download RNA genes for decontamination (such as [this set](http://sgd-archive.yeastgenome.org/sequence/S288C_reference/rna/)).

## Running
Running the workflow requires:
1. The python path that contains anadama2 (the path is usually identified when running `pip install anadama2`)
2. An input folder with fastq or fastq.gz read files
3. An output folder
4. The adaptor sequence for the read library
5. A fasta file containing the yeast genome to align mRNA against
6. A fasta file with rRNA genes for removing non-mRNA reads

The workflow also allows extra arguments to TopHat2 as a string in quotes, an option to keep intermediate files, a minimum and maximum read length for calculating coverage, a chromosome for which to subset the coverage (the chromosome's yeast genome fasta header line up to the first space), the number of CPUs to use per task if the task allows multithreading, and the maximum number of tasks to run at once.
```
[PATH TO PYTHON3 WITH ANADAMA2] workflow.py \
  -i [FOLDER WITH FASTQS] \
  -o [OUTPUT FOLDER] \
  --adapter [ADAPTER SEQUENCE] \
  --genome [PATH TO YEAST GENOME] \
  --rna-coding [PATH TO RNA CODING FILE] \
  --input-extension [fastq/fastq.gz] \
  (--mitochrom [DEFAULT:None]) \
  (--tophat-args="[DEFAULT:None]") \
  (--keep-intermediates) \
  (--min-read-length [DEFAULT:23]) \
  (--max-read-length [DEFAULT:41]) \
  (--cores [DEFAULT:1]) \
  (--local-jobs [DEFAULT:1])
```

An example run is:
```
/home/czh/anaconda3/bin/python3 workflow.py \
  -i /media/czh/Zanlin_Chuankai001/MitoRibo_20230726/input/ \
  -o /media/czh/Zanlin_Chuankai001/MitoRibo_20230726/output/ \
  --adapter CTGTAGGCACCATCAAT \
  --genome /media/czh/Zanlin_Chuankai001/MitoRibo_20230726/refseq/genome.fna \
  --rna-coding /media/czh/Zanlin_Chuankai001/MitoRibo_20230726/refseq/rna_coding.fasta \
  --gene-file /media/czh/Zanlin_Chuankai001/MitoRibo_20230726/refseq/saccharomyces_cerevisiae.gff \
  --mitochrom chrmt \
  --input-extension fastq.gz \
  --cores 4 \
  --local-jobs 8
```

Based on [Couvillion et al. 2016](https://doi.org/10.1038/nature18015) and [Friedrich et al. 2021](https://doi.org/10.1016/j.celrep.2021.108711), the workflow does the following:
1. Trims read adapters with cutadapt
2. Removes the first 5' nucleotide because it frequently represents untemplated addition by reverse transcriptase Supercript III.
3. Removes reads that align to the non-tRNA RNA genes with bowtie1.
4. Removes reads that align to the tRNA genes with bowtie1 after trimming 3 nucleotides from the 3' end to account for non-templated CCA addition on mature tRNAs. (The 3' trimming is only temporary for this step.)
5. Aligns the cleaned reads to the yeast genome with TopHat2.
6. Reports read coverage over the chromosomes with samtools.

## Output
The output folder will contain the following items:
1. `anadama.log` which gives a log of the workflow commands run and any standard out or standard error output generated
2. `cleaned_reads` which contains a fastq for each input file once the reads are trimmed and decontaminated
3. `tophat_out` which contains one folder for each input file with the TopHat2 results
4. `coverage` which contains a `_coverage.txt` file for each input file with the cleaned read coverage of each position. This also contains a separate `_chromo_coverage.txt` file for the chromosome of interest if a chromosome of interest is specified.

## SRA Downloader
The `scripts/download_files.py` script can be used to download many fastq files in parallel like:
```
[PATH TO PYTHON3 WITH ANADAMA2] scripts/download_files.py -i to_download.txt -o [OUTPUT_FOLDER] --local-jobs [N]
```
where `to_download.txt` has one SRR identifier per line and N is the number of parallel downloads. If some downloads fail, you can run the script again, and already-downloaded files will be skipped.

Make sure you have downloaded and installed SRA tools according to the instructions [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit).
