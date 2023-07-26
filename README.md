# mitoribo
Bioinformatics workflow for ribosome sequencing of yeast

## Installation
```
pip install anadama2
git clone https://github.com/WillNickols/mitoribo
cd mitoribo
conda env create -f mitoribo.yml
conda activate mitoribo
```

## Databases
Create a folder called `refseqs` and download a yeast genome into it (such as [this one](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/)). Also, download RNA genes for decontamination (such as [this set](http://sgd-archive.yeastgenome.org/sequence/S288C_reference/rna/)).

## Running
```
[PATH TO PYTHON3 WITH ANADAMA2] workflow.py \
  -i [FOLDER WITH FASTQS] \
  -o [OUTPUT FOLDER] \
  --adapter [ADAPTER SEQUENCE] \
  --genome [PATH TO YEAST GENOME] \
  --rna-coding [PATH TO RNA CODING FILE] \
  --mitochrom [FIRST LINE UP TO A SPACE OF MITOCHONDRION FASTA ENTRY] \
  --input-extension [fastq/fastq.gz] \
  (--tophat-args="[EXTRA TOPHAT ARGUMENTS]") \
  (--keep-intermediates) \
  (--min-read-length) \
  (--max-read-length) \
  (--cores [N]) \
  (--local-jobs [N])
```

For example,
```
/home/czh/anaconda3/bin/python3 workflow.py \
  -i /media/czh/Zanlin_Chuankai001/MitoRibo_20230726/input/ \
  -o /media/czh/Zanlin_Chuankai001/MitoRibo_20230726/output/ \
  --adapter CTGTAGGCACCATCAAT \
  --genome /media/czh/Zanlin_Chuankai001/MitoRibo_20230726/refseq/GCF_000146045.2_R64_genomic.fna \
  --rna-coding /media/czh/Zanlin_Chuankai001/MitoRibo_20230726/refseq/rna_coding.fasta \
  --mitochrom NC_001224.1 \
  --input-extension fastq \
  --keep-intermediates \
  --cores 4 \
  --local-jobs 8
```
