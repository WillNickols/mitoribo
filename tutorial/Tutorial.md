## Installation and downloads
If you do not already have it installed, install Conda [according to its instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). 

Open the terminal and navigate to a directory from which you'll run the tutorial. Make sure there are no spaces in the path to the directory (this will cause issues with Tophat).
```
cd YOUR_DIRECTORY
```

The following commands will install Anadama2 (a workflow management library) and this GitHub repository to run the workflow. The GitHub repository contains a `.yml` file to create a conda environment with all the necessary tools. Because TopHat2 requires Python 2 and Anadama2 requires Python 3, we will need to run the workflow with the same python path that contains Anadama2. Everything else will then use Python 2.
```
pip install anadama2
```
After you run this, you'll see a message showing where anadama2 was installed (something like `/Users/williamnickols/opt/anaconda3/lib/python3.9/site-packages`). Copy this location to a word document for temporary storage.

Continue installing the rest of the software:
```
git clone https://github.com/WillNickols/mitoribo
cd mitoribo
conda env create -f mitoribo.yml
conda activate mitoribo
```

If you open a file viewer, you will now see `mitoribo` as a folder in the directory you chose at first. In the terminal, create a new subdirectory for downloading the yeast genomes and move into it:
```
mkdir genomes
cd genomes
```

Now, we will install the necessary databases.

1. Download the `saccharomyces_cerevisiae.gff.gz` file from [the SGD database](http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/). This contains the yeast genes and genome. Once downloaded into the `genomes` folder, unzip this with
```
gunzip saccharomyces_cerevisiae.gff.gz
```

2. Split the genes file so that the appended chromosomes are their own genomic fasta file:

```
grep -n ">chrI$" saccharomyces_cerevisiae.gff | cut -d':' -f1 | xargs -I {} tail -n +{} saccharomyces_cerevisiae.gff > genome.fna
```

3. Download the `rna_coding.fasta.gz` file from [the SGD database](http://sgd-archive.yeastgenome.org/sequence/S288C_reference/rna/). This contains RNA genes for decontamination. Once downloaded into the `genomes` folder, unzip this with
```
gunzip rna_coding.fasta.gz
```

Finally, we need to download the sequencing files we'll be using. Back up to the `mitoribo` folder, create a directory called `inputs`, and place these four files in it (they can stay gzipped): [SRR2889825](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2889825&display=download), [SRR2889829](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2889829&display=download), [SRR2889830](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2889830&display=download), [SRR2889834](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2889834&display=download).
```
cd ..
mkdir inputs
```

## Running
Running the workflow requires:
1. The python path that contains anadama2 (if your anadama2 installation was in `../anaconda3/lib/python3.9/site-packages`, this will be something like `../anaconda3/bin/python3.9`)
2. An input folder with fastq or fastq.gz read files
3. An output folder
4. The adaptor sequence for the read library
5. A fasta file containing the yeast genome to align mRNA against
6. A fasta file with rRNA genes for removing non-mRNA reads

There are also many optional arguments, and the relevant ones for us are `--keep-intermediates` which keeps intermediate files, `--mitochrom` which is the mitochondrial chromosome's fasta header line, `--cores` which is the number of CPUs to use per task, and `--local-jobs` which is the maximum number of tasks to run simultaneously. On MacOS, the `--local-jobs` currently must be set to 1, but on a Linux machine, it can (and should) be set to the number of CPUs available.

We can now run the workflow with the following command:

```
<PATH TO PYTHON3 WITH ANADAMA2>python3 workflow.py \
  -i inputs \
  -o outputs \
  --adapter CTGTAGGCACCATCAAT \
  --genome genomes/genome.fna \
  --rna-coding genomes/rna_coding.fasta \
  --gene-file genomes/saccharomyces_cerevisiae.gff \
  --mitochrom chrmt \
  --input-extension fastq.gz \
  --cores 4 \
  --local-jobs 1 \
  --keep-intermediates
```

## Output
The `outputs` folder will contain the following items:
1. `anadama.log` which gives a log of the workflow commands run and any standard out or standard error output generated
2. `cleaned_reads` which contains a fastq for each input file once the reads are trimmed and decontaminated
3. `tophat_out` which contains one folder for each input file with the TopHat2 results
4. `coverage` which contains a `_coverage.txt` file for each input file with the cleaned read coverage of each position. This also contains a separate `_chromo_coverage.txt` file for the chromosome of interest if a chromosome of interest is specified.

## Analysis
The `analysis.py` script in the tutorials folder will run a simple analysis for this set of ribosome sequencing data. It takes in the outputs folder from the workflow we ran and produces plots of gene coverage in the `tutorial/coverage/` directory.
```
python tutorial/analysis.py --input outputs --output tutorial/coverage/
```






