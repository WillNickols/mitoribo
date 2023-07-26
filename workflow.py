#!/usr/bin/env python3

from anadama2 import Workflow
import glob
import os
import math
import shutil
from anadama2.tracked import TrackedDirectory
from pathlib import Path

# Parse arguments
workflow = Workflow(version="0.1", description="MAG and SGB workflow")
workflow.add_argument("adapter", desc="Adapter sequence to trim")
workflow.add_argument("refs", desc="File directory with genome.fasta and rna_coding.fasta")
workflow.add_argument("mitochrom", desc="Name of mitochondrion chromosome in reference sequence")
workflow.add_argument("input-extension", desc="The input file extension", default="fastq.gz")
workflow.add_argument("tophat-args", desc="Extra tophat arguments", default="")
workflow.add_argument("keep-intermediates", desc="Don't delete intermediate files", action="store_true")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
args = workflow.parse_args()

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

scripts_folder = this_folder + "scripts/"

# list the input fastq files
in_dir = args.input
output = os.path.abspath(args.output.rstrip("/")) + "/"
adapter_seq = args.adapter
mitochrom = args.mitochrom
input_extension = args.input_extension
ref_seq_dir = args.refs
tophat_args = args.tophat_args

if input_extension not in ["fastq", "fastq.gz"]:
	raise ValueError("Input extension must be fastq or fastq.gz")

# Get filepath without paired/unpaired ending
paths = glob.glob(os.path.abspath(in_dir.rstrip("/")) + "/" + '*.' + input_extension)
names = set((file.split("." + input_extension)[0]).split("/")[-1] for file in paths)

if len(paths) == 0:
	raise ValueError("No input files")

tmp_dir = output + "tmp/"
os.makedirs(tmp_dir, exist_ok=True)

gunzip_dir = tmp_dir + "gunzip_dir/"
os.makedirs(gunzip_dir, exist_ok=True)

decon_dir = tmp_dir + "decon/"
os.makedirs(decon_dir, exist_ok=True)

non_tRNA = decon_dir + "non_tRNA/"
os.makedirs(non_tRNA, exist_ok=True)

tRNA = decon_dir + "tRNA/"
os.makedirs(tRNA, exist_ok=True)

genome_bt = tmp_dir + "genome_bt/"
os.makedirs(genome_bt, exist_ok=True)

cut_reads = tmp_dir + "cut_reads/"
os.makedirs(cut_reads, exist_ok=True)

cut_reads_5 = tmp_dir + "5_cut_reads/"
os.makedirs(cut_reads_5, exist_ok=True)

cleaned_reads = output + "cleaned_reads/"
os.makedirs(cleaned_reads, exist_ok=True)

decon_nontRNA = tmp_dir + "decon_nontRNA/"
os.makedirs(decon_nontRNA, exist_ok=True)

decon_tRNA = tmp_dir + "decon_tRNA/"
os.makedirs(decon_tRNA, exist_ok=True)

last_3_nucs_gone = tmp_dir + "last_3_nucs_gone/"
os.makedirs(last_3_nucs_gone, exist_ok=True)

tophat_out = output + "tophat_out/"
os.makedirs(tophat_out, exist_ok=True)

cores = args.cores
local_jobs = args.jobs

#################################
# function to list dependencies #
#################################

def list_depends(name, step):
	if step == "gunzip":
		depends_list = [os.path.abspath(in_dir.rstrip("/")) + "/" + name + "." + input_extension]
	if step == "cutadapt":
		if input_extension == "fastq":
			depends_list = [os.path.abspath(in_dir.rstrip("/")) + "/" + name + ".fastq"]
		else:
			depends_list = [gunzip_dir + name + ".fastq"]
		return depends_list
	elif step == "cut_5_reads":
		depends_list = [cut_reads + name + ".fastq"]
	elif step == "split_tRNAs":
		depends_list = [ref_seq_dir + "rna_coding.fasta"]
	elif step == "decon_non_tRNA":
		depends_list = [cut_reads_5 + name + ".fastq",
		decon_dir + "nontRNA.1.ebwt", 
		decon_dir + "nontRNA.2.ebwt", 
		decon_dir + "nontRNA.3.ebwt", 
		decon_dir + "nontRNA.4.ebwt", 
		decon_dir + "nontRNA.rev.1.ebwt", 
		decon_dir + "nontRNA.rev.2.ebwt"]
	elif step == "last_3_nucs_gone":
		depends_list = [non_tRNA + name + ".fastq"]
	elif step == "decon_tRNA":
		depends_list = [last_3_nucs_gone + name + ".fastq",
		decon_dir + "tRNA.1.ebwt", 
		decon_dir + "tRNA.2.ebwt", 
		decon_dir + "tRNA.3.ebwt", 
		decon_dir + "tRNA.4.ebwt", 
		decon_dir + "tRNA.rev.1.ebwt", 
		decon_dir + "tRNA.rev.2.ebwt"]
	elif step == "subset_reads":
		depends_list = [non_tRNA + name + ".fastq", tRNA + name + ".fastq"]
	elif step == "build_ref_bowtie_db":
		depends_list = [ref_seq_dir + "genome.fasta"]
	elif step == "tophat":
		depends_list = [cleaned_reads + name + ".fastq",
		genome_bt + "full_genome.1.bt2", 
		genome_bt + "full_genome.2.bt2", 
		genome_bt + "full_genome.3.bt2", 
		genome_bt + "full_genome.4.bt2", 
		genome_bt + "full_genome.rev.1.bt2", 
		genome_bt + "full_genome.rev.2.bt2"]
	elif step == "samtools":
		depends_list = [tophat_out + name + "/accepted_hits.bam"]
	elif step == "get_chromo_coverage":
		depends_list = [output + name + "/coverage.txt"]
	elif step == "delete_tmp":
		depends_list = [tophat_out + name + "/accepted_hits.bam"]
	else:
		raise ValueError("Invalid step")
	return depends_list

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if step == "gunzip":
		target_list = [gunzip_dir + name + ".fastq"]
	if step == "cutadapt":
		target_list = [cut_reads + name + ".fastq"]
	elif step == "cut_5_reads":
		target_list = [cut_reads_5 + name + ".fastq"]
	elif step == "split_tRNAs":
		target_list = [decon_dir + "nontRNA.1.ebwt", 
		decon_dir + "nontRNA.2.ebwt", 
		decon_dir + "nontRNA.3.ebwt", 
		decon_dir + "nontRNA.4.ebwt", 
		decon_dir + "nontRNA.rev.1.ebwt", 
		decon_dir + "nontRNA.rev.2.ebwt",
		decon_dir + "tRNA.1.ebwt", 
		decon_dir + "tRNA.2.ebwt", 
		decon_dir + "tRNA.3.ebwt", 
		decon_dir + "tRNA.4.ebwt", 
		decon_dir + "tRNA.rev.1.ebwt", 
		decon_dir + "tRNA.rev.2.ebwt"]
	elif step == "decon_non_tRNA":
		target_list = [non_tRNA + name + ".fastq"]
	elif step == "last_3_nucs_gone":
		target_list = [last_3_nucs_gone + name + ".fastq"]
	elif step == "decon_tRNA":
		target_list = [tRNA + name + ".fastq"]
	elif step == "subset_reads":
		target_list = [cleaned_reads + name + ".fastq"]
	elif step == "build_ref_bowtie_db":
		target_list = [genome_bt + "full_genome.1.bt2", 
		genome_bt + "full_genome.2.bt2", 
		genome_bt + "full_genome.3.bt2", 
		genome_bt + "full_genome.4.bt2", 
		genome_bt + "full_genome.rev.1.bt2", 
		genome_bt + "full_genome.rev.2.bt2"]
	elif step == "tophat":
		target_list = [tophat_out + name + "/accepted_hits.bam"]
	elif step == "samtools":
		target_list = [output + name + "/coverage.txt"]
	elif step == "get_chromo_coverage":
		target_list = [output + name + "/chromo_coverage.txt"]
	else:
		raise ValueError("Invalid step")
	return target_list

############################
# unzip files if necessary #
############################

def gunzip(name):
	command = "gunzip -c [depends[0]] > [targets[0]]"
	return str(command)

if input_extension == "fastq.gz":
	for name in names:
		workflow.add_task(actions=gunzip(name),
			depends=list_depends(name=name, step="gunzip"),
			targets=list_targets(name=name, step="gunzip"),
			name="Decompress " + name.split("/")[-1]
			)

################
# run cutadapt #
################

def cutadapt(name):
	command = "cutadapt -a " + adapter_seq + " -o [targets[0]] [depends[0]] --cores=" + str(cores)
	return str(command)

for name in names:
	workflow.add_task(actions=cutadapt(name),
		depends=list_depends(name=name, step="cutadapt"),
		targets=list_targets(name=name, step="cutadapt"),
		name="Cutadapt for " + name.split("/")[-1]
		)

###################
# cut 5-prime end #
###################

def cut_5_reads(name):
	command = "python " + scripts_folder + "chop_first_nuc.py -i [depends[0]] -o [targets[0]]"
	return str(command)

for name in names:
	workflow.add_task(actions=cut_5_reads(name),
		depends=list_depends(name=name, step="cut_5_reads"),
		targets=list_targets(name=name, step="cut_5_reads"),
		name="Cut 1 nt on 5-prime end for " + name.split("/")[-1]
		)

#################################
# build decontamination indexes #
#################################

def split_tRNAs(name):
	command = '''{a} && {b} && {c}'''.format(
		a = "python " + scripts_folder + "split_tRNAs.py -i [depends[0]] -o " + decon_dir + "/rna_coding",
		b = "bowtie-build " + decon_dir + "/rna_coding_nontRNA.fasta " + decon_dir + "/nontRNA",
		c = "bowtie-build " + decon_dir + "/rna_coding_tRNA.fasta " + decon_dir + "/tRNA"
		)

	return str(command)

workflow.add_task(actions=split_tRNAs(name),
	depends=list_depends(name=name, step="split_tRNAs"),
	targets=list_targets(name=name, step="split_tRNAs"),
	name="Build decontamination indexes"
	)

##########################
# decontaminate non-tRNA #
##########################

def decon_non_tRNA(name):
	command = "bowtie --threads " + str(cores) + " --sam " + decon_dir + "/nontRNA [depends[0]] /dev/null --un [targets[0]]"
	return str(command)

for name in names:
	workflow.add_task(actions=decon_non_tRNA(name),
		depends=list_depends(name=name, step="decon_non_tRNA"),
		targets=list_targets(name=name, step="decon_non_tRNA"),
		name="Decontaminate the non-tRNA for " + name.split("/")[-1]
		)

##############################
# remove 3 nt for tRNA decon #
##############################

def last_3_nucs_gone(name):
	command = "python " + scripts_folder + "chop_last_3_nuc.py -i [depends[0]] -o [targets[0]]"
	return str(command)

for name in names:
	workflow.add_task(actions=last_3_nucs_gone(name),
		depends=list_depends(name=name, step="last_3_nucs_gone"),
		targets=list_targets(name=name, step="last_3_nucs_gone"),
		name="Remove last 3 nucleotides for decontaminating " + name.split("/")[-1]
		)

######################
# decontaminate tRNA #
######################

def decon_tRNA(name):
	command = "bowtie --threads " + str(cores) + " --sam " + decon_dir + "/nontRNA [depends[0]] /dev/null --un [targets[0]]"
	return str(command)

for name in names:
	workflow.add_task(actions=decon_tRNA(name),
		depends=list_depends(name=name, step="decon_tRNA"),
		targets=list_targets(name=name, step="decon_tRNA"),
		name="Decontaminate the tRNA for " + name.split("/")[-1]
		)

################################################
# subset reads to only the decontaminated ones #
################################################

def subset_reads(name):
	command = "python " + scripts_folder + "keep_decon_reads.py -i [depends[0]] --decon [depends[1]] -o [targets[0]]"
	return str(command)

for name in names:
	workflow.add_task(actions=subset_reads(name),
		depends=list_depends(name=name, step="subset_reads"),
		targets=list_targets(name=name, step="subset_reads"),
		name="Finish decontaminating for " + name.split("/")[-1]
		)

################################
# build tophat reference index #
################################

def build_ref_bowtie_db(name):
	command = "bowtie2-build [depends[0]] " + genome_bt + "full_genome"
	return str(command)

workflow.add_task(actions=build_ref_bowtie_db(name),
	depends=list_depends(name=name, step="build_ref_bowtie_db"),
	targets=list_targets(name=name, step="build_ref_bowtie_db"),
	name="Build tophat reference index"
	)

##############
# run tophat #
##############

def tophat(name):
	command = "tophat -p " + str(cores) + " " + str(tophat_args) + " -o " + tophat_out + name + " " + genome_bt + "full_genome [depends[0]]"
	return str(command)

for name in names:
	workflow.add_task(actions=tophat(name),
		depends=list_depends(name=name, step="tophat"),
		targets=list_targets(name=name, step="tophat"),
		name="Run tophat for " + name.split("/")[-1]
		)

################
# run samtools #
################

def samtools(name):
	command = '''{a} && {b}'''.format(
		a = "samtools index [depends[0]]",
		b = "samtools depth -a [depends[0]] > [targets[0]]"
		)

	return str(command)

for name in names:
	workflow.add_task(actions=samtools(name),
		depends=list_depends(name=name, step="samtools"),
		targets=list_targets(name=name, step="samtools"),
		name="Get coverage depth for " + name.split("/")[-1]
		)

###########################################
# get coverage of a particular chromosome #
###########################################

def get_chromo_coverage(name):
	command = "grep " + mitochrom + " [depends[0]] > [targets[0]]"
	return str(command)

for name in names:
	workflow.add_task(actions=get_chromo_coverage(name),
		depends=list_depends(name=name, step="get_chromo_coverage"),
		targets=list_targets(name=name, step="get_chromo_coverage"),
		name="Get chromosome coverage depth for " + name.split("/")[-1]
		)

##########################
# remove temporary files #
##########################

if not args.keep_intermediates:
	workflow.add_task(actions="rm -r " + tmp_dir,
		depends=list_depends(name=name, step="delete_tmp"),
		name="Delete temporary files"
		)

####################
# run the workflow #
####################

workflow.go()