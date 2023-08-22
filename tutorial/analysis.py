#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import glob
import os
import numpy as np
import argparse

# Read arguments
parser = argparse.ArgumentParser(description="Analysis of glucose/glycerol mitoribosome sequencing")
parser.add_argument("-i", "--input", help="Path to workflow output directory")
parser.add_argument("-o", "--output", help="Path to new output directory")
args = parser.parse_args()

input_dir = os.path.abspath(args.input.rstrip("/")) + "/"
output_dir = os.path.abspath(args.output.rstrip("/")) + "/"

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

# Read in glycerol or glucose
this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"
dataset = pd.read_csv(this_folder + "datasets.csv")

# Read in coverage files
files = glob.glob(os.path.join(input_dir, "coverage/*chromo_coverage.txt"))
tmp_dict = {file.split('/')[-1].split('_')[0]: file for file in files}
df = pd.DataFrame(tmp_dict.items(), columns=['SRR', 'file'])
dataset = pd.merge(dataset, df, on='SRR')
conditions = dataset['Condition'][dataset['study'] == 'PRJNA300880']

# Create dictionaries to store coverage data for each condition
condition_coverage = {condition: [] for condition in set(conditions)}

for index, row in dataset.iterrows():
    file = row['file']
    condition = row['Condition']

    # Read the coverage data from the file
    positions = []
    coverage = []
    with open(file, "r") as f:
        for line in f:
            chrom, pos, cov = line.strip().split("\t")
            positions.append(int(pos))
            coverage.append(int(cov))

    coverage = np.array(coverage).astype(float) / np.sum(coverage) * 1000000

    condition_coverage[condition].append(coverage)

# Plot the coverage data for each condition
plt.figure(figsize=(10, 6))
colors = ["blue", "red"]
for condition, color in zip(set(conditions), colors):
    mean_coverage = np.mean(condition_coverage[condition], axis=0)
    for coverage in condition_coverage[condition]:
        plt.plot(positions, coverage, color=color, alpha=0.2)
    plt.plot(positions, mean_coverage, color=color, label=condition)
plt.xlabel("Genomic Position")
plt.ylabel("Coverage (depth per million mitochondrial reads)")
plt.title("Whole Chromosome Coverage Map")
plt.grid(True)
plt.legend()
plt.savefig(output_dir + 'whole_chromosome_coverage.png', dpi=500)

genes_and_positions = pd.read_csv(this_folder + "mt_genes.csv")

for gene_index, gene_row in genes_and_positions.iterrows():
    # Create dictionaries to store coverage data for each condition
    condition_coverage = {condition: [] for condition in set(conditions)}

    n_panels = 1
    if ',' in gene_row['position']:
        n_panels = len(gene_row['position'].split(','))
    start_pos = int(gene_row['position'].split("-")[0])
    end_pos = int(gene_row['position'].split("-")[-1])
    for index, row in dataset.iterrows():
        file = row['file']
        condition = row['Condition']
        
        # Read the coverage data from the file
        positions = []
        coverage = []
        with open(file, "r") as f:
            for line in f:
                chrom, pos, cov = line.strip().split("\t")
                positions.append(int(pos))
                coverage.append(int(cov))

        coverage = np.array(coverage).astype(float) / np.sum(coverage) * 1000000
        positions = positions[start_pos:end_pos]
        coverage = coverage[start_pos:end_pos]
        
        condition_coverage[condition].append(coverage)
        
    # Plot the coverage data for each condition
    colors = ["blue", "red"]
    
    if n_panels > 1:
        f, axs = plt.subplots(1, n_panels, sharey=True, facecolor='w', figsize=(10, 6))
        plt.suptitle(gene_row['gene'] + " Coverage Map")
        
        splits = [(int(tmp.split('-')[0]), int(tmp.split('-')[1])) for tmp in gene_row['position'].split(',')]
        for i, axes_lim in enumerate(splits):
            for condition, color in zip(set(conditions), colors):
                mean_coverage = np.mean(condition_coverage[condition], axis=0)
                for coverage in condition_coverage[condition]:
                    axs[i].plot(positions, coverage, color=color, alpha=0.2)
                axs[i].plot(positions, mean_coverage, color=color, label=condition)
            axs[i].set_xlim(axes_lim[0],axes_lim[1])

            # hide the spines between ax and ax2
            if i < len(splits) - 1:
                axs[i].spines['right'].set_visible(False)
                axs[i + 1].spines['left'].set_visible(False)
                axs[i + 1].yaxis.tick_right()
            axs[i].grid(True)
            axs[i].tick_params(right = False)
            axs[i].tick_params(axis="x",labelrotation=90)

        axs[0].set_ylabel("Coverage (depth per million mitochondrial reads)")
        plt.legend()
        plt.savefig(output_dir + gene_row['gene'] + '_coverage.png', dpi=500)
        
    else:
        plt.figure(figsize=(10, 6))
        for condition, color in zip(set(conditions), colors):
            mean_coverage = np.mean(condition_coverage[condition], axis=0)
            for coverage in condition_coverage[condition]:
                plt.plot(positions, coverage, color=color, alpha=0.2)
            plt.plot(positions, mean_coverage, color=color, label=condition)
        
        plt.legend()
        plt.xlabel("Genomic Position")
        plt.ylabel("Coverage (depth per million mitochondrial reads)")
        plt.title(gene_row['gene'] + " Coverage Map")
        plt.grid(True)
        plt.xticks(rotation=90)
        plt.savefig(output_dir + gene_row['gene'] + '_coverage.png', dpi=500)
