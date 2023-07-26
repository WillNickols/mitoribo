from Bio import SeqIO
import argparse

def separate_sequences(input_file, trna_output_file, non_trna_output_file):
    with open(trna_output_file, 'w') as f_trna, open(non_trna_output_file, 'w') as f_non_trna:
        for record in SeqIO.parse(input_file, "fasta"):
            if "tRNA_gene" in record.description:
                SeqIO.write(record, f_trna, "fasta")
            else:
                SeqIO.write(record, f_non_trna, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove the first nucleotide from each sequence in a FASTQ file.")
    parser.add_argument("-i", "--input_file", help="Path to the input FASTQ file.")
    parser.add_argument("-o", "--output_file", help="Path to the output FASTQ file with modified headers.")
    args = parser.parse_args()

    separate_sequences(args.input_file, args.output_file + "_tRNA.fasta", args.output_file + "_nontRNA.fasta")
