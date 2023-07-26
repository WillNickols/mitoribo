import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator #Biopython 1.51 or later

def main():
    parser = argparse.ArgumentParser(description="Remove the first nucleotide from each sequence in a FASTQ file.")
    parser.add_argument("-i", "--input_file", help="Path to the input FASTQ file.")
    parser.add_argument("-o", "--output_file", help="Path to the output FASTQ file with modified headers.")
    args = parser.parse_args()
    
    new_file = open(args.output_file, "w")

    with open(args.input_file, "rt") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            if len(seq) > 3:
                seq = seq[:-3]
                qual = qual[:-3]
                new_file.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

    new_file.close()

if __name__ == "__main__":
    main()
