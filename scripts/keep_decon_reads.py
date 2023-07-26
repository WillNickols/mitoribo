import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator #Biopython 1.51 or later

def main():
    parser = argparse.ArgumentParser(description="Keep only reads passing decontamination.")
    parser.add_argument("-i", "--input_file", help="Path to the input FASTQ file.")
    parser.add_argument("--decon", help="Path to the chopped tRNA decon file.")
    parser.add_argument("--min", help="Minimum read length.", default = 23)
    parser.add_argument("--max", help="Maximum read length.", default = 41)
    parser.add_argument("-o", "--output_file", help="Path to the output FASTQ file with modified headers.")
    args = parser.parse_args()

    new_file = open(args.output_file, "w")

    with open(args.decon, "rt") as handle:
        decon_set = set([title for title, seq, qual in FastqGeneralIterator(handle)])

    with open(args.input_file, "rt") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            if title in decon_set and len(seq) >= int(args.min) and len(seq) <= int(args.max):
                new_file.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                
    new_file.close()

if __name__ == "__main__":
    main()