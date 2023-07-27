import pysam
import csv
import argparse

def main():
    parser = argparse.ArgumentParser(description="Get the length distribution of reads for alignment.")
    parser.add_argument("-i", "--input_file", help="Path to the input bam file.")
    parser.add_argument("-o", "--output_file", help="Path to the output csv file.")
    args = parser.parse_args()

    read_lengths = []
    with pysam.AlignmentFile(args.input_file, "rb") as bam:
        for read in bam.fetch():
            read_lengths.append(read.query_length)

    with open(args.output_file, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        for item in read_lengths:
            csv_writer.writerow([item])
                

if __name__ == "__main__":
    main()
