#!/usr/bin/python

import sys
import argparse

def parse_arguments():
    """Read arguments from the console"""
    parser = argparse.ArgumentParser(description="Note: A simple conversiton of fastq file to fasta.")
    parser.add_argument("-i", "--input", help='fastq file as input')
    parser.add_argument("-o", "--output", help='fasta file as output')

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    with open(args.input, "r") as fq:
        with open(args.output, "w") as fa:
            for line_num, line in enumerate(fq):
                if (line_num % 4 ==0):
                    fa.write(line.replace("@", ">").rstrip() + "\n")
                if (line_num % 4 ==1):
                    fa.write(line.rstrip() + "\n")

if __name__ == "__main__":
    main()