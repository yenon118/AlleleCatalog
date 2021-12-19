#!/usr/bin/env python3

import sys
import os
import re
import argparse
import pathlib
import gzip
import threading

from joblib import Parallel, delayed


# Generate allele catalog with imputation data
def combine_imputation_data(line, imputation_dictionary, output_file_path, lock):
    # Parse line to get variant
    line_array = str(line).strip("\n").split("\t")

    accession = str(line_array[5])
    gene = str(line_array[7])
    imputation_notation = "-"

    if accession in imputation_dictionary.keys():
        if gene in imputation_dictionary[accession].keys():
            imputation_notation = imputation_dictionary[accession][gene]

    lock.acquire()
    with open(output_file_path, "a") as writer:
        writer.write(line.strip("\n") + "\t" + imputation_notation + "\n")
    lock.release()


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_path = args.input_file
    reference_file_path = args.reference_file
    output_file_path = args.output_file

    n_jobs = args.threads

    #######################################################################
    # Create a threading lock
    #######################################################################
    lock = threading.Lock()

    #######################################################################
    # Check if output parent folder exists
    # If not, create the output parent folder
    #######################################################################
    if not output_file_path.parent.exists():
        try:
            output_file_path.parent.mkdir(parents=True)
        except FileNotFoundError as e:
            pass
        except FileExistsError as e:
            pass
        except Exception as e:
            pass
        if not output_file_path.parent.exists():
            sys.exit(1)

    #######################################################################
    # Read reference imputation file then make a imputation dictionary
    #######################################################################
    imputation_dictionary = {}
    with open(reference_file_path, "r") as reader:
        header = reader.readline()
        for line in reader:
            values = line.strip("\n").split("\t")
            accession = str(values[0])
            gene = str(values[1])
            imputation_notation = str(values[2])
            if accession in imputation_dictionary.keys():
                if gene in imputation_dictionary[accession].keys():
                    if imputation_dictionary[accession][gene] == "-" and imputation_notation == "+":
                        imputation_dictionary[accession][gene] = imputation_notation
                else:
                    imputation_dictionary[accession][gene] = imputation_notation
            else:
                imputation_dictionary[accession] = {}
                imputation_dictionary[accession][gene] = imputation_notation

    #######################################################################
    # Read input file then generate allele catalog
    #######################################################################
    if str(input_file_path).endswith('gz'):
        with gzip.open(input_file_path, 'rt') as reader:
            header = ""
            while not header.strip().startswith("Classification"):
                header = reader.readline()
            if header != "":
                with open(output_file_path, "w") as writer:
                    writer.write(header.strip() + "\t" + "Imputation" + "\n")
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(combine_imputation_data)(line, imputation_dictionary, output_file_path, lock)
                for line in reader
            )
    else:
        with open(input_file_path, "r") as reader:
            header = ""
            while not header.strip().startswith("Classification"):
                header = reader.readline()
            if header != "":
                with open(output_file_path, "w") as writer:
                    writer.write(header.strip() + "\t" + "Imputation" + "\n")
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(combine_imputation_data)(line, imputation_dictionary, output_file_path, lock)
                for line in reader
            )


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='combine_imputation_data', description='combine imputation data')

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
    parser.add_argument('-r', '--reference_file', help='Reference file', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    parser.add_argument('-t', '--threads', help='Number of threads', type=int, default=10)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)