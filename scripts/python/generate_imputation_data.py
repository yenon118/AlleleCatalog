#!/usr/bin/env python3

import sys
import os
import re
import argparse
import pathlib
import gzip
import threading

from joblib import Parallel, delayed


# Generate gff dictionary
# This function is for joblib parallel
def generate_gff_dictionary(line, gff_category, gff_key):
    gff_dictionary = {}
    if (line.startswith("#")) or (re.search(gff_category, line) is None) or \
            (re.search(gff_key, line) is None) or (re.search("scaffold", line) is not None):
        return None
    else:
        values = line.strip("\n").split("\t")
        chromosome = values[0]
        category = values[2]
        start = values[3]
        stop = values[4]
        key = ""
        key_array = values[8].strip().split(";")
        if category != gff_category:
            return None
        else:
            for i in range(len(key_array)):
                if str(key_array[i]).startswith(gff_key):
                    key = re.sub("(.*=)|(.*:)|(.*-)", "", re.sub(gff_key, "", str(key_array[i])))
                    break
            if (key != "") and (key not in gff_dictionary.keys()):
                gff_dictionary[key.upper()] = {
                    "Chromosome": chromosome,
                    "Category": category,
                    "Start": start,
                    "Stop": stop,
                    "Gene": key
                }
            elif (key != "") and (key in gff_dictionary.keys()):
                gff_dictionary[key.upper()]["Chromosome"] = chromosome
                gff_dictionary[key.upper()]["Category"] = category
                gff_dictionary[key.upper()]["Start"] = start
                gff_dictionary[key.upper()]["Stop"] = stop
                gff_dictionary[key.upper()]["Gene"] = key
            return gff_dictionary


# Generate imputation data for allele catalog
def generate_imputation_data(header, line, gff_dictionary, output_file_path, lock):
    # Parse header and line to get variant
    header_array = str(header).strip().split("\t")
    line_array = str(line).strip().split("\t")

    output_array = []

    chromosome = line_array[0]
    position = line_array[1]

    genes = []

    for key in gff_dictionary.keys():
        if (gff_dictionary[key]["Chromosome"] == chromosome) and (float(gff_dictionary[key]["Start"].strip()) <= float(position.strip()) <= float(gff_dictionary[key]["Stop"].strip())):
            if gff_dictionary[key]["Gene"] not in genes:
                genes.append(gff_dictionary[key]["Gene"])

    if len(genes) > 0:
        for gene in genes:
            for i in range(9, len(header_array)):
                genotype = re.sub(":.*", "", line_array[i])
                if re.fullmatch('(\./\.)|(\.\|\.)', genotype):
                    output_array.append(
                        header_array[i] + "\t" + gene + "\t" + chromosome + "\t" + position + "\t+\n" 
                    )
                else:
                    output_array.append(
                        header_array[i] + "\t" + gene + "\t" + chromosome + "\t" + position + "\t-\n" 
                    )
        if len(output_array) > 0:
            # Append the allele catalog string to the output file
            lock.acquire()
            with open(output_file_path, "a") as writer:
                writer.write(''.join(output_array))
            lock.release()


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_path = args.input_file
    gff_file_path = args.gff_file
    output_file_path = args.output_file

    n_jobs = args.threads
    gff_category = args.gff_category
    gff_key = args.gff_key

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
    # Read gff file then make a gff dictionary
    #######################################################################
    with open(gff_file_path, "r") as reader:
        results = Parallel(n_jobs=n_jobs)(
            delayed(generate_gff_dictionary)(line.strip(), gff_category, gff_key) for line in reader
        )
    gff_dictionary = {}
    for i in range(len(results)):
        if results[i] is not None:
            for key in results[i].keys():
                if key not in gff_dictionary.keys():
                    gff_dictionary[key] = results[i][key]

    #######################################################################
    # Write header to the output file
    #######################################################################
    writer = open(output_file_path, "w")
    writer.write(
        "Accession\tGene\tChromosome\tPosition\tImputation\n"
    )
    writer.close()

    #######################################################################
    # Read input file then generate allele catalog
    #######################################################################
    if str(input_file_path).endswith('gz'):
        with gzip.open(input_file_path, 'rt') as reader:
            header = ""
            while not header.strip().startswith("#CHROM"):
                header = reader.readline()
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(generate_imputation_data)(header, line, gff_dictionary, output_file_path, lock)
                for line in reader
            )
    else:
        with open(input_file_path, "r") as reader:
            header = ""
            while not header.strip().startswith("#CHROM"):
                header = reader.readline()
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(generate_imputation_data)(header, line, gff_dictionary, output_file_path, lock)
                for line in reader
            )

    #######################################################################
    # Create imputation dictionary
    #######################################################################
    imputation_dictionary = {}
    with open(output_file_path, "r") as reader:
        header = reader.readline()
        for line in reader:
            values = line.strip("\n").split("\t")
            accession = str(values[0])
            gene = str(values[1])
            imputation_notation = str(values[4])
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
    # Write header to the output file
    #######################################################################
    writer = open(output_file_path, "w")
    writer.write(
        "Accession\tGene\tImputation\n"
    )
    writer.close()

    #######################################################################
    # Write inputation data to the output file
    #######################################################################
    output_array = []
    for accession in imputation_dictionary.keys():
        for gene in imputation_dictionary[accession].keys():
            output_array.append(accession + "\t" + gene + "\t" + imputation_dictionary[accession][gene] + "\n")

            if len(output_array) > 10000:
                with open(output_file_path, "a") as writer:
                    writer.write(''.join(output_array))
                    output_array.clear()
    if len(output_array) != 0:
        with open(output_file_path, "a") as writer:
            writer.write(''.join(output_array))
            output_array.clear()


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='generate_imputation_data', description='generate imputation data')

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
    parser.add_argument('-g', '--gff_file', help='GFF file', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    parser.add_argument('-t', '--threads', help='Number of threads', type=int, default=10)
    parser.add_argument('-c', '--gff_category', help='Gff category', type=str, default='gene')
    parser.add_argument('-k', '--gff_key', help='Gff key', type=str, default='Name')

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)