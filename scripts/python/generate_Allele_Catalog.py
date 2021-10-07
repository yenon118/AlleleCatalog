#!/usr/bin/env python3

import sys
import os
import re
import argparse
import pathlib
import gzip
import threading

from joblib import Parallel, delayed

from Variant.Variant import Variant


# Generate reference dictionary
# This function is for joblib parallel
def generate_reference_dictionary(line):
    reference_dictionary = {}
    values = line.strip("\n").split("\t")
    if values[0] not in reference_dictionary.keys():
        reference_dictionary[values[0]] = {
            "Sample": values[0],
            "Maturity_Group": str(values[1]) if (
                    isinstance(values[1], int) or
                    str(values[1]).isdigit() or
                    isinstance(values[1], float) or
                    str(values[1]).replace('.', '', 1).isdigit()) else "",
            "Country": values[2],
            "State": values[3],
            "Improvement_Status": values[4],
            "Classification": values[5]
        }
    else:
        reference_dictionary[values[0]]["Sample"] = values[0]
        reference_dictionary[values[0]]["Maturity_Group"] = int(values[1])
        reference_dictionary[values[0]]["Country"] = values[2]
        reference_dictionary[values[0]]["State"] = values[3]
        reference_dictionary[values[0]]["Improvement_Status"] = values[4]
        reference_dictionary[values[0]]["Classification"] = values[5]
    return reference_dictionary


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


# Generate allele catalog
def generate_allele_catalog(header, line, reference_dictionary, gff_dictionary, output_file_path, lock):
    # Parse header and line to get variant
    variant = Variant(header, line)

    # Get chromosome and position
    chromosome = variant.get_chromosome()
    position = variant.get_position()

    # Fix gene because SnpEff gene is not always accurate
    gene_array = variant.get_gene().split("&")
    for g in gene_array:
        try:
            if g.upper() in gff_dictionary.keys():
                if (gff_dictionary[g.upper()]["Chromosome"] == chromosome) and \
                        (float(gff_dictionary[g.upper()]["Start"].strip()) <= float(position.strip()) <= float(
                            gff_dictionary[g.upper()]["Stop"].strip())):
                    variant.set_gene(g)
        except Exception as e:
            print("Allele Catalog - generate_allele_catalog function error !!!")
            print("Key " + g.upper() + " is not in the gff dictionary.")

    gene_array = variant.get_gene().split("&")
    if len(gene_array) > 1:
        for key in gff_dictionary.keys():
            if (gff_dictionary[key]["Chromosome"] == chromosome) and \
                    (float(gff_dictionary[key]["Start"].strip()) <= float(position.strip()) <= float(
                        gff_dictionary[key]["Stop"].strip())):
                variant.set_gene(gff_dictionary[key]["Gene"])

    gene_array = variant.get_gene().split("&")
    if len(gene_array) > 1:
        variant.set_gene("")

    # Get gene and genotypes dictionary
    gene = variant.get_gene()
    genotypes_dictionary = variant.get_genotypes_dictionary()

    # Append the allele catalog string to the output file
    lock.acquire()
    with open(output_file_path, "a") as writer:
        for key in genotypes_dictionary.keys():
            if (key in reference_dictionary.keys()) and (gene != "") and (gene is not None) and \
                    (genotypes_dictionary[key]["Genotype"] != "<INS>") and \
                    (genotypes_dictionary[key]["Genotype"] != "<DEL>"):
                writer.write(
                    reference_dictionary[key]["Classification"] + "\t" +
                    reference_dictionary[key]["Improvement_Status"] + "\t" +
                    reference_dictionary[key]["Maturity_Group"] + "\t" +
                    reference_dictionary[key]["Country"] + "\t" +
                    reference_dictionary[key]["State"] + "\t" +
                    key + "\t" +
                    chromosome + "\t" +
                    gene + "\t" +
                    position + "\t" +
                    genotypes_dictionary[key]["Genotype"] + "\t" +
                    genotypes_dictionary[key]["Annotated_Genotype"] + "\n"
                )
            elif (key not in reference_dictionary.keys()) and (gene != "") and (gene is not None) and \
                    (genotypes_dictionary[key]["Genotype"] != "<INS>") and \
                    (genotypes_dictionary[key]["Genotype"] != "<DEL>"):
                writer.write(
                    "" + "\t" +
                    "" + "\t" +
                    "" + "\t" +
                    "" + "\t" +
                    "" + "\t" +
                    key + "\t" +
                    chromosome + "\t" +
                    gene + "\t" +
                    position + "\t" +
                    genotypes_dictionary[key]["Genotype"] + "\t" +
                    genotypes_dictionary[key]["Annotated_Genotype"] + "\n"
                )
    lock.release()


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_path = args.input_file
    reference_file_path = args.reference_file
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
    # Read reference file then make a reference dictionary
    #######################################################################
    reference_dictionary = {}
    if reference_file_path is not None:
        with open(reference_file_path, "r") as reader:
            header = reader.readline()
            reference_column_names = header.strip().split("\t")
            # The results is a dictionary list, unwind the dictionary list and create dictionary with all samples.
            results = Parallel(n_jobs=n_jobs)(
                delayed(generate_reference_dictionary)(line) for line in reader
            )
        for i in range(len(results)):
            for key in results[i].keys():
                if key not in reference_dictionary.keys():
                    reference_dictionary[key] = results[i][key]

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
        "Classification\tImprovement_Status\tMaturity_Group\tCountry\tState\t" +
        "Accession\tChromosome\tGene\tPosition\tGenotype\tGenotype_with_Description\n"
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
                delayed(generate_allele_catalog)(header, line, reference_dictionary, gff_dictionary, output_file_path, lock)
                for line in reader
            )
    else:
        with open(input_file_path, "r") as reader:
            header = ""
            while not header.strip().startswith("#CHROM"):
                header = reader.readline()
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(generate_allele_catalog)(header, line, reference_dictionary, gff_dictionary, output_file_path, lock)
                for line in reader
            )


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='generate_Allele_Catalog', description='generate Allele Catalog')

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
    parser.add_argument('-g', '--gff_file', help='GFF file', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    parser.add_argument('-r', '--reference_file', help='Reference file', type=pathlib.Path, default=None)
    parser.add_argument('-t', '--threads', help='Number of threads', type=int, default=10)
    parser.add_argument('-c', '--gff_category', help='Gff category', type=str, default='gene')
    parser.add_argument('-k', '--gff_key', help='Gff key', type=str, default='Name')

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
