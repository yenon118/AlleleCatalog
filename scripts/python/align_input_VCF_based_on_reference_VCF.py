#!/usr/bin/env python3

import sys
import os
import re
import argparse
import pathlib

from joblib import Parallel, delayed

import pandas as pd


def process_reference_file(header, line, reference_output_file_path):
    line = str(line).strip()
    line_array = line.split("\t")
    line_array[2] = "."
    line_array[5] = "."
    line_array[6] = "PASS"
    line_array[7] = "."
    line_array[8] = "GT"
    for i in range(9, len(line_array)):
        genotype_indexes = re.sub("(:.*)", "", line_array[i])
        genotype_indexes_array = genotype_indexes.split("/")
        if len(genotype_indexes_array) == 2:
            if genotype_indexes_array[0] != genotype_indexes_array[1]:
                line_array[i] = "./."
            else:
                line_array[i] = genotype_indexes
        else:
            line_array[i] = "./."
    with open(reference_output_file_path, 'a') as writer:
        writer.write(str('\t'.join(line_array)) + '\n')
    return [line_array[0], line_array[1], line_array[3], line_array[4]]


def align_input_file_based_on_reference_dictionary(header, line, chrom_pos_ref_alt_dict, input_output_file_path):
    line = str(line).strip()
    line_array = line.split("\t")
    line_array[2] = "."
    line_array[5] = "."
    line_array[6] = "PASS"
    line_array[7] = "."
    line_array[8] = "GT"
    if line_array[0] in chrom_pos_ref_alt_dict.keys():
        if line_array[1] in chrom_pos_ref_alt_dict[line_array[0]].keys():
            for i in range(9, len(line_array)):
                genotype_indexes = re.sub(":.*", "", line_array[i])
                genotype_indexes_array = genotype_indexes.split("/")
                if len(genotype_indexes_array) == 2:
                    if genotype_indexes_array[0] != genotype_indexes_array[1]:
                        line_array[i] = "./."
                    else:
                        line_array[i] = genotype_indexes
                else:
                    line_array[i] = "./."
            with open(input_output_file_path, 'a') as writer:
                writer.write(str('\t'.join(line_array)) + '\n')
            return [line_array[0], line_array[1], line_array[3], line_array[4]]
        else:
            return None
    else:
        return None


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_path = args.input_file
    reference_file_path = args.reference_file
    output_folder_path = args.output_folder

    n_jobs = args.threads

    #######################################################################
    # Check if output parent folder exists
    # If not, create the output parent folder
    #######################################################################
    if not output_folder_path.exists():
        try:
            output_folder_path.mkdir(parents=True)
        except FileNotFoundError as e:
            pass
        except FileExistsError as e:
            pass
        except Exception as e:
            pass
        if not output_folder_path.exists():
            sys.exit(1)

    #######################################################################
    # Make sure output files exist
    #######################################################################
    with open(output_folder_path.joinpath(input_file_path.name), 'w') as writer:
        writer.write("")
    with open(output_folder_path.joinpath(reference_file_path.name), 'w') as writer:
        writer.write("")

    #######################################################################
    # Collect chrom, pos, ref, and alt from reference file
    # Process reference file
    #######################################################################
    reference_metadata_array = []
    header = ""
    with open(reference_file_path, 'r') as reader:
        # Get metadata and header
        while not header.strip().startswith("#CHROM"):
            header = reader.readline()
            if header.startswith("##"):
                reference_metadata_array.append(header)
        # Process reference file and return chrom, pos, ref, and alt
        chrom_pos_ref_alt_array = Parallel(n_jobs=n_jobs)(
            delayed(process_reference_file)(header, line, output_folder_path.joinpath(reference_file_path.name))
            for line in reader
        )

    #######################################################################
    # Write processed reference file
    #######################################################################
    # Read reference file without header then sort chromosomes and positions
    dat = pd.read_table(
        filepath_or_buffer=output_folder_path.joinpath(reference_file_path.name),
        header=None
    )
    dat = dat.sort_values(by=[0, 1])

    # Write metadata and header
    with open(output_folder_path.joinpath(reference_file_path.name), 'w') as writer:
        for text in reference_metadata_array:
            writer.write(text)
        writer.write(header)
    dat.to_csv(
        path_or_buf=output_folder_path.joinpath(reference_file_path.name),
        sep='\t',
        index=False,
        header=False,
        doublequote=False,
        mode='a'
    )

    #######################################################################
    # Make chrom, pos, ref, alt dictionary
    #######################################################################
    chrom_pos_ref_alt_dict = {}
    for i in range(len(chrom_pos_ref_alt_array)):
        if chrom_pos_ref_alt_array[i] is not None:
            if chrom_pos_ref_alt_array[i][0] not in chrom_pos_ref_alt_dict.keys():
                chrom_pos_ref_alt_dict[chrom_pos_ref_alt_array[i][0]] = {
                    chrom_pos_ref_alt_array[i][1]: {
                        'REF': chrom_pos_ref_alt_array[i][2],
                        'ALT': chrom_pos_ref_alt_array[i][3]
                    }
                }
            else:
                if chrom_pos_ref_alt_array[i][1] not in chrom_pos_ref_alt_dict[chrom_pos_ref_alt_array[i][0]].keys():
                    chrom_pos_ref_alt_dict[chrom_pos_ref_alt_array[i][0]][chrom_pos_ref_alt_array[i][1]] = {
                        'REF': chrom_pos_ref_alt_array[i][2],
                        'ALT': chrom_pos_ref_alt_array[i][3]
                    }

    #######################################################################
    # Align input file based on chrom, pos, ref, alt dictionary
    #######################################################################
    input_metadata_array = []
    header = ""
    with open(input_file_path, 'r') as reader:
        # Get metadata and header
        while not header.strip().startswith("#CHROM"):
            header = reader.readline()
            if header.startswith("##"):
                input_metadata_array.append(header)
        # Process input file and return chrom, pos, ref, and alt
        chrom_pos_ref_alt_array = Parallel(n_jobs=n_jobs)(
            delayed(align_input_file_based_on_reference_dictionary)(
                header,
                line,
                chrom_pos_ref_alt_dict,
                output_folder_path.joinpath(input_file_path.name)
            )
            for line in reader
        )

    #######################################################################
    # Remove positions from chrom_pos_ref_alt_dict
    #######################################################################
    for i in range(len(chrom_pos_ref_alt_array)):
        if chrom_pos_ref_alt_array[i] is not None:
            if chrom_pos_ref_alt_array[i][0] in chrom_pos_ref_alt_dict.keys():
                if chrom_pos_ref_alt_array[i][1] in chrom_pos_ref_alt_dict[chrom_pos_ref_alt_array[i][0]].keys():
                    temp = chrom_pos_ref_alt_dict[chrom_pos_ref_alt_array[i][0]].pop(chrom_pos_ref_alt_array[i][1])
    for key in chrom_pos_ref_alt_dict.keys():
        if len(chrom_pos_ref_alt_dict[key].keys()) == 0:
            temp = chrom_pos_ref_alt_dict.pop(key)

    #######################################################################
    # Insert positions into input file
    #######################################################################
    number_of_columns = len(header.strip().split("\t"))
    with open(output_folder_path.joinpath(input_file_path.name), 'a') as writer:
        for chrom_key in chrom_pos_ref_alt_dict.keys():
            for pos_key in chrom_pos_ref_alt_dict[chrom_key].keys():
                line_array = ["./."] * number_of_columns
                line_array[0] = str(chrom_key)
                line_array[1] = str(pos_key)
                line_array[2] = "."
                line_array[3] = chrom_pos_ref_alt_dict[chrom_key][pos_key]['REF']
                line_array[4] = chrom_pos_ref_alt_dict[chrom_key][pos_key]['ALT']
                line_array[5] = "."
                line_array[6] = "PASS"
                line_array[7] = "."
                line_array[8] = "GT"
                writer.write(str("\t".join(line_array))+"\n")

    #######################################################################
    # Write processed input file
    #######################################################################
    # Read input file without header then sort chromosomes and positions
    dat = pd.read_table(
        filepath_or_buffer=output_folder_path.joinpath(input_file_path.name),
        header=None
    )
    dat = dat.sort_values(by=[0, 1])

    # Write metadata and header
    with open(output_folder_path.joinpath(input_file_path.name), 'w') as writer:
        for text in input_metadata_array:
            writer.write(text)
        writer.write(header)
    dat.to_csv(
        path_or_buf=output_folder_path.joinpath(input_file_path.name),
        sep='\t',
        index=False,
        header=False,
        doublequote=False,
        mode='a'
    )


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(
        prog='align_input_VCF_based_on_reference_VCF',
        description='align input VCF based on reference VCF'
    )

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
    parser.add_argument('-r', '--reference_file', help='Reference file', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_folder', help='Output folder', type=pathlib.Path, required=True)

    parser.add_argument('-t', '--threads', help='Number of threads', type=int, default=10)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
