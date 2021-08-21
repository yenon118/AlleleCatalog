#!/usr/bin/env python3

import sys
import os
import re
import argparse
import pathlib

from collections import OrderedDict
from joblib import Parallel, delayed


def append_ancestry_binaries(header, line, ancestor_dict_arr, output_file_path):
    line = str(line).strip("\n")
    line_array = line.split("\t")
    ancestry_binaries_arr = []
    # Check genotypes with ancestors' genotypes
    if line_array[7] in ancestor_dict_arr.keys():
        for key in ancestor_dict_arr[line_array[7]].keys():
            if line_array[9] == ancestor_dict_arr[line_array[7]][key]:
                ancestry_binaries_arr.append(1)
            else:
                ancestry_binaries_arr.append(0)
    # Make ancestry binaries a string
    if len(ancestry_binaries_arr) == 0 or sum(ancestry_binaries_arr) == 0:
        ancestry_binaries = ""
    else:
        ancestry_binaries = "".join([str(e) for e in ancestry_binaries_arr])
    # Append to file
    with open(output_file_path, "a") as writer:
        writer.write(line + "\t" + ancestry_binaries + "\n")


def collect_all_ancestor_entries(header, line, ancestors):
    line = str(line).strip()
    line_array = line.split("\t")
    if ancestors is not None:
        if len(ancestors) != 0:
            if line_array[5] in ancestors:
                return {
                    (line_array[5], line_array[7]): line_array[9]
                }
            else:
                return None
        else:
            return None
    else:
        return None


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_path = args.input_file
    output_file_path = args.output_file

    n_jobs = args.threads
    ancestors = args.ancestor

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
    # Read input file then collect all ancestor entries
    #######################################################################
    with open(input_file_path, "r") as reader:
        header = ""
        while not header.strip().startswith("Classification"):
            header = reader.readline()
        temp_ancestor_dict_arr = Parallel(n_jobs=n_jobs)(
            delayed(collect_all_ancestor_entries)(header, line, ancestors)
            for line in reader
        )
    ancestor_dict_arr = OrderedDict()
    if temp_ancestor_dict_arr is not None:
        if len(temp_ancestor_dict_arr) > 0:
            # Parse to get ancestry dictionary
            for i in range(len(temp_ancestor_dict_arr)):
                if temp_ancestor_dict_arr[i] is not None:
                    for key in temp_ancestor_dict_arr[i].keys():
                        if key[1] not in ancestor_dict_arr.keys():
                            ancestor_dict_arr[key[1]] = OrderedDict()
                            ancestor_dict_arr[key[1]][key[0]] = temp_ancestor_dict_arr[i][key]
                        else:
                            ancestor_dict_arr[key[1]][key[0]] = temp_ancestor_dict_arr[i][key]

    #######################################################################
    # Read input file then append ancestry binaries
    #######################################################################
    with open(input_file_path, "r") as reader:
        header = ""
        while not header.strip().startswith("Classification"):
            header = reader.readline()
        if header != "":
            with open(output_file_path, "w") as writer:
                writer.write(header.strip() + "\t" + "Ancestry_Binary" + "\n")
        Parallel(n_jobs=n_jobs, backend="threading")(
            delayed(append_ancestry_binaries)(header, line, ancestor_dict_arr, output_file_path)
            for line in reader
        )


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='combine_Allele_Catalog_wide', description='combine Allele Catalog (wide)')

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    parser.add_argument('-t', '--threads', help='Number of threads', type=int, default=10)
    parser.add_argument('-a', '--ancestor', help='Ancestor Accession', type=str, action='append', default=None)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
