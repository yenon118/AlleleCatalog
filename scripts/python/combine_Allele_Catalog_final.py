#!/usr/bin/env python3

import sys
import os
import re
import argparse
import pathlib

import pandas as pd


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_paths = args.input_file
    output_file_path = args.output_file

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
    # Merge all Allele Catalog together
    #######################################################################
    header = True
    for i in range(len(input_file_paths)):
        dat = pd.read_table(
            filepath_or_buffer=input_file_paths[i],
            dtype={
                "Classification": str,
                "Improvement_Status": str,
                "Maturity_Group": str,
                "Country": str,
                "State": str,
                "Accession": str,
                "Chromosome": str,
                "Gene": str,
                "Position": str,
                "Genotype": str,
                "Genotype_with_Description": str,
                "Ancestry_Binary": str
            }
        )
        dat = dat.sort_values(by=['Accession', 'Chromosome', 'Gene'])
        if header:
            dat.to_csv(
                path_or_buf=output_file_path, sep='\t', index=False, header=header, doublequote=False, mode='w'
            )
            header = False
        else:
            dat.to_csv(
                path_or_buf=output_file_path, sep='\t', index=False, header=header, doublequote=False, mode='a'
            )


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='combine_Allele_Catalog_wide', description='combine Allele Catalog (wide)')

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, action='append', required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
