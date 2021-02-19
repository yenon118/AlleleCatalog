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
    input_file_path = args.input_file
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

    dat = pd.read_table(filepath_or_buffer=input_file_path)

    if dat.shape[0] > 0:
        dat = dat.drop_duplicates(subset=['Accession', 'Chromosome', 'Gene', 'Position'])

        dat = dat.sort_values(by=['Accession', 'Chromosome', 'Gene', 'Position'])

        column_names = dat.columns.tolist()

        subset_of_column_names = dat.columns.tolist()
        subset_of_column_names.remove('Position')
        subset_of_column_names.remove('Genotype')
        subset_of_column_names.remove('Genotype_with_Description')

        dat = dat.groupby(subset_of_column_names, dropna=False).aggregate({
            'Position': lambda x: ' '.join(map(str, x)),
            'Genotype': lambda x: ' '.join(map(str, x)),
            'Genotype_with_Description': lambda x: ' '.join(map(str, x))
        }).reset_index()

        dat = dat.drop_duplicates(subset=['Accession', 'Chromosome', 'Gene', 'Position'])

        dat = dat.loc[:, column_names]

        dat = dat.sort_values(by=['Accession', 'Chromosome', 'Gene'])

    dat.to_csv(path_or_buf=output_file_path, sep='\t', index=False, header=True, doublequote=False)


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='generate_Allele_Catalog_wide', description='generate Allele Catalog (wide)')

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
