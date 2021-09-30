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
    if output_file_path.exists():
        output_file_path.unlink()

    #######################################################################
    # Collect genes in the input file
    #######################################################################
    accessions = []
    genes = []
    chunksize = 100000
    for dat in pd.read_table(filepath_or_buffer=input_file_path, usecols=['Accession', 'Gene'], chunksize=chunksize):
        for accession in dat['Accession']:
            if accession not in accessions:
                accessions.append(accession)
        for gene in dat['Gene']:
            if gene not in genes:
                genes.append(gene)
    accessions.sort()
    genes.sort()

    #######################################################################
    # Generate allele catalog wide
    #######################################################################
    for gene in genes:
        # Create a data frame that has header only
        df = pd.DataFrame(
            columns=[
                'Classification', 'Improvement_Status', 'Maturity_Group',
                'Country', 'State', 'Accession', 'Chromosome', 'Gene',
                'Position', 'Genotype', 'Genotype_with_Description'
            ]
        )
        # Read data that has that specific gene
        chunksize = 100000
        for dat in pd.read_table(filepath_or_buffer=input_file_path, chunksize=chunksize):
            dat = dat[dat['Gene'] == gene]
            if dat.shape[0] > 0:
                df = pd.concat([df, dat], ignore_index=True)
        # Restructure the data frame
        if df.shape[0] > 0:
            df = df.drop_duplicates(subset=['Accession', 'Chromosome', 'Gene', 'Position'])
            df = df.sort_values(by=['Accession', 'Chromosome', 'Gene', 'Position'])
            column_names = df.columns.tolist()
            subset_of_column_names = df.columns.tolist()
            subset_of_column_names.remove('Position')
            subset_of_column_names.remove('Genotype')
            subset_of_column_names.remove('Genotype_with_Description')
            df = df.groupby(subset_of_column_names, dropna=False).aggregate({
                'Position': lambda x: ' '.join(map(str, x)),
                'Genotype': lambda x: ' '.join(map(str, x)),
                'Genotype_with_Description': lambda x: ' '.join(map(str, x))
            }).reset_index()
            df = df.drop_duplicates(subset=['Accession', 'Chromosome', 'Gene', 'Position'])
            df = df.loc[:, column_names]
            df = df.sort_values(by=['Accession', 'Chromosome', 'Gene'])
            if (df is not None) and (df.shape[0] > 0) and (not output_file_path.exists()):
                df.to_csv(path_or_buf=output_file_path, sep='\t', index=False, header=True, doublequote=False, mode='w')
            elif (df is not None) and (df.shape[0] > 0) and (output_file_path.exists()) and (
            output_file_path.is_file()):
                df.to_csv(path_or_buf=output_file_path, sep='\t', index=False, header=False, doublequote=False,
                          mode='a')
    # If the output file still does not exists, write an empty table with header only
    if not output_file_path.exists():
        df = pd.DataFrame(
            columns=['Classification', 'Improvement_Status', 'Maturity_Group', 'Country', 'State', 'Accession',
                     'Chromosome', 'Gene', 'Position', 'Genotype', 'Genotype_with_Description'])
        df.to_csv(path_or_buf=output_file_path, sep='\t', index=False, header=True, doublequote=False, mode='w')


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
