#!/usr/bin/env python3

import sys
import os
import re
import math
import argparse
import pathlib
import threading

from joblib import Parallel, delayed

import pandas as pd


# Generate gene array
# This function is for joblib parallel
def generate_gene_array(line, gene_array):
    values = line.strip("\n").split("\t")
    gene = values[7]
    if gene not in gene_array:
        gene_array.append(gene)


# Generate wide allele catalog
# This function is for joblib parallel
def generate_allele_catalog_wide(input_file_path, genes, output_file_path, lock):
    # Read data of specific genes
    df = pd.DataFrame(columns = ['Classification', 'Improvement_Status', 'Maturity_Group', 'Country', 'State', 'Accession', 'Chromosome', 'Gene', 'Position', 'Genotype', 'Genotype_with_Description'])
    chunksize = 100000
    for dat in pd.read_table(filepath_or_buffer=input_file_path, chunksize=chunksize):
        dat = dat[dat['Gene'].isin(genes)]
        if dat.shape[0] > 0:
            df = pd.concat([df, dat])

    # Process data if table has at least one row
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

        # Save table to output file
        lock.acquire()
        if (df is not None) and (df.shape[0] > 0) and (not output_file_path.exists()):
            df.to_csv(path_or_buf=output_file_path, sep='\t', index=False, header=True, doublequote=False, mode='w')
        elif (df is not None) and (df.shape[0] > 0) and (output_file_path.exists()) and (output_file_path.is_file()):
            df.to_csv(path_or_buf=output_file_path, sep='\t', index=False, header=False, doublequote=False, mode='a')
        lock.release()


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_path = args.input_file
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
    if output_file_path.exists():
        output_file_path.unlink()

    #######################################################################
    # Collect genes in the input file
    #######################################################################
    grouped_gene_array = []
    gene_array = []
    n = 10
    with open(input_file_path, 'r') as reader:
        header = reader.readline()
        Parallel(n_jobs=n_jobs, backend="threading")(
            delayed(generate_gene_array)(line, gene_array)
            for line in reader
        )
    if len(gene_array) > 0:
        for i in range(math.ceil(len(gene_array)/n)):
            i_start = i*n
            i_stop = i*n+n if i*n+n < len(gene_array) else len(gene_array)
            grouped_gene_array.append(gene_array[i_start:i_stop])
    else:
        df = pd.DataFrame(
            columns=['Classification', 'Improvement_Status', 'Maturity_Group', 'Country', 'State', 'Accession',
                     'Chromosome', 'Gene', 'Position', 'Genotype', 'Genotype_with_Description'])
        df.to_csv(path_or_buf=output_file_path, sep='\t', index=False, header=True, doublequote=False, mode='w')
        sys.exit(0)

    #######################################################################
    # Generate allele catalog wide
    #######################################################################
    Parallel(n_jobs=n_jobs, backend="threading")(
        delayed(generate_allele_catalog_wide)(input_file_path, genes, output_file_path, lock)
        for genes in grouped_gene_array 
    )
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

    parser.add_argument('-t', '--threads', help='Number of threads', type=int, default=10)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)