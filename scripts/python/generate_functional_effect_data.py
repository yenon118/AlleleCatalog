#!/usr/bin/env python3

# python3 /scratch/yenc/projects/2022_05_25_AlleleCatalog/scripts/generate_functional_effect_data.py \
# -i /scratch/yenc/projects/2022_05_25_AlleleCatalog/data/Bcftools_view_reheader/Soy1066_Chr01.vcf.gz \
# -i /scratch/yenc/projects/2022_05_25_AlleleCatalog/data/Bcftools_view_reheader/Soy1066_Chr02.vcf.gz \
# -i /scratch/yenc/projects/2022_05_25_AlleleCatalog/data/Bcftools_view_reheader/Soy1066_Chr03.vcf.gz \
# -o /scratch/yenc/projects/2022_05_25_AlleleCatalog/output/Soy1066_Allele_Catalog/Soy1066_Allele_Catalog_functional_effect_ac2.txt

import sys
import os
import re
import argparse
import pathlib
import gzip

import pandas as pd


# These are the effects we need to focus on.
EFFECTS = [
    'frameshift_variant',
    'exon_loss_variant',
    'duplication',
    'inversion',
    'feature_ablation',
    'gene_fusion',
    'rearranged_at_DNA_level',
    'missense_variant',
    'protein_protein_contact',
    'structural_interaction_variant',
    'rare_amino_acid_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_lost',
    'start_lost',
    'stop_gained',
    'inframe_insertion',
    'disruptive_inframe_insertion',
    'inframe_deletion',
    'disruptive_inframe_deletion'
]


# Clean amino acid string
def clean_amino_acid_change_string(amino_acid_change):
    amino_acid_change = re.sub('Gly', 'G', amino_acid_change)
    amino_acid_change = re.sub('Pro', 'P', amino_acid_change)
    amino_acid_change = re.sub('Val', 'V', amino_acid_change)
    amino_acid_change = re.sub('Leu', 'L', amino_acid_change)
    amino_acid_change = re.sub('Met', 'M', amino_acid_change)
    amino_acid_change = re.sub('Cys', 'C', amino_acid_change)
    amino_acid_change = re.sub('Phe', 'F', amino_acid_change)
    amino_acid_change = re.sub('Tyr', 'Y', amino_acid_change)
    amino_acid_change = re.sub('Trp', 'W', amino_acid_change)
    amino_acid_change = re.sub('His', 'H', amino_acid_change)
    amino_acid_change = re.sub('Lys', 'K', amino_acid_change)
    amino_acid_change = re.sub('Arg', 'R', amino_acid_change)
    amino_acid_change = re.sub('Gln', 'Q', amino_acid_change)
    amino_acid_change = re.sub('Asp', 'D', amino_acid_change)
    amino_acid_change = re.sub('Ser', 'S', amino_acid_change)
    amino_acid_change = re.sub('Thr', 'T', amino_acid_change)
    amino_acid_change = re.sub('Asn', 'N', amino_acid_change)
    amino_acid_change = re.sub('Ile', 'I', amino_acid_change)
    amino_acid_change = re.sub('Glu', 'E', amino_acid_change)
    amino_acid_change = re.sub('Ala', 'A', amino_acid_change)
    amino_acid_change = re.sub('p\\.', '', amino_acid_change)
    return amino_acid_change


# Process line in VCF file
def process_line(header_array, line_array, output_array):
    chromosome = line_array[0]
    position = line_array[1]

    annotation_string = re.sub('(.*ANN=)|(;.*)', '', line_array[7])

    annotation_array = [re.split('\\|', annotation) for annotation in re.split(';', annotation_string)]

    for i in range(len(annotation_array)):
        functional_effect = str(annotation_array[i][1]).strip()
        amino_acid_change = str(annotation_array[i][10]).strip()

        if any([str(functional_effect).find(effect) != -1 for effect in EFFECTS]):
            output_array.append(
                chromosome + "\t" +
                position + "\t" +
                str(annotation_array[i][0]).strip() + "\t" +
                functional_effect + "\t" +
                str(annotation_array[i][3]).strip() + "\t" +
                str(annotation_array[i][4]).strip() + "\t" +
                str(annotation_array[i][5]).strip() + "\t" +
                str(annotation_array[i][6]).strip() + "\t" +
                clean_amino_acid_change_string(amino_acid_change) + "\n"
            )


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
    # Write output file header
    #######################################################################
    with open(output_file_path, 'w') as writer:
        writer.write("Chromosome\tPosition\tAllele\tFunctional_Effect\tGene_Name\tGene\tFeature_Type\tFeature\tAmino_Acid_Change\n")

    #######################################################################
    # Process VCF files
    #######################################################################
    chunksize = 10000
    output_array = []

    for input_file_path in input_file_paths:
        if str(input_file_path).endswith('gz'):
            with gzip.open(input_file_path, 'rt') as reader:
                header = ""
                while not header.strip().startswith("#CHROM"):
                    header = reader.readline()
                    header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
                for line in reader:
                    line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
                    process_line(header_array, line_array, output_array)
                    # Check and write data
                    if (len(output_array) > chunksize):
                        with open(output_file_path, 'a') as writer:
                            writer.write("".join(output_array))
                            output_array.clear()
        else:
            with open(input_file_path, "r") as reader:
                header = ""
                while not header.strip().startswith("#CHROM"):
                    header = reader.readline()
                    header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
                for line in reader:
                    line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
                    process_line(header_array, line_array, output_array)
                    # Check and write data
                    if (len(output_array) > chunksize):
                        with open(output_file_path, 'a') as writer:
                            writer.write("".join(output_array))
                            output_array.clear()

    # Write the remaining data
    if (len(output_array) > 0):
        with open(output_file_path, 'a') as writer:
            writer.write("".join(output_array))
            output_array.clear()

    #######################################################################
    # Process functional effects
    #######################################################################
    dat = pd.read_table(
        filepath_or_buffer=output_file_path
    )

    dat = dat.sort_values(by=['Feature'])

    dat = dat.drop_duplicates(subset=['Chromosome', 'Position', 'Allele'])

    dat = dat.sort_values(by=['Chromosome', 'Position', 'Allele'])

    dat.to_csv(
        path_or_buf=output_file_path,
        sep='\t',
        index=False,
        doublequote=False,
        mode='w'
    )


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='generate_functional_effect_data', description='generate_functional_effect_data')

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, action='append', required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
