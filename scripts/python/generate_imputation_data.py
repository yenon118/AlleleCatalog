#!/usr/bin/env python3

# python3 /scratch/yenc/projects/2022_05_25_AlleleCatalog/scripts/generate_imputation_data.py \
# -i /scratch/yenc/projects/2022_05_25_AlleleCatalog/data/reheader_unimpute_data/Soy1066_Chr01.vcf.gz \
# -i /scratch/yenc/projects/2022_05_25_AlleleCatalog/data/reheader_unimpute_data/Soy1066_Chr02.vcf.gz \
# -i /scratch/yenc/projects/2022_05_25_AlleleCatalog/data/reheader_unimpute_data/Soy1066_Chr03.vcf.gz \
# -o /scratch/yenc/projects/2022_05_25_AlleleCatalog/output/Soy1066_Allele_Catalog/Soy1066_Allele_Catalog_imputation_ac2.txt

import sys
import os
import re
import argparse
import pathlib
import gzip

import pandas as pd


# Process line in VCF file
def process_line(header_array, line_array, output_array):
	chromosome = line_array[0]
	position = line_array[1]

	genotype_index_array = [str(re.split('/|\\|', genotype_index)[0]).strip() for genotype_index in line_array[9:]]

	missing_value_indexes = [i+9 for i in range(len(genotype_index_array)) if genotype_index_array[i] == '.']

	imputed_accessions = [header_array[missing_value_index] for missing_value_index in missing_value_indexes]

	for i in range(len(imputed_accessions)):
		output_array.append(
			chromosome + "\t" +
			position + "\t" +
			imputed_accessions[i] + "\t" +
			"+" + "\n"
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
		writer.write("Chromosome\tPosition\tAccession\tImputation\n")

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


if __name__ == "__main__":
	#######################################################################
	# Parse arguments
	#######################################################################
	parser = argparse.ArgumentParser(prog='generate_imputation_data', description='generate_imputation_data')

	parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, action='append', required=True)
	parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

	args = parser.parse_args()

	#######################################################################
	# Call main function
	#######################################################################
	main(args)
