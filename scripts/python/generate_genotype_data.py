#!/usr/bin/env python3

# python3 /scratch/yenc/projects/2022_05_25_AlleleCatalog/scripts/generate_genotype_data.py \
# -i /scratch/yenc/projects/2022_05_25_AlleleCatalog/data/Bcftools_view_reheader/Soy1066_Chr01.vcf.gz \
# -m /scratch/yenc/projects/2022_05_25_AlleleCatalog/output/Soy1066_Allele_Catalog/Soy1066_Allele_Catalog_imputation_ac2.txt \
# -o /scratch/yenc/projects/2022_05_25_AlleleCatalog/output/Soy1066_Allele_Catalog/Soy1066_Allele_Catalog_Chr01_genotype_ac2.txt

import sys
import os
import re
import argparse
import pathlib
import gzip


# Process line in VCF file
def process_line(header_array, line_array, imputation_dict, output_array):
	chromosome = line_array[0]
	position = line_array[1]
	reference_allele = line_array[3]
	alternate_alleles = line_array[4]

	alleles = re.split(',', reference_allele) + re.split(',', alternate_alleles)

	accession_array = header_array[9:]

	genotype_index_array = [int(str(re.split('/|\\|', genotype_index)[0]).strip()) for genotype_index in line_array[9:]]

	genotype_array = [alleles[genotype_index] for genotype_index in genotype_index_array]

	category_array = ['Ref' if genotype_index == 0 else 'Alt' for genotype_index in genotype_index_array]

	category_index_array = [0 if genotype_index == 0 else 1 for genotype_index in genotype_index_array]


	for i in range(len(accession_array)):

		accession = str(accession_array[i]).strip()
		genotype = str(genotype_array[i]).strip()
		category = str(category_array[i]).strip()
		category_index = str(category_index_array[i]).strip()
		imputation = ''

		if (chromosome in imputation_dict.keys()):
			if (position in imputation_dict[chromosome].keys()):
				if (accession in imputation_dict[chromosome][position]):
					imputation = '+'

		output_array.append(
			chromosome + "\t" +
			position + "\t" +
			accession + "\t" +
			genotype + "\t" +
			category + "\t" +
			category_index + "\t" +
			imputation + "\n"
		)


def main(args):
	#######################################################################
	# Get arguments
	#######################################################################
	input_file_path = args.input_file
	imputation_file_path = args.imputation_file
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
		writer.write("Chromosome\tPosition\tAccession\tGenotype\tCategory\tCategory_Index\tImputation\n")

	#######################################################################
	# Load imputation file
	#######################################################################
	imputation_dict = {}

	with open(imputation_file_path, 'r') as reader:
		header = reader.readline()
		header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
		for line in reader:
			line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")

			chromosome = str(line_array[0]).strip()
			position = str(line_array[1]).strip()
			accession = str(line_array[2]).strip()

			if (chromosome not in imputation_dict.keys()):
				imputation_dict[chromosome] = {}

			if (position not in imputation_dict[chromosome].keys()):
				imputation_dict[chromosome][position] = []

			if (accession not in imputation_dict[chromosome][position]):
				imputation_dict[chromosome][position].append(accession)

	#######################################################################
	# Process VCF files
	#######################################################################
	chunksize = 10000
	output_array = []

	if str(input_file_path).endswith('gz'):
		with gzip.open(input_file_path, 'rt') as reader:
			header = ""
			while not header.strip().startswith("#CHROM"):
				header = reader.readline()
				header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
			for line in reader:
				line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
				process_line(header_array, line_array, imputation_dict, output_array)
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
				process_line(header_array, line_array, imputation_dict, output_array)
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
	parser = argparse.ArgumentParser(prog='generate_genotype_data', description='generate_genotype_data')

	parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
	parser.add_argument('-m', '--imputation_file', help='Imputation file', type=pathlib.Path, required=True)
	parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

	args = parser.parse_args()

	#######################################################################
	# Call main function
	#######################################################################
	main(args)
