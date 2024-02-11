#!/usr/bin/env python3

# python3 /scratch/yenc/projects/2022_05_25_AlleleCatalog/scripts/generate_Allele_Catalog.py \
# -i /scratch/yenc/projects/2022_05_25_AlleleCatalog/output/Soy1066_Allele_Catalog/Soy1066_Allele_Catalog_Chr01_genotype_ac2.txt \
# -f /scratch/yenc/projects/2022_05_25_AlleleCatalog/output/Soy1066_Allele_Catalog/Soy1066_Allele_Catalog_functional_effect_ac2.txt \
# -m /scratch/yenc/projects/2022_05_25_AlleleCatalog/data/Soy1066_allele_line_info.txt \
# -g /scratch/yenc/datasets/soybean_reference_genome/Wm82.a2.v1.genes.gff \
# -o /scratch/yenc/projects/2022_05_25_AlleleCatalog/output/Soy1066_Allele_Catalog/Soy1066_Allele_Catalog_Chr01_final.txt \
# -c Chr01

import sys
import os
import re
import argparse
import pathlib
import gzip


def generate_allele_catalog(chromosome, gene_id, gene_start, gene_end, input_file_path, functional_effect_header_array, functional_effect_dict, metadata_header_array, metadata_dict, output_array):
	genotype_dict = {}

	with open(input_file_path, 'r') as reader:
		header = reader.readline()
		header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
		for line in reader:
			line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")

			# Check chromosome
			if (str(line_array[0]).strip() == chromosome):
				if (int(gene_start) <= int(str(line_array[1]).strip()) <= int(gene_end)):

					position = str(line_array[1]).strip()
					accession = str(line_array[2]).strip()
					genotype = str(line_array[3]).strip()
					category = str(line_array[4]).strip()
					imputation = str(line_array[6]).strip()

					# Add functional effect
					functional_effect = category
					amino_acid_change = ""
					if (gene_id in functional_effect_dict.keys()):
						if (chromosome in functional_effect_dict[gene_id].keys()):
							if (position in functional_effect_dict[gene_id][chromosome].keys()):
								if (genotype in functional_effect_dict[gene_id][chromosome][position].keys()):
									if functional_effect_dict[gene_id][chromosome][position][genotype]["Functional_Effect"] != "":
										functional_effect = functional_effect_dict[gene_id][chromosome][position][genotype]["Functional_Effect"]
									if functional_effect_dict[gene_id][chromosome][position][genotype]["Amino_Acid_Change"] != "":
										amino_acid_change = functional_effect_dict[gene_id][chromosome][position][genotype]["Amino_Acid_Change"]

					if (accession not in genotype_dict.keys()):
						genotype_dict[accession] = {}

					if (position not in genotype_dict[accession].keys()):
						genotype_dict[accession][position] = {
							"Genotype": genotype,
							"Functional_Effect": functional_effect,
							"Amino_Acid_Change": amino_acid_change,
							"Imputation": imputation
						}
					else:
						genotype_dict[accession][position]["Genotype"] = genotype
						genotype_dict[accession][position]["Functional_Effect"] = functional_effect,
						genotype_dict[accession][position]["Amino_Acid_Change"] = amino_acid_change,
						genotype_dict[accession][position]["Imputation"] = imputation

	output_string = ""
	if genotype_dict:
		for accession in genotype_dict.keys():
			output_string = ""

			# Add metadata
			if (accession in metadata_dict.keys()):
				for metadata_key in metadata_header_array[1:]:
					if (metadata_key in metadata_dict[accession].keys()):
						output_string = output_string + metadata_dict[accession][metadata_key] + '\t'
					else:
						output_string = output_string + '' + '\t'
			else:
				output_string = output_string + str('\t'.join(['']*len(metadata_header_array[1:]))) + '\t'

			# Add accession, gene, and chromosome
			output_string = output_string + accession + '\t' + gene_id + '\t' + chromosome + '\t'

			if (accession in genotype_dict.keys()):

				# Add positions
				position_array = list(genotype_dict[accession].keys())
				position_array.sort()

				output_string = output_string + str(' '.join(position_array)) + '\t'

				# Add genotypes
				output_string = output_string + str(' '.join([genotype_dict[accession][position]['Genotype'] for position in position_array])) + '\t'

				# Add genotypes with description
				for position in position_array:
					output_string = output_string + str(
						genotype_dict[accession][position]['Genotype'] + '|' +
						genotype_dict[accession][position]['Functional_Effect'] + '|' +
						genotype_dict[accession][position]['Amino_Acid_Change'] + '|' +
						genotype_dict[accession][position]['Imputation']
					).strip('|') + ' '
				output_string = output_string + '\t'

				# Add imputation
				if '+' in [genotype_dict[accession][position]['Imputation'] for position in position_array]:
					output_string = output_string + '+' + '\n'
				else:
					output_string = output_string + '-' + '\n'

				# Append output_string to output_array
				output_array.append(output_string)


def main(args):
	#######################################################################
	# Get arguments
	#######################################################################
	input_file_path = args.input_file
	functional_effect_file_path = args.functional_effect_file
	metadata_file_path = args.metadata_file
	gff_file_path = args.gff_file
	output_file_path = args.output_file

	chromosome_array = args.chromosome
	gff_category = args.gff_category
	gff_key = args.gff_key

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
	# Load functional effect file
	#######################################################################
	functional_effect_dict = {}

	with open(functional_effect_file_path, 'r') as reader:
		functional_effect_header = reader.readline()
		functional_effect_header_array = str(functional_effect_header).strip("\n").strip("\r").strip("\r\n").split("\t")
		for line in reader:
			line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")

			chromosome = str(line_array[0]).strip()
			position = str(line_array[1]).strip()
			allele = str(line_array[2]).strip()
			functional_effect = str(line_array[3]).strip()
			gene_name = str(line_array[4]).strip()
			amino_acid_change = str(line_array[8]).strip()

			if (gene_name not in functional_effect_dict.keys()):
				functional_effect_dict[gene_name] = {}

			if (chromosome not in functional_effect_dict[gene_name].keys()):
				functional_effect_dict[gene_name][chromosome] = {}

			if (position not in functional_effect_dict[gene_name][chromosome].keys()):
				functional_effect_dict[gene_name][chromosome][position] = {}

			if (allele not in functional_effect_dict[gene_name][chromosome][position].keys()):
				functional_effect_dict[gene_name][chromosome][position][allele] = {
					"Functional_Effect": functional_effect,
					"Amino_Acid_Change": amino_acid_change
				}
			else:
				if functional_effect_dict[gene_name][chromosome][position][allele]["Functional_Effect"] == "":
					functional_effect_dict[gene_name][chromosome][position][allele]["Functional_Effect"] = functional_effect
				if functional_effect_dict[gene_name][chromosome][position][allele]["Amino_Acid_Change"] == "":
					functional_effect_dict[gene_name][chromosome][position][allele]["Amino_Acid_Change"] = amino_acid_change

	#######################################################################
	# Load metadata file
	#######################################################################
	metadata_dict = {}

	with open(metadata_file_path, 'r') as reader:
		metadata_header = reader.readline()
		metadata_header_array = str(metadata_header).strip("\n").strip("\r").strip("\r\n").split("\t")
		for line in reader:
			line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")

			accession = str(line_array[0]).strip()

			if (accession not in metadata_dict.keys()):
				metadata_dict[accession] = {}

			for i in range(1, len(metadata_header_array)):
				metadata_dict[accession][str(metadata_header_array[i]).strip()] = str(line_array[i]).strip()

	#######################################################################
	# Load GFF file
	#######################################################################
	gff_dict = {}

	with open(gff_file_path, 'r') as reader:
		for line in reader:
			if not line.startswith('#'):
				line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")

				if (str(str(line_array[2]).strip()) == gff_category):

					chromosome = str(line_array[0]).strip()
					start = str(line_array[3]).strip()
					end = str(line_array[4]).strip()
					attributes = str(line_array[8]).strip()

					gene_id = re.sub('(;.*)', '', re.sub(re.compile(".*"+str(gff_key)), '', attributes))

					if (chromosome_array is None):
						if (chromosome not in gff_dict.keys()):
							gff_dict[chromosome] = {}

						if (gene_id not in gff_dict[chromosome].keys()):
							gff_dict[chromosome][gene_id] = {
								"Start": start,
								"End": end
							}
						else:
							gff_dict[chromosome][gene_id]['Start'] = start
							gff_dict[chromosome][gene_id]['End'] = end
					else:
						if (chromosome in chromosome_array):
							if (chromosome not in gff_dict.keys()):
								gff_dict[chromosome] = {}

							if (gene_id not in gff_dict[chromosome].keys()):
								gff_dict[chromosome][gene_id] = {
									"Start": start,
									"End": end
								}
							else:
								gff_dict[chromosome][gene_id]['Start'] = start
								gff_dict[chromosome][gene_id]['End'] = end

	#######################################################################
	# Write output file header
	#######################################################################
	with open(output_file_path, 'w') as writer:
		writer.write(str('\t'.join(metadata_header_array[1:])) + "\tAccession\tGene\tChromosome\tPosition\tGenotype\tGenotype_with_Description\tImputation\n")

	#######################################################################
	# Generate Allele Catalog
	#######################################################################
	chunksize = 10000
	output_array = []

	for chromosome in gff_dict.keys():
		for gene_id in gff_dict[chromosome].keys():
			generate_allele_catalog(
				chromosome,
				gene_id,
				gff_dict[chromosome][gene_id]['Start'],
				gff_dict[chromosome][gene_id]['End'],
				input_file_path,
				functional_effect_header_array,
				functional_effect_dict,
				metadata_header_array,
				metadata_dict,
				output_array
			)
			# Check and write data
			if (len(output_array) > chunksize):
				with open(output_file_path, 'a') as writer:
					writer.write(''.join(output_array))
					output_array.clear()

	# Check and write data
	if (len(output_array) > 0):
		with open(output_file_path, 'a') as writer:
			writer.write(''.join(output_array))
			output_array.clear()


if __name__ == "__main__":
	#######################################################################
	# Parse arguments
	#######################################################################
	parser = argparse.ArgumentParser(prog='generate_Allele_Catalog', description='generate_Allele_Catalog')

	parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
	parser.add_argument('-f', '--functional_effect_file', help='Functional effect file', type=pathlib.Path, required=True)
	parser.add_argument('-m', '--metadata_file', help='Metadata file', type=pathlib.Path, required=True)
	parser.add_argument('-g', '--gff_file', help='GFF file', type=pathlib.Path, required=True)
	parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

	parser.add_argument('-c', '--chromosome', help='Chromosome', type=str, action='append')
	parser.add_argument('-a', '--gff_category', help='Gff category', type=str, default='gene')
	parser.add_argument('-k', '--gff_key', help='Gff key', type=str, default='Name')

	args = parser.parse_args()

	#######################################################################
	# Call main function
	#######################################################################
	main(args)
