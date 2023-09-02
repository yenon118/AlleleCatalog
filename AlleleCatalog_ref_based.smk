import sys
import os
import re


project_name = config['project_name']
workflow_path = config['workflow_path']
input_files = config['input_files']
reference_files = config['reference_files']
chromosomes = config['chromosomes']
beagle_window = config['beagle_window']
metadata_file = config['metadata_file']
genome_version = config['genome_version']
gff_file = config['gff_file']
gff_category = config['gff_category']
gff_key = config['gff_key']
output_folder = config['output_folder']
memory = config['memory']
threads = config['threads']


input_folder = ''
input_sample = ''
input_samples = []
input_extension = ''


reference_folder = ''
reference_sample = ''
reference_samples = []
reference_extension = ''


for i in range(len(input_files)):
	if os.path.dirname(input_files[i]) != input_folder:
		input_folder = os.path.dirname(input_files[i])
	possible_sample = re.sub('(\\.vcf.*)', '', str(os.path.basename(input_files[i])))
	if possible_sample not in input_samples:
		input_samples.append(possible_sample)
	possible_extension = re.sub(possible_sample,'',str(os.path.basename(input_files[i])))
	if possible_extension != input_extension:
		input_extension = possible_extension
	possible_sample = re.sub(re.compile('_'+str(chromosomes[i])+'$'), '', possible_sample)
	if input_sample != possible_sample:
		input_sample = possible_sample


for i in range(len(reference_files)):
	if os.path.dirname(reference_files[i]) != reference_folder:
		reference_folder = os.path.dirname(reference_files[i])
	possible_sample = re.sub('(\\.vcf.*)', '', str(os.path.basename(reference_files[i])))
	if possible_sample not in reference_samples:
		reference_samples.append(possible_sample)
	possible_extension = re.sub(possible_sample,'',str(os.path.basename(reference_files[i])))
	if possible_extension != reference_extension:
		reference_extension = possible_extension
	possible_sample = re.sub(re.compile('_'+str(chromosomes[i])+'$'), '', possible_sample)
	if reference_sample != possible_sample:
		reference_sample = possible_sample


rule all:
	input:
		expand(os.path.join(os.path.abspath(output_folder), 'beagle_impute_input_file', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'snpEff_input_file', input_sample+'_{chromosome}.html'), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'snpEff_input_file', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'generate_functional_effect_data', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'generate_imputation_data', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'generate_genotype_data', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog', input_sample+'_{chromosome}.txt'), chromosome=chromosomes)


include: './rules/java/beagle_impute_input_file_with_reference_file.smk'
include: './rules/java/snpEff_input_file.smk'

include: './rules/bash/grep_input_effects.smk'

include: './rules/python/generate_functional_effect_data.smk'
include: './rules/python/generate_imputation_data.smk'
include: './rules/python/generate_genotype_data.smk'
include: './rules/python/generate_Allele_Catalog.smk'
