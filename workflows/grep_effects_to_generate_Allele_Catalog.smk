import sys
import os
import re

project_name = config['project_name']
workflow_path = config['workflow_path']
input_files = config['input_files']
chromosomes = config['chromosomes']
unimputed_input_files = config['unimputed_input_files']
metadata_file = config['metadata_file']
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


unimputed_input_folder = ''
unimputed_input_sample = ''
unimputed_input_samples = []
unimputed_input_extension = ''


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


for i in range(len(unimputed_input_files)):
	if os.path.dirname(unimputed_input_files[i]) != unimputed_input_folder:
		unimputed_input_folder = os.path.dirname(unimputed_input_files[i])
	possible_sample = re.sub('(\\.vcf.*)', '', str(os.path.basename(unimputed_input_files[i])))
	if possible_sample not in unimputed_input_samples:
		unimputed_input_samples.append(possible_sample)
	possible_extension = re.sub(possible_sample,'',str(os.path.basename(unimputed_input_files[i])))
	if possible_extension != unimputed_input_extension:
		unimputed_input_extension = possible_extension
	possible_sample = re.sub(re.compile('_'+str(chromosomes[i])+'$'), '', possible_sample)
	if unimputed_input_sample != possible_sample:
		unimputed_input_sample = possible_sample


rule all:
	input:
		expand(os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'generate_imputation_data', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'generate_functional_effect_data', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'generate_genotype_data', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog', input_sample+'_{chromosome}.txt'), chromosome=chromosomes)


rule grep_effects:
	input:
		in_file = os.path.join(os.path.abspath(input_folder), input_sample+'_{chromosome}'+input_extension)
	output:
		out_file = os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension)
	priority: 70
	shell:
		"""
		zgrep -e "^#" -e "^#CHROM" -e "frameshift_variant" -e "exon_loss_variant" -e "duplication" -e "inversion" \
		-e "feature_ablation" -e "gene_fusion" -e "rearranged_at_DNA_level" -e "missense_variant" -e "protein_protein_contact" \
		-e "structural_interaction_variant" -e "rare_amino_acid_variant" -e "splice_acceptor_variant" -e "splice_donor_variant" \
		-e "stop_lost" -e "start_lost" -e "stop_gained" -e "inframe_insertion" -e "disruptive_inframe_insertion" -e "inframe_deletion" \
		-e "disruptive_inframe_deletion" {input.in_file} | grep -v -e "DEL" -e "INS" | bgzip > {output.out_file}
		"""


rule generate_imputation_data:
	input:
		in_file = os.path.join(os.path.abspath(unimputed_input_folder), unimputed_input_sample+'_{chromosome}'+unimputed_input_extension)
	output:
		out_file = os.path.join(os.path.abspath(output_folder), 'generate_imputation_data', input_sample+'_{chromosome}.txt')
	log:
		os.path.join(os.path.abspath(output_folder), 'generate_imputation_data_log', input_sample+'_{chromosome}.log')
	threads: threads
	shell:
		"""
		python3 {workflow_path}/scripts/python/generate_imputation_data.py -i {input.in_file} -o {output.out_file} 2>&1 > {log}
		"""


include: './../rules/python/generate_functional_effect_data.smk'
include: './../rules/python/generate_genotype_data.smk'
include: './../rules/python/generate_Allele_Catalog.smk'
