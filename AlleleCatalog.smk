import sys
import os
import re


project_name = config['project_name']
workflow_path = config['workflow_path']
input_files = config['input_files']
chromosomes = config['chromosomes']
reference_file = config['reference_file']
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


for i in range(len(input_files)):
    if os.path.dirname(input_files[i]) != input_folder:
        input_folder = os.path.dirname(input_files[i])
    possible_sample = re.sub('(\\.vcf.*)', '', str(os.path.basename(input_files[i])))
    if possible_sample not in input_samples:
        input_samples.append(possible_sample)
    possible_extension = re.sub(possible_sample,'',str(os.path.basename(input_files[i])))
    if possible_extension != input_extension:
        input_extension = possible_extension
    for chromosome in chromosomes:
        possible_sample = re.sub("(_$)", "", re.sub(chromosome, "", possible_sample))
        if input_sample != possible_sample:
            input_sample = possible_sample


rule all:
    input:
        os.path.join(os.path.abspath(output_folder), 'process_gff', 'processed_gff.gff'),
        expand(os.path.join(os.path.abspath(output_folder), 'beagle_impute_input_file', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'snpEff_input_file', input_sample+'_{chromosome}.html'), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'snpEff_input_file', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_wide', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
        os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_wide', '{project_name}.txt'.format(project_name=project_name))


include: './rules/bash/grep_gff.smk'

include: './rules/java/beagle_impute_input_file.smk'
include: './rules/java/snpEff_input_file.smk'

include: './rules/bash/grep_input_effects.smk'

include: './rules/python/generate_Allele_Catalog.smk'
include: './rules/python/generate_Allele_Catalog_wide.smk'
include: './rules/python/combine_Allele_Catalog_wide.smk'
