import sys
import os
import re


project_name = config['project_name']
workflow_path = config['workflow_path']
input_files = config['input_files']
reference_files = config['reference_files']
chromosomes = config['chromosomes']
ancestors = config['ancestors']
reference_file = config['reference_file']
genome_version = config['genome_version']
gff_file = config['gff_file']
gff_category = config['gff_category']
gff_key = config['gff_key']
output_folder = config['output_folder']
memory = config['memory']
threads = config['threads']


samples = []
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
    for chromosome in chromosomes:
        possible_sample = re.sub("(_$)", "", re.sub(chromosome, "", possible_sample))
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
    for chromosome in chromosomes:
        possible_sample = re.sub("(_$)", "", re.sub(chromosome, "", possible_sample))
        if reference_sample != possible_sample:
            reference_sample = possible_sample

for i, j in zip(input_samples, reference_samples):
    samples.append(i)
    samples.append(j)


rule all:
    input:
        os.path.join(os.path.abspath(output_folder), 'process_gff', 'processed_gff.gff'),
        expand(os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}'), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}', reference_sample+'_{chromosome}'+reference_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'beagle_impute_reference_file', reference_sample+'_{chromosome}'+reference_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'beagle_impute_input_file', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'snpEff_reference_file', reference_sample+'_{chromosome}.html'), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'snpEff_reference_file', reference_sample+'_{chromosome}'+reference_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'snpEff_input_file', input_sample+'_{chromosome}.html'), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'snpEff_input_file', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'grep_effects', reference_sample+'_{chromosome}.txt'), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog', '{sample}.txt'), sample=samples),
        expand(os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_wide', '{sample}.txt'), sample=samples),
        os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_wide', '{project_name}.txt'.format(project_name=project_name)),
        expand(os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_wide_split', '{project_name}_{{chromosome}}.txt'.format(project_name=project_name)), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'generate_Ancestry_Binary', '{project_name}_{{chromosome}}.txt'.format(project_name=project_name)), chromosome=chromosomes),
        os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_final', '{project_name}.txt'.format(project_name=project_name))


include: './rules/bash/grep_gff.smk'

include: './rules/python/align_input_VCF_based_on_reference_VCF.smk'

include: './rules/java/beagle_impute_reference_file.smk'
include: './rules/java/beagle_impute_input_file.smk'
include: './rules/java/snpEff_reference_file.smk'
include: './rules/java/snpEff_input_file.smk'

include: './rules/bash/grep_reference_effects.smk'
include: './rules/bash/grep_input_effects.smk'

include: './rules/python/generate_Allele_Catalog.smk'
include: './rules/python/generate_Allele_Catalog_wide.smk'
include: './rules/python/combine_Allele_Catalog_wide.smk'

include: './rules/bash/grep_header_and_chromosome.smk'

include: './rules/python/generate_Ancestry_Binary.smk'
include: './rules/python/combine_Allele_Catalog_final.smk'
