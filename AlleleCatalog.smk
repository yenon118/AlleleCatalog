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
gff_file = config['gff_file']
gff_category = config['gff_category']
gff_key = config['gff_key']
output_folder = config['output_folder']
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


print(input_folder)
print(input_sample)
print(input_samples)
print(input_extension)

print(reference_folder)
print(reference_sample)
print(reference_samples)
print(reference_extension)

print(output_folder)
print(samples)


rule all:
    input:
        os.path.join(os.path.abspath(output_folder), 'process_gff', 'processed_gff.gff'),
        expand(os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}'), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}', reference_sample+'_{chromosome}'+reference_extension), chromosome=chromosomes)


include: './rules/bash/grep_gff.smk'

include: './rules/python/align_input_VCF_based_on_reference_VCF.smk'
