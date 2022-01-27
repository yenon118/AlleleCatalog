import sys
import os
import re

project_name = config['project_name']
workflow_path = config['workflow_path']
input_files = config['input_files']
chromosomes = config['chromosomes']
ancestors = config['ancestors']
reference_file = config['reference_file']
gff_file = config['gff_file']
gff_category = config['gff_category']
gff_key = config['gff_key']
output_folder = config['output_folder']
threads = config['threads']


input_folder = ''
input_sample = ''
input_samples = []
input_extension = ''


for i in range(len(input_files)):
    if os.path.dirname(input_files[i]) != input_folder:
        input_folder = os.path.dirname(input_files[i])
    possible_sample = re.sub('(\\.vcf.*)', '', str(os.path.basename(input_files[i])))
    if not possible_sample in input_samples:
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
        expand(os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_wide', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
        os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_wide', '{project_name}.txt'.format(project_name=project_name)),
        expand(os.path.join(os.path.abspath(output_folder), 'generate_Ancestry_Binary', input_sample+'_{chromosome}.txt'), chromosome=chromosomes),
        os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_final', '{project_name}.txt'.format(project_name=project_name))


rule process_gff:
    input:
        in_file = gff_file
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'process_gff', 'processed_gff.gff')
    priority: 90
    shell:
        """
        zgrep -e "{gff_category}" -e "{gff_key}" {input.in_file} > {output.out_file}
        """


rule grep_effects:
    input:
        in_file = os.path.join(os.path.abspath(input_folder), input_sample+'_{chromosome}'+input_extension)
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension)
    priority: 70
    shell:
        """
        zgrep -e "^#CHROM" -e "frameshift_variant" -e "exon_loss_variant" -e "duplication" -e "inversion" -e "feature_ablation" -e "gene_fusion" -e "rearranged_at_DNA_level" -e "missense_variant" -e "protein_protein_contact" -e "structural_interaction_variant" -e "rare_amino_acid_variant" -e "splice_acceptor_variant" -e "splice_donor_variant" -e "stop_lost" -e "start_lost" -e "stop_gained" -e "inframe_insertion" -e "disruptive_inframe_insertion" -e "inframe_deletion" -e "disruptive_inframe_deletion" {input.in_file} > {output.out_file}
        """

include: './../rules/python/generate_Allele_Catalog.smk'
include: './../rules/python/generate_Allele_Catalog_wide.smk'
include: './../rules/python/combine_Allele_Catalog_wide.smk'
include: './../rules/python/generate_Ancestry_Binary.smk'
include: './../rules/python/combine_Allele_Catalog_final.smk'
