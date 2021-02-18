rule grep_header_and_chromosome:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_wide', '{project_name}.txt'.format(project_name=project_name))
    params:
        chromosome = "{chromosome}"
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_wide_split', '{project_name}_{{chromosome}}.txt'.format(project_name=project_name))
    shell:
        """
        grep -e "^Classification" -e "{params.chromosome}" {input.in_file} > {output.out_file}
        """
