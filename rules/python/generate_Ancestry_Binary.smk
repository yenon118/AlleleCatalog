rule generate_Ancestry_Binary:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_wide_split', '{project_name}_{{chromosome}}.txt'.format(project_name=project_name))
    params:
        ancestors_params = ' '.join(['-a '+str(ancestor) for ancestor in ancestors]) if len(ancestors) > 0 else ""
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'generate_Ancestry_Binary', '{project_name}_{{chromosome}}.txt'.format(project_name=project_name))
    threads: threads
    shell:
        """
        python3 {workflow_path}/scripts/python/generate_Ancestry_Binary.py -i {input.in_file} -o {output.out_file} -t {threads} {params.ancestors_params}
        """
