rule generate_Ancestry_Binary:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_wide', input_sample+'_{chromosome}.txt')
    params:
        ancestors_params = ' '.join(['-a "'+str(ancestor)+'"' for ancestor in ancestors]) if len(ancestors) > 0 else ""
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'generate_Ancestry_Binary', input_sample+'_{chromosome}.txt')
    threads: threads
    shell:
        """
        python3 {workflow_path}/scripts/python/generate_Ancestry_Binary.py -i {input.in_file} -o {output.out_file} -t {threads} {params.ancestors_params}
        """
