rule generate_Allele_Catalog_wide:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog', input_sample+'_{chromosome}.txt')
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_wide', input_sample+'_{chromosome}.txt')
    log:
        os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_wide_log', input_sample+'_{chromosome}.log')
    shell:
        """
        python3 {workflow_path}/scripts/python/generate_Allele_Catalog_wide.py -i {input.in_file} -o {output.out_file} 2>&1 > {log}
        """
