rule generate_Allele_Catalog:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'generate_genotype_data', input_sample+'_{chromosome}.txt'),
        metadata_file = metadata_file,
        gff_file = gff_file
    params:
        chromosome = '{chromosome}'
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog', input_sample+'_{chromosome}.txt')
    log:
        os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_log', input_sample+'_{chromosome}.log')
    threads: threads
    shell:
        """
        python3 {workflow_path}/scripts/python/generate_Allele_Catalog.py \
        -i {input.in_file} \
        -m {input.metadata_file} \
        -g {input.gff_file} \
        -o {output.out_file} \
        -c {params.chromosome} \
        -a {gff_category} \
        -k {gff_key} 2>&1 > {log}
        """
