rule generate_Allele_Catalog:
    input:
        reference_file = reference_file,
        in_file = os.path.join(os.path.abspath(output_folder), 'grep_effects', '{sample}.txt'),
        gff_file = os.path.join(os.path.abspath(output_folder), 'process_gff', 'processed_gff.gff')
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog', '{sample}.txt')
    log:
        os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_log', '{sample}.log')
    message: "Executing generate_Allele_Catalog rule with {threads} threads."
    threads: threads
    shell:
        """
        python3 {workflow_path}/scripts/python/generate_Allele_Catalog.py -i {input.in_file} -r {input.reference_file} -g {input.gff_file} -o {output.out_file} -t {threads} -c {gff_category} -k {gff_key} 2>&1 > {log}
        """
