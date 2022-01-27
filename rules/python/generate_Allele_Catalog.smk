rule generate_Allele_Catalog:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension),
        gff_file = os.path.join(os.path.abspath(output_folder), 'process_gff', 'processed_gff.gff')
    params:
        '-r '+reference_file if (os.path.exists(reference_file) and os.path.isfile(reference_file)) else ""
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog', input_sample+'_{chromosome}.txt')
    log:
        os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_log', input_sample+'_{chromosome}.log')
    message: "Executing generate_Allele_Catalog rule with {threads} threads."
    threads: threads
    shell:
        """
        python3 {workflow_path}/scripts/python/generate_Allele_Catalog.py -i {input.in_file} {params} -g {input.gff_file} -o {output.out_file} -t {threads} -c {gff_category} -k {gff_key} 2>&1 > {log}
        """
