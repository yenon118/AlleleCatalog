rule grep_gff:
    input:
        in_file = gff_file
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'process_gff', 'processed_gff.gff')
    priority: 90
    shell:
        """
        zgrep -e "{gff_category}" -e "{gff_key}" {input.in_file} > {output.out_file}
        """
