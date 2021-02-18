rule combine_Allele_Catalog_wide:
    input:
        in_file = expand(os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_wide', '{sample}.txt'), sample=samples)
    params:
        ' '.join(['-i '+os.path.join(os.path.abspath(output_folder), 'generate_Allele_Catalog_wide', '{sample}.txt'.format(sample=sample))  for sample in samples])
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_wide', '{project_name}.txt'.format(project_name=project_name))
    log:
        os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_wide_log', '{project_name}.log'.format(project_name=project_name))
    shell:
        """
        python3 {workflow_path}/scripts/python/combine_Allele_Catalog_wide.py {params} -o {output.out_file} 2> {log}
        """
