rule combine_Allele_Catalog_final:
    input:
        in_file = expand(os.path.join(os.path.abspath(output_folder), 'generate_Ancestry_Binary', '{project_name}_{{chromosome}}.txt'.format(project_name=project_name)), chromosome=chromosomes)
    params:
        ' '.join(['-i '+os.path.join(os.path.abspath(output_folder), 'generate_Ancestry_Binary', '{project_name}_{chromosome}.txt'.format(project_name=project_name, chromosome=chromosome))  for chromosome in chromosomes])
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_final', '{project_name}.txt'.format(project_name=project_name))
    log:
        os.path.join(os.path.abspath(output_folder), 'combine_Allele_Catalog_final_log', '{project_name}.log'.format(project_name=project_name))
    shell:
        """
        python3 {workflow_path}/scripts/python/combine_Allele_Catalog_final.py {params} -o {output.out_file} 2> {log}
        """
