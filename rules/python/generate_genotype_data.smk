rule generate_genotype_data:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension),
        functional_effect_file = os.path.join(os.path.abspath(output_folder), 'generate_functional_effect_data', input_sample+'_{chromosome}.txt'),
        imputation_file = os.path.join(os.path.abspath(output_folder), 'generate_imputation_data', input_sample+'_{chromosome}.txt')
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'generate_genotype_data', input_sample+'_{chromosome}.txt')
    log:
        os.path.join(os.path.abspath(output_folder), 'generate_genotype_data_log', input_sample+'_{chromosome}.log')
    threads: threads
    shell:
        """
        python3 {workflow_path}/scripts/python/generate_genotype_data.py \
        -i {input.in_file} \
        -f {input.functional_effect_file} \
        -m {input.imputation_file} \
        -o {output.out_file} 2>&1 > {log}
        """
