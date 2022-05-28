rule generate_functional_effect_data:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension)
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'generate_functional_effect_data', input_sample+'_{chromosome}.txt')
    log:
        os.path.join(os.path.abspath(output_folder), 'generate_functional_effect_data_log', input_sample+'_{chromosome}.log')
    threads: threads
    shell:
        """
        python3 {workflow_path}/scripts/python/generate_functional_effect_data.py -i {input.in_file} -o {output.out_file} 2>&1 > {log}
        """
