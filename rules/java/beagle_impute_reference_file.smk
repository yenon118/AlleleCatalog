rule beagle_impute_reference_file:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}', reference_sample+'_{chromosome}'+reference_extension)
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'beagle_impute_reference_file', reference_sample+'_{chromosome}'+reference_extension),
        out_tmp_dir = temp(directory(os.path.join(os.path.abspath(output_folder), 'beagle_impute_reference_file', 'tmp', '{chromosome}')))
    params:
        out_file_param = os.path.join(os.path.abspath(output_folder), 'beagle_impute_reference_file', reference_sample+'_{chromosome}')
    log:
        os.path.join(os.path.abspath(output_folder), 'beagle_impute_reference_file_log', reference_sample+'_{chromosome}.log')
    threads: threads
    resources:
        memory = memory
    shell:
        """
        mkdir -p {output.out_tmp_dir};
        java -Xmx{resources.memory}G -Djava.io.tmpdir={output.out_tmp_dir} -jar {workflow_path}/tools/beagle.jar gt={input.in_file} out={params.out_file_param} nthreads={threads} 2> {log}
        {workflow_path}/tools/bgzip -d {output.out_file}.gz
        """
