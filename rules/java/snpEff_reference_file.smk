rule snpEff_reference_file:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'beagle_impute_reference_file', reference_sample+'_{chromosome}'+reference_extension)
    output:
        summary_file = os.path.join(os.path.abspath(output_folder), 'snpEff_reference_file', reference_sample+'_{chromosome}.html'),
        out_file = os.path.join(os.path.abspath(output_folder), 'snpEff_reference_file', reference_sample+'_{chromosome}'+reference_extension),
        out_tmp_dir = temp(directory(os.path.join(os.path.abspath(output_folder), 'snpEff_reference_file', 'tmp', '{chromosome}')))
    params:
        genome_version = genome_version
    log:
        os.path.join(os.path.abspath(output_folder), 'snpEff_reference_file_log', reference_sample+'_{chromosome}.log')
    threads: threads
    resources:
        memory = memory
    shell:
        """
        mkdir -p {output.out_tmp_dir};
        java -Xmx{resources.memory}G -Djava.io.tmpdir={output.out_tmp_dir} -jar {workflow_path}/tools/snpEff/snpEff.jar -s {output.summary_file} -v {params.genome_version} {input.in_file} 2> {log} 1> {output.out_file}
        """
