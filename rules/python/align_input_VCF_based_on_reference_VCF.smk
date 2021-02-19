rule align_input_VCF_based_on_reference_VCF:
    input:
        in_file = os.path.join(os.path.abspath(input_folder), input_sample+'_{chromosome}'+input_extension),
        ref_file = os.path.join(os.path.abspath(reference_folder), reference_sample+'_{chromosome}'+reference_extension)
    output:
        out_dir = directory(os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}')),
        out_file_in = os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}', input_sample+'_{chromosome}'+input_extension),
        out_file_ref = os.path.join(os.path.abspath(output_folder), 'align_input_VCF_based_on_reference_VCF', '{chromosome}', reference_sample+'_{chromosome}'+reference_extension)
    shell:
        """
        mkdir -p {output.out_dir};
        python3 {workflow_path}/scripts/python/align_input_VCF_based_on_reference_VCF.py -i {input.in_file} -r {input.ref_file} -o {output.out_dir}
        """
