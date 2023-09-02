rule beagle_impute_input_file:
	input:
		in_file = os.path.join(os.path.abspath(input_folder), input_sample+'_{chromosome}'+input_extension),
		ref_file = os.path.join(os.path.abspath(reference_folder), reference_sample+'_{chromosome}'+reference_extension)
	output:
		out_file = os.path.join(os.path.abspath(output_folder), 'beagle_impute_input_file', input_sample+'_{chromosome}'+input_extension),
		out_tmp_dir = temp(directory(os.path.join(os.path.abspath(output_folder), 'beagle_impute_input_file', 'tmp', '{chromosome}')))
	params:
		beagle_window = beagle_window,
		out_file_param = os.path.join(os.path.abspath(output_folder), 'beagle_impute_input_file', input_sample+'_{chromosome}')
	threads: threads
	resources:
		memory = memory
	shell:
		"""
		mkdir -p {output.out_tmp_dir};

		java -Xmx{resources.memory}G \
		-Djava.io.tmpdir={output.out_tmp_dir} \
		-jar {workflow_path}/tools/beagle.jar \
		gt={input.in_file} \
		ref={input.ref_file} \
		out={params.out_file_param} \
		window={params.beagle_window} \
		nthreads={threads}
		"""
