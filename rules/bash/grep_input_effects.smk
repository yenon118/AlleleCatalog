rule grep_input_effects:
	input:
		in_file = os.path.join(os.path.abspath(output_folder), 'snpEff_input_file', input_sample+'_{chromosome}'+input_extension)
	output:
		out_file = os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension)
	shell:
		"""
		zgrep -e "^#" -e "^#CHROM" -e "frameshift_variant" -e "exon_loss_variant" -e "duplication" -e "inversion" \
		-e "feature_ablation" -e "gene_fusion" -e "rearranged_at_DNA_level" -e "missense_variant" -e "protein_protein_contact" \
		-e "structural_interaction_variant" -e "rare_amino_acid_variant" -e "splice_acceptor_variant" -e "splice_donor_variant" \
		-e "stop_lost" -e "start_lost" -e "stop_gained" -e "inframe_insertion" -e "disruptive_inframe_insertion" -e "inframe_deletion" \
		-e "disruptive_inframe_deletion" {input.in_file} | grep -v -e "DEL" -e "INS" | bgzip > {output.out_file}
		"""
