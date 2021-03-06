rule grep_reference_effects:
    input:
        in_file = os.path.join(os.path.abspath(output_folder), 'snpEff_reference_file', reference_sample+'_{chromosome}'+reference_extension)
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'grep_effects', reference_sample+'_{chromosome}.txt')
    shell:
        """
        zgrep -e "^#CHROM" -e "frameshift_variant" -e "exon_loss_variant" -e "duplication" -e "inversion" -e "feature_ablation" -e "gene_fusion" -e "rearranged_at_DNA_level" -e "missense_variant" -e "protein_protein_contact" -e "structural_interaction_variant" -e "rare_amino_acid_variant" -e "splice_acceptor_variant" -e "splice_donor_variant" -e "stop_lost" -e "start_lost" -e "stop_gained" -e "inframe_insertion" -e "disruptive_inframe_insertion" -e "inframe_deletion" -e "disruptive_inframe_deletion" {input.in_file} > {output.out_file}
        """
