import re


class FunctionalEffectAnnotationObject:
    # These are the effects we need to focus on.
    effects = [
        'frameshift_variant',
        'exon_loss_variant',
        'duplication',
        'inversion',
        'feature_ablation',
        'gene_fusion',
        'rearranged_at_DNA_level',
        'missense_variant',
        'protein_protein_contact',
        'structural_interaction_variant',
        'rare_amino_acid_variant',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'stop_lost',
        'start_lost',
        'stop_gained',
        'inframe_insertion',
        'disruptive_inframe_insertion',
        'inframe_deletion',
        'disruptive_inframe_deletion'
    ]

    def __init__(self, allele, functional_effect="Alt", amino_acid_change=""):
        self.allele = allele
        self.functional_effects = []
        self.functional_effects.append(functional_effect)
        self.amino_acid_changes = []
        if amino_acid_change != "":
            self.amino_acid_changes.append(amino_acid_change)

    def get_allele(self):
        return self.allele

    def get_functional_effects(self):
        return self.functional_effects

    def get_amino_acid_changes(self):
        return self.amino_acid_changes

    def get_functional_effects_string(self):
        # Return functional effect joined string
        return '&'.join(self.functional_effects)

    def get_amino_acid_changes_string(self):
        # Create amino acid changes joined string
        amino_acid_changes_string = '&'.join(self.amino_acid_changes)
        # Replace amino acid three letters code to amino acid one letter code
        amino_acid_changes_string = re.sub('Gly', 'G', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Pro', 'P', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Val', 'V', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Leu', 'L', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Met', 'M', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Cys', 'C', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Phe', 'F', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Tyr', 'Y', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Trp', 'W', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('His', 'H', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Lys', 'K', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Arg', 'R', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Gln', 'Q', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Asp', 'D', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Ser', 'S', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Thr', 'T', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Asn', 'N', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Ile', 'I', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Glu', 'E', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('Ala', 'A', amino_acid_changes_string)
        amino_acid_changes_string = re.sub('p\\.', '', amino_acid_changes_string)
        # Return amino acid changes joined string
        return amino_acid_changes_string

    def get_functional_effect_annotation_string(self):
        # If functional effects and amino acid changes are not empty, return everything
        # If functional effects are not empty but amino acid changes are empty, return allele and functional effects
        # If functional effects and amino acid changes are empty, return allele with "Alt"
        if len(self.functional_effects) > 0 and len(self.amino_acid_changes) > 0:
            return str(
                self.allele) + "|" + self.get_functional_effects_string() + "|" + self.get_amino_acid_changes_string()
        elif len(self.functional_effects) > 0 and len(self.amino_acid_changes) == 0:
            return str(self.allele) + "|" + self.get_functional_effects_string()
        else:
            return str(self.allele) + "|Alt"

    def append_functional_effect(self, functional_effect):
        # If the functional effect is in the functional effect array then it is not going to be appended
        if (functional_effect != "") and (functional_effect not in self.functional_effects):
            # If the functional effect is in the functional effect joined string then it is not going to be appended
            if re.search(functional_effect, self.get_functional_effects_string()) is None:
                # It is only going to be appended if it is the effect we are interested in
                if any(
                        [(re.search(e, functional_effect) is not None)
                         for e in FunctionalEffectAnnotationObject.effects]
                ):
                    self.functional_effects.append(functional_effect)
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def append_amino_acid_changes(self, amino_acid_change):
        amino_acid_change_pattern = ""
        # If the amino acid change is found in the amino acid change array then it is not going to be appended
        if (amino_acid_change != "") and (amino_acid_change not in self.amino_acid_changes):
            # If the amino acid change is found in the amino acid change joined string then it is not to be appended
            try:
                amino_acid_change_pattern = re.sub("\\\\", "\\\\\\\\", amino_acid_change)
                amino_acid_change_pattern = re.sub("\\+", "\\+", amino_acid_change_pattern)
                amino_acid_change_pattern = re.sub("\\*", "\\*", amino_acid_change_pattern)
                amino_acid_change_pattern = re.sub("\\?", "\\?", amino_acid_change_pattern)
                amino_acid_change_pattern = re.sub("\\(", "\\(", amino_acid_change_pattern)
                amino_acid_change_pattern = re.sub("\\)", "\\)", amino_acid_change_pattern)
                amino_acid_change_pattern = re.sub("\\[", "\\[", amino_acid_change_pattern)
                if (amino_acid_change_pattern != "") and \
                        (re.search(amino_acid_change_pattern, self.get_amino_acid_changes_string()) is None):
                    self.amino_acid_changes.append(amino_acid_change)
                    return True
                else:
                    return False
            except Exception as e:
                print("FunctionalEffectAnnotationObject - append_amino_acid_changes function error !!!")
                print(
                    amino_acid_change_pattern + " extracted from " + amino_acid_change +
                    " is causing error when it match with " + self.get_amino_acid_changes_string()
                )
                return False
        else:
            return False
