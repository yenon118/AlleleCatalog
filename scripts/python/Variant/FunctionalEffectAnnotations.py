import re

from Variant.FunctionalEffectAnnotationObject import FunctionalEffectAnnotationObject


class FunctionalEffectAnnotations:

    def __init__(self, alternate_alleles, gene="", functional_effect_annotation_string=""):
        self.alternate_alleles = alternate_alleles
        self.gene = gene
        self.functional_effect_annotation_string = str(functional_effect_annotation_string)
        if self.functional_effect_annotation_string != "":
            self.functional_effect_annotation_string_array = self.functional_effect_annotation_string.split(",")
        else:
            self.functional_effect_annotation_string_array = []
        self.annotated_alternate_alleles = self.generate_annotated_alternate_alleles()

    def set_gene(self, gene):
        self.gene = gene

    def get_gene(self):
        return self.gene

    def get_annotated_alternate_alleles(self):
        return self.annotated_alternate_alleles

    def generate_annotated_alternate_alleles(self):
        # Initialize variables
        annotated_alternate_alleles_dict = {}
        annotated_alternate_alleles = []
        genes = []

        # If functional effect annotation string is not empty, process the functional effect string
        # to get annotated alternate alleles and genes.
        # If functional effect annotation string is empty, use alternate alleles to generate
        # annotated alternate alleles.
        if self.functional_effect_annotation_string != "":
            # For every functional effect annotation string array, look into the
            # allele, functional effect, gene, and amino acid change.
            # Use allele, functional effect, and amino acid change to create
            # annotated alternate alleles dictionary.
            for i in range(len(self.functional_effect_annotation_string_array)):
                single_functional_effect_annotation = re.split('/|\\|', str(self.functional_effect_annotation_string_array[i]))
                # Check if it is a primary transcript.
                # We only consider primary transcript.
                pattern = re.sub("\\\\", "\\\\\\\\", single_functional_effect_annotation[3])
                pattern = re.sub("\\+", "\\+", pattern)
                pattern = re.sub("\\*", "\\*", pattern)
                pattern = re.sub("\\?", "\\?", pattern)
                pattern = re.sub("\\(", "\\(", pattern)
                pattern = re.sub("\\)", "\\)", pattern)
                pattern = re.sub("\\[", "\\[", pattern)
                if re.search("\\.1", re.sub(pattern, "", single_functional_effect_annotation[6])) is not None:
                    allele = single_functional_effect_annotation[0]
                    functional_effect = single_functional_effect_annotation[1]
                    amino_acid_change = single_functional_effect_annotation[10]
                    # If a gene is not in the gene array, append that gene
                    if re.sub("(.*-)", "", single_functional_effect_annotation[3]) not in genes:
                        genes.append(
                            re.sub("(.*-)", "", single_functional_effect_annotation[3])
                        )
                        self.gene = '&'.join(genes)
                    # If allele is not in the annotated alternate alleles dictionary keys, create a new record
                    # in the annotated alternate alleles dictionary with allele as the key.
                    # If allele is already in the annotated alternate alleles dictionary keys, use the
                    # key to add new information into the record in the dictionary.
                    if allele not in annotated_alternate_alleles_dict.keys():
                        if amino_acid_change != "":
                            annotated_alternate_alleles_dict[allele] = FunctionalEffectAnnotationObject(
                                allele, functional_effect, amino_acid_change
                            )
                        else:
                            annotated_alternate_alleles_dict[allele] = FunctionalEffectAnnotationObject(
                                allele, functional_effect
                            )
                    else:
                        flag = annotated_alternate_alleles_dict[allele].append_functional_effect(functional_effect)
                        if amino_acid_change != "":
                            flag = annotated_alternate_alleles_dict[allele].append_amino_acid_changes(amino_acid_change)
            # From the annotated alternate alleles dictionary, use alternate alleles as key to
            # generate annotated alternate alleles which has order just like the alternate alleles.
            for i in range(len(self.alternate_alleles)):
                if self.alternate_alleles[i] in annotated_alternate_alleles_dict.keys():
                    annotated_alternate_alleles.append(
                        annotated_alternate_alleles_dict[
                            self.alternate_alleles[i]].get_functional_effect_annotation_string()
                    )
                else:
                    annotated_alternate_alleles.append(
                        FunctionalEffectAnnotationObject(
                            self.alternate_alleles[i]
                        ).get_functional_effect_annotation_string()
                    )
        else:
            for i in range(len(self.alternate_alleles)):
                annotated_alternate_alleles.append(
                    FunctionalEffectAnnotationObject(
                        self.alternate_alleles[i]
                    ).get_functional_effect_annotation_string()
                )
        return annotated_alternate_alleles
