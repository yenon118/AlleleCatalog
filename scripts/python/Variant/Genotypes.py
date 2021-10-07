import re

from Variant.GenotypeObject import GenotypeObject


class Genotypes:

    def __init__(self, reference_allele, annotated_reference_allele, alternate_alleles, annotated_alternate_alleles,
                 sample_array, indexes_string_array):
        self.reference_allele = reference_allele
        self.annotated_reference_allele = annotated_reference_allele
        self.alternate_alleles = alternate_alleles
        self.annotated_alternate_alleles = annotated_alternate_alleles
        self.sample_array = sample_array
        self.indexes_string_array = indexes_string_array
        self.indexes_array = self.generate_indexes_array_from_indexes_string_array()
        self.genotype_object_dictionary = self.generate_genotype_object_dictionary()

    # Get genotype object dictionary
    def get_genotype_object_dictionary(self):
        return self.genotype_object_dictionary

    # Generate indexes array
    def generate_indexes_array_from_indexes_string_array(self):
        indexes_array = []
        for i in range(len(self.indexes_string_array)):
            indexes_array.append(re.split('/|\\|', str(re.sub(":.*", "", self.indexes_string_array[i]))))
        return indexes_array

    # Generate genotype object dictionary based on samples and all related information
    def generate_genotype_object_dictionary(self):
        genotype_object_dictionary = {}
        for i in range(len(self.sample_array)):
            if self.sample_array[i] not in genotype_object_dictionary.keys():
                genotype_object_dictionary[self.sample_array[i]] = GenotypeObject(
                    self.reference_allele,
                    self.annotated_reference_allele,
                    self.alternate_alleles,
                    self.annotated_alternate_alleles,
                    self.sample_array[i],
                    self.indexes_array[i]
                )
        return genotype_object_dictionary
