import re


class GenotypeObject:

    def __init__(self, reference_allele, annotated_reference_allele, alternate_alleles, annotated_alternate_alleles,
                 sample, genotype_indexes):
        self.reference_allele = reference_allele
        self.annotated_reference_allele = annotated_reference_allele
        self.alternate_alleles = alternate_alleles
        self.annotated_alternate_alleles = annotated_alternate_alleles
        self.sample = sample
        self.genotype_indexes = genotype_indexes
        self.genotype = self.map_genotype_indexes_for_genotype()
        self.annotated_genotype = self.map_genotype_indexes_for_annotated_genotype()

    # Get sample
    def get_sample(self):
        return self.sample

    # Get first genotype
    def get_genotype0(self):
        if len(self.genotype) > 0:
            return self.genotype[0]
        else:
            return ""

    # Get second genotype
    def get_genotype1(self):
        if len(self.genotype) > 1:
            return self.genotype[1]
        else:
            return ""

    # Get the genotype that is in alternate alleles if possible
    def get_genotype(self):
        for g in self.genotype:
            if g in self.alternate_alleles:
                return g
        if len(self.genotype) > 1:
            return self.genotype[1]
        elif len(self.genotype) > 0:
            return self.genotype[0]
        else:
            return ""

    # Get joined genotype
    def get_genotypes(self):
        return '/'.join(self.genotype)

    # Get first annotated genotype
    def get_annotated_genotype0(self):
        if len(self.annotated_genotype) > 0:
            return self.annotated_genotype[0]
        else:
            return ""

    # Get second annotated genotype
    def get_annotated_genotype1(self):
        if len(self.annotated_genotype) > 1:
            return self.annotated_genotype[1]
        else:
            return ""

    # Get annotated genotype that matches any in the annotated alternate alleles if possible
    def get_annotated_genotype(self):
        for g in self.annotated_genotype:
            if g in self.annotated_alternate_alleles:
                return g
        if len(self.annotated_genotype) > 1:
            return self.annotated_genotype[1]
        elif len(self.annotated_genotype) > 0:
            return self.annotated_genotype[0]
        else:
            return ""

    # Get joined annotated genotype
    def get_annotated_genotypes(self):
        return '/'.join(self.annotated_genotype)

    # Map genotype indexes to get genotype
    def map_genotype_indexes_for_genotype(self):
        genotype = []
        for i in range(len(self.genotype_indexes)):
            if int(self.genotype_indexes[i]) == 0:
                genotype.append(self.reference_allele)
            elif int(self.genotype_indexes[i]) > 0:
                genotype.append(self.alternate_alleles[int(self.genotype_indexes[i]) - 1])
        return genotype

    # Map genotype indexes to get annotated genotype
    def map_genotype_indexes_for_annotated_genotype(self):
        annotated_genotype = []
        for i in range(len(self.genotype_indexes)):
            if int(self.genotype_indexes[i]) == 0:
                annotated_genotype.append(self.annotated_reference_allele)
            elif int(self.genotype_indexes[i]) > 0:
                annotated_genotype.append(self.annotated_alternate_alleles[int(self.genotype_indexes[i]) - 1])
        return annotated_genotype
