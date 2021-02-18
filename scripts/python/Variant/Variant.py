import re

from Variant.FunctionalEffectAnnotationObject import FunctionalEffectAnnotationObject
from Variant.FunctionalEffectAnnotations import FunctionalEffectAnnotations
from Variant.Genotypes import Genotypes
from Variant.GenotypeObject import GenotypeObject


class Variant:

    def __init__(self, header, line):
        self.header = str(header).strip()
        self.line = str(line).strip()
        self.chromosome = self.acquire_chromosome()
        self.position = self.acquire_position()
        self.reference_allele, self.annotated_reference_allele = self.acquire_reference_allele()
        self.alternate_alleles, self.annotated_alternate_alleles, self.gene = self.acquire_alternate_alleles_and_gene()
        self.genotypes_dictionary = self.acquire_genotypes_dictionary()

    # Print the info in this class
    def print(self):
        print(
            self.chromosome + " " +
            self.position + " " +
            self.gene + " " +
            self.reference_allele + " " +
            self.annotated_reference_allele + " " +
            str(self.alternate_alleles) + " " +
            str(self.annotated_alternate_alleles) + " " +
            str(self.genotypes_dictionary)
        )

    # Convert to allele catalog format string
    def to_allele_catalog_string(self):
        array = []
        for key in self.genotypes_dictionary.keys():
            array.append(
                self.genotypes_dictionary[key]["Sample"] + "\t" +
                self.chromosome + "\t" +
                self.gene + "\t" +
                self.position + "\t" +
                self.genotypes_dictionary[key]["Genotype"] + "\t" +
                self.genotypes_dictionary[key]["Annotated_Genotype"]
            )
        return "\n".join(array)

    def set_gene(self, gene):
        self.gene = gene

    def get_chromosome(self):
        return self.chromosome

    def get_position(self):
        return self.position

    def get_gene(self):
        return self.gene

    def get_genotypes_dictionary(self):
        return self.genotypes_dictionary

    def acquire_chromosome(self):
        # Get chromosome from line
        return self.line.split("\t")[0].strip()

    def acquire_position(self):
        # Get position from line
        return self.line.split("\t")[1].strip()

    def acquire_reference_allele(self):
        # Get reference allele from line
        reference_allele = self.line.split("\t")[3].strip()
        # Return reference allele and annotated reference allele
        return reference_allele, reference_allele + "|Ref"

    def acquire_alternate_alleles_and_gene(self):
        # Initialize gene as empty string
        gene = ""

        # Get alternate alleles in array format
        alternate_alleles = self.line.split("\t")[4].split(",")

        # Look for functional effect annotation string in the info column
        info = self.line.split("\t")[7].split(";")
        functional_effect_annotation_string = ""
        for i in range(len(info)):
            if info[i].startswith("EFF="):
                functional_effect_annotation_string = re.sub("EFF=", "", info[i])
                functional_effect_annotation_string = functional_effect_annotation_string.strip()
            if info[i].startswith("ANN="):
                functional_effect_annotation_string = re.sub("ANN=", "", info[i])
                functional_effect_annotation_string = functional_effect_annotation_string.strip()
                break

        # Create functional effect annotation object
        # This functional effect annotation object can automatically process functional effect annotation string
        functional_effect_annotation = FunctionalEffectAnnotations(
            alternate_alleles,
            gene,
            functional_effect_annotation_string
        )

        # Get annotated alternate alleles and gene
        # The gene here is based on SnpEff output. It might not be accurate
        annotated_alternate_alleles = functional_effect_annotation.get_annotated_alternate_alleles()
        gene = functional_effect_annotation.get_gene()

        return alternate_alleles, annotated_alternate_alleles, gene

    def acquire_genotypes_dictionary(self):
        # Initialize variable
        genotypes_dictionary = {}

        # Get samples and genotype indexes
        sample_array = self.header.split("\t")[10:]
        indexes_string_array = self.line.split("\t")[10:]

        # Use genotypes class to generate genotype object dictionary
        genotypes = Genotypes(
            self.reference_allele,
            self.annotated_reference_allele,
            self.alternate_alleles,
            self.annotated_alternate_alleles,
            sample_array,
            indexes_string_array
        )

        # Get the genotype object dictionary
        genotype_object_dictionary = genotypes.get_genotype_object_dictionary()

        # Re-organize the genotype object dictionary to make it more accessible
        for key in genotype_object_dictionary.keys():
            if key not in genotypes_dictionary.keys():
                genotypes_dictionary[key] = {
                    "Sample": genotype_object_dictionary[key].get_sample(),
                    "Genotype": genotype_object_dictionary[key].get_genotype(),
                    "Annotated_Genotype": genotype_object_dictionary[key].get_annotated_genotype()
                }

        return genotypes_dictionary
