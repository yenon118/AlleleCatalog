#!/usr/bin/env python3

import sys
import os
import re
import argparse
import pathlib
import gzip
import threading

from joblib import Parallel, delayed


# These are the effects we need to focus on.
EFFECTS = [
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


def clean_amino_acid_change_string(amino_acid_change):
    amino_acid_change = re.sub('Gly', 'G', amino_acid_change)
    amino_acid_change = re.sub('Pro', 'P', amino_acid_change)
    amino_acid_change = re.sub('Val', 'V', amino_acid_change)
    amino_acid_change = re.sub('Leu', 'L', amino_acid_change)
    amino_acid_change = re.sub('Met', 'M', amino_acid_change)
    amino_acid_change = re.sub('Cys', 'C', amino_acid_change)
    amino_acid_change = re.sub('Phe', 'F', amino_acid_change)
    amino_acid_change = re.sub('Tyr', 'Y', amino_acid_change)
    amino_acid_change = re.sub('Trp', 'W', amino_acid_change)
    amino_acid_change = re.sub('His', 'H', amino_acid_change)
    amino_acid_change = re.sub('Lys', 'K', amino_acid_change)
    amino_acid_change = re.sub('Arg', 'R', amino_acid_change)
    amino_acid_change = re.sub('Gln', 'Q', amino_acid_change)
    amino_acid_change = re.sub('Asp', 'D', amino_acid_change)
    amino_acid_change = re.sub('Ser', 'S', amino_acid_change)
    amino_acid_change = re.sub('Thr', 'T', amino_acid_change)
    amino_acid_change = re.sub('Asn', 'N', amino_acid_change)
    amino_acid_change = re.sub('Ile', 'I', amino_acid_change)
    amino_acid_change = re.sub('Glu', 'E', amino_acid_change)
    amino_acid_change = re.sub('Ala', 'A', amino_acid_change)
    amino_acid_change = re.sub('p\\.', '', amino_acid_change)
    return amino_acid_change


# Generate reference dictionary
# This function is for joblib parallel
def generate_reference_dictionary(line):
    reference_dictionary = {}
    values = line.strip("\n").split("\t")
    if values[0] not in reference_dictionary.keys():
        reference_dictionary[values[0]] = {
            "Sample": values[0],
            "Maturity_Group": str(values[1]) if (
                    isinstance(values[1], int) or
                    str(values[1]).isdigit() or
                    isinstance(values[1], float) or
                    str(values[1]).replace('.', '', 1).isdigit()) else "",
            "Country": values[2],
            "State": values[3],
            "Improvement_Status": values[4],
            "Classification": values[5]
        }
    else:
        reference_dictionary[values[0]]["Sample"] = values[0]
        reference_dictionary[values[0]]["Maturity_Group"] = int(values[1])
        reference_dictionary[values[0]]["Country"] = values[2]
        reference_dictionary[values[0]]["State"] = values[3]
        reference_dictionary[values[0]]["Improvement_Status"] = values[4]
        reference_dictionary[values[0]]["Classification"] = values[5]
    return reference_dictionary


# Generate gff dictionary
# This function is for joblib parallel
def generate_gff_dictionary(line, gff_category, gff_key):
    gff_dictionary = {}
    if (line.startswith("#")) or (re.search(gff_category, line) is None) or \
            (re.search(gff_key, line) is None) or (re.search("scaffold", line) is not None):
        return None
    else:
        values = line.strip("\n").split("\t")
        chromosome = values[0]
        category = values[2]
        start = values[3]
        stop = values[4]
        key = ""
        key_array = values[8].strip().split(";")
        if category != gff_category:
            return None
        else:
            for i in range(len(key_array)):
                if str(key_array[i]).startswith(gff_key):
                    key = re.sub("(.*=)|(.*:)|(.*-)", "", re.sub(gff_key, "", str(key_array[i])))
                    break
            if (key != "") and (key not in gff_dictionary.keys()):
                gff_dictionary[key.upper()] = {
                    "Chromosome": chromosome,
                    "Category": category,
                    "Start": start,
                    "Stop": stop,
                    "Gene": key
                }
            elif (key != "") and (key in gff_dictionary.keys()):
                gff_dictionary[key.upper()]["Chromosome"] = chromosome
                gff_dictionary[key.upper()]["Category"] = category
                gff_dictionary[key.upper()]["Start"] = start
                gff_dictionary[key.upper()]["Stop"] = stop
                gff_dictionary[key.upper()]["Gene"] = key
            return gff_dictionary


# Generate allele catalog
def generate_allele_catalog(header, line, reference_dictionary, gff_dictionary, output_file_path, lock):

    # Create an output array
    output_array = []

    # Parse header and line to get variant
    accessions = str(header).strip().split("\t")[9:]
    line_array = str(line).strip().split("\t")

    chromosome = line_array[0]
    position = line_array[1]

    reference_allele = line_array[3]
    annotated_reference_allele = reference_allele + "|Ref"

    alternate_alleles = line_array[4].split(",")

    genotypes = [re.split('/|\\|', genotype)[0] for genotype in line_array[9:]]

    info_dict = {}
    info = line_array[7].split(";")
    functional_effect_annotation_string = ""
    for i in range(len(info)):
        if info[i].startswith("EFF="):
            functional_effect_annotation_string = re.sub("EFF=", "", info[i])
            functional_effect_annotation_string = functional_effect_annotation_string.strip()
            break
        if info[i].startswith("ANN="):
            functional_effect_annotation_string = re.sub("ANN=", "", info[i])
            functional_effect_annotation_string = functional_effect_annotation_string.strip()
            break
    if functional_effect_annotation_string != "":
        functional_effect_annotation_string_array = functional_effect_annotation_string.split(",")
        for i in range(len(functional_effect_annotation_string_array)):
            single_functional_effect_annotation = re.split('\\|', str(functional_effect_annotation_string_array[i]))
            # Check if it is a primary transcript.
            # We only consider primary transcript.
            pattern = re.sub("\\\\", "\\\\\\\\", single_functional_effect_annotation[3])
            pattern = re.sub("\\+", "\\+", pattern)
            pattern = re.sub("\\*", "\\*", pattern)
            pattern = re.sub("\\?", "\\?", pattern)
            pattern = re.sub("\\(", "\\(", pattern)
            pattern = re.sub("\\)", "\\)", pattern)
            pattern = re.sub("\\[", "\\[", pattern)
            if str(single_functional_effect_annotation[6]).find(str(pattern+".1")) != -1:
                allele = str(single_functional_effect_annotation[0]).strip()
                gene = str(single_functional_effect_annotation[3]).strip()
                functional_effect = str(single_functional_effect_annotation[1]).strip()
                amino_acid_change = str(single_functional_effect_annotation[10]).strip()
                if any([str(functional_effect).find(effect) != -1 for effect in EFFECTS]):
                    if amino_acid_change != "":
                        amino_acid_change = clean_amino_acid_change_string(amino_acid_change)
                    if gene.upper() not in info_dict.keys():
                        info_dict[gene.upper()] = {}
                        if allele not in info_dict[gene.upper()].keys():
                            info_dict[gene.upper()][allele] = {
                                'Functional_Effects': functional_effect,
                                'Amino_Acid_Changes': amino_acid_change
                            }
                        else:
                            if info_dict[gene.upper()][allele]['Functional_Effects'] == "":
                                info_dict[gene.upper()][allele]['Functional_Effects'] = functional_effect
                            else:
                                info_dict[gene.upper()][allele]['Functional_Effects'] = info_dict[gene.upper()][allele]['Functional_Effects']+'&'+functional_effect
                            if info_dict[gene.upper()][allele]['Amino_Acid_Changes'] == "":
                                info_dict[gene.upper()][allele]['Amino_Acid_Changes'] = amino_acid_change
                            else:
                                info_dict[gene.upper()][allele]['Amino_Acid_Changes'] = info_dict[gene.upper()][allele]['Amino_Acid_Changes']+'&'+amino_acid_change
                    else:
                        if allele not in info_dict[gene.upper()].keys():
                            info_dict[gene.upper()][allele] = {
                                'Functional_Effects': functional_effect,
                                'Amino_Acid_Changes': amino_acid_change
                            }
                        else:
                            if info_dict[gene.upper()][allele]['Functional_Effects'] == "":
                                info_dict[gene.upper()][allele]['Functional_Effects'] = functional_effect
                            else:
                                info_dict[gene.upper()][allele]['Functional_Effects'] = info_dict[gene.upper()][allele]['Functional_Effects']+'&'+functional_effect
                            if info_dict[gene.upper()][allele]['Amino_Acid_Changes'] == "":
                                info_dict[gene.upper()][allele]['Amino_Acid_Changes'] = amino_acid_change
                            else:
                                info_dict[gene.upper()][allele]['Amino_Acid_Changes'] = info_dict[gene.upper()][allele]['Amino_Acid_Changes']+'&'+amino_acid_change
    
    genes = []
    for key in gff_dictionary.keys():
        if (gff_dictionary[key]["Chromosome"] == chromosome) and (float(gff_dictionary[key]["Start"].strip()) <= float(position.strip()) <= float(gff_dictionary[key]["Stop"].strip())):
            gene = gff_dictionary[key]["Gene"]
            if gene is not None:
                if gene != "":
                    genes.append(gene)

    # Collect all genotype information and write to the output file
    if ((len(genes) > 0) and (len(accessions) > 0) and (len(genotypes) > 0) and (len(accessions) == len(genotypes))):
        for i in range(len(genes)):
            for j in range(len(accessions)):
                try:
                    genotype = ""
                    genotype_with_description = ""
                    if genotypes[j] is not None:
                        if genotypes[j] != "":
                            genotype_index = int(genotypes[j])
                            if genotype_index == 0:
                                genotype = reference_allele
                                genotype_with_description = annotated_reference_allele
                            elif genotype_index > 0:
                                genotype = alternate_alleles[genotype_index-1]
                                flags = [str(info_key).find(genes[i].upper()) != -1 for info_key in info_dict.keys()]
                                flags = [flag_index for flag_index, flag in enumerate(flags) if flag]
                                if (len(flags) == 1) and (list(info_dict.keys())[flags[0]] in info_dict.keys()):
                                    if genotype in info_dict[genes[i].upper()].keys():
                                        if info_dict[genes[i].upper()][genotype]['Functional_Effects'] != "":
                                            genotype_with_description = genotype+"|"+info_dict[genes[i].upper()][genotype]['Functional_Effects']
                                            if info_dict[genes[i].upper()][genotype]['Amino_Acid_Changes'] != "":
                                                genotype_with_description = genotype_with_description+"|"+info_dict[genes[i].upper()][genotype]['Amino_Acid_Changes']
                                        else:
                                            genotype_with_description = genotype+"|Alt"
                                    else:
                                        genotype_with_description = genotype+"|Alt"
                                else:
                                    genotype_with_description = genotype+"|Alt"

                    # Check accession and append result to output array
                    if (accessions[j] in reference_dictionary.keys()) and (genes[i] != "") and (genotype != "") and (genotype_with_description != "") and (genotype != "<INS>") and (genotype != "<DEL>"):
                        output_array.append(
                            reference_dictionary[accessions[j]]["Classification"] + "\t" +
                            reference_dictionary[accessions[j]]["Improvement_Status"] + "\t" +
                            reference_dictionary[accessions[j]]["Maturity_Group"] + "\t" +
                            reference_dictionary[accessions[j]]["Country"] + "\t" +
                            reference_dictionary[accessions[j]]["State"] + "\t" +
                            accessions[j] + "\t" +
                            chromosome + "\t" +
                            genes[i] + "\t" +
                            position + "\t" +
                            genotype + "\t" +
                            genotype_with_description + "\n"
                        )
                    elif (accessions[j] not in reference_dictionary.keys()) and (genes[i] != "") and (genotype != "") and (genotype_with_description != "") and (genotype != "<INS>") and (genotype != "<DEL>"):
                        output_array.append(
                            "" + "\t" +
                            "" + "\t" +
                            "" + "\t" +
                            "" + "\t" +
                            "" + "\t" +
                            accessions[j] + "\t" +
                            chromosome + "\t" +
                            genes[i] + "\t" +
                            position + "\t" +
                            genotype + "\t" +
                            genotype_with_description + "\n"
                        )
                except Exception as e:
                    print(e)

    # Write output array to output file
    if len(output_array) > 0:
        lock.acquire()
        try:
            with open(output_file_path, "a") as writer:
                writer.write(''.join(output_array))
        except Exception as e:
            print(e)
        lock.release()


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_path = args.input_file
    reference_file_path = args.reference_file
    gff_file_path = args.gff_file
    output_file_path = args.output_file

    n_jobs = args.threads
    gff_category = args.gff_category
    gff_key = args.gff_key

    #######################################################################
    # Create a threading lock
    #######################################################################
    lock = threading.Lock()

    #######################################################################
    # Check if output parent folder exists
    # If not, create the output parent folder
    #######################################################################
    if not output_file_path.parent.exists():
        try:
            output_file_path.parent.mkdir(parents=True)
        except FileNotFoundError as e:
            pass
        except FileExistsError as e:
            pass
        except Exception as e:
            pass
        if not output_file_path.parent.exists():
            sys.exit(1)

    #######################################################################
    # Read reference file then make a reference dictionary
    #######################################################################
    reference_dictionary = {}
    if reference_file_path is not None:
        with open(reference_file_path, "r") as reader:
            header = reader.readline()
            reference_column_names = header.strip().split("\t")
            # The results is a dictionary list, unwind the dictionary list and create dictionary with all samples.
            results = Parallel(n_jobs=n_jobs)(
                delayed(generate_reference_dictionary)(line) for line in reader
            )
        for i in range(len(results)):
            for key in results[i].keys():
                if key not in reference_dictionary.keys():
                    reference_dictionary[key] = results[i][key]

    #######################################################################
    # Read gff file then make a gff dictionary
    #######################################################################
    with open(gff_file_path, "r") as reader:
        results = Parallel(n_jobs=n_jobs)(
            delayed(generate_gff_dictionary)(line.strip(), gff_category, gff_key) for line in reader
        )
    gff_dictionary = {}
    for i in range(len(results)):
        if results[i] is not None:
            for key in results[i].keys():
                if key not in gff_dictionary.keys():
                    gff_dictionary[key] = results[i][key]

    #######################################################################
    # Write header to the output file
    #######################################################################
    writer = open(output_file_path, "w")
    writer.write(
        "Classification\tImprovement_Status\tMaturity_Group\tCountry\tState\t" +
        "Accession\tChromosome\tGene\tPosition\tGenotype\tGenotype_with_Description\n"
    )
    writer.close()

    #######################################################################
    # Read input file then generate allele catalog
    #######################################################################
    if str(input_file_path).endswith('gz'):
        with gzip.open(input_file_path, 'rt') as reader:
            header = ""
            while not header.strip().startswith("#CHROM"):
                header = reader.readline()
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(generate_allele_catalog)(header, line, reference_dictionary, gff_dictionary, output_file_path, lock)
                for line in reader
            )
    else:
        with open(input_file_path, "r") as reader:
            header = ""
            while not header.strip().startswith("#CHROM"):
                header = reader.readline()
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(generate_allele_catalog)(header, line, reference_dictionary, gff_dictionary, output_file_path, lock)
                for line in reader
            )


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='generate_Allele_Catalog', description='generate Allele Catalog')

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
    parser.add_argument('-g', '--gff_file', help='GFF file', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    parser.add_argument('-r', '--reference_file', help='Reference file', type=pathlib.Path, default=None)
    parser.add_argument('-t', '--threads', help='Number of threads', type=int, default=10)
    parser.add_argument('-c', '--gff_category', help='Gff category', type=str, default='gene')
    parser.add_argument('-k', '--gff_key', help='Gff key', type=str, default='Name')

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
