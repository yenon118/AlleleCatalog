# AlleleCatalog

<!-- badges: start -->
<!-- badges: end -->

The AlleleCatalog is a pipeline built for generating allele catalog using next-generation sequencing (NGS) whole-genome sequencing datasets.

## Requirements

In order to run the AlleleCatalog, users need to install software, programming languages, and packages in their computing systems.
The software, programming languages, and packages include: 

```
Python3>=3.7.0
Snakemake>=5.31.0
Pandas>=1.1.3
Beagle>=5.2
SnpEff>=4.3
``` 

## Installation

You can install the AlleleCatalog from [Github](https://github.com/yenon118/AlleleCatalog.git) with:

```
git clone https://github.com/yenon118/AlleleCatalog.git
```

## Usage

#### Write a configuration file in json format

Please save the file with .json extension.

```
{
  "project_name": "Test_Soy1066",
  "workflow_path": "/scratch/yenc/projects/AlleleCatalog",
  "input_files": [
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr01.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr02.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr03.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr04.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr05.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr06.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr07.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr08.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr09.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr10.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr11.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr12.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr13.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr14.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr15.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr16.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr17.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr18.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr19.vcf.gz",
    "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/unimpute_data/test_Soy1066_Chr20.vcf.gz"
  ],
  "chromosomes": [
    "Chr01",
    "Chr02",
    "Chr03",
    "Chr04",
    "Chr05",
    "Chr06",
    "Chr07",
    "Chr08",
    "Chr09",
    "Chr10",
    "Chr11",
    "Chr12",
    "Chr13",
    "Chr14",
    "Chr15",
    "Chr16",
    "Chr17",
    "Chr18",
    "Chr19",
    "Chr20"
  ],
  "metadata_file": "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/Soy1066_allele_line_info.txt",
  "genome_version": "Wm82.a2.v1",
  "gff_file": "/scratch/yenc/projects/AlleleCatalog/data/test_Soy1066/Wm82.a2.v1.genes.tiny.gff",
  "gff_category": "gene",
  "gff_key": "Name",
  "output_folder": "/scratch/yenc/projects/AlleleCatalog/output/Test_Soy1066/",
  "memory": 100,
  "threads": 4
}
```

#### Run workflow with the Snakemake workflow management system

```
snakemake -j NUMBER_OF_JOBS --configfile CONFIGURATION_FILE --snakefile SNAKEMAKE_FILE

Mandatory Positional Argumants:
    NUMBER_OF_JOBS                          - the number of jobs
    CONFIGURATION_FILE                      - a configuration file
    SNAKEMAKE_FILE                          - the AlleleCatalog.smk file that sit inside this repository 
```

## Examples

These are a few basic examples which show you how to use the AlleleCatalog:

```
cd /path/to/AlleleCatalog

snakemake -pj 3 --configfile inputs.json --snakefile AlleleCatalog.smk
```

```
cd /path/to/AlleleCatalog

snakemake --cluster "sbatch --account=xulab --cpus-per-task=10 --time=0-02:00 \
--partition=Lewis,BioCompute,hpc5,General --mem=64G" \
--jobs 25 --latency-wait 30 \
--configfile /storage/htc/joshilab/yenc/projects/AlleleCatalog/inputs.json \
--snakefile /storage/htc/joshilab/yenc/projects/AlleleCatalog/AlleleCatalog.smk
```