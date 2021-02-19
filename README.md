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
  "project_name": "Test",
  "workflow_path": "/storage/htc/joshilab/yenc/projects/AlleleCatalog",
  "input_files": [
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr01.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr02.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr03.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr04.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr05.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr06.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr07.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr08.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr09.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr10.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr11.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr12.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr13.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr14.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr15.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr16.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr17.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr18.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr19.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Nebraska/Nebraska_Chr20.vcf"
  ],
  "reference_files": [
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr01.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr02.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr03.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr04.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr05.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr06.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr07.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr08.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr09.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr10.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr11.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr12.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr13.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr14.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr15.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr16.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr17.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr18.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr19.vcf",
    "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/Soy775/Soy775_Chr20.vcf"
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
  "ancestors": [
    "PI_634883",
    "PI_548348",
    "PI_594301",
    "PI_578457A",
    "PI_548190",
    "ZJ-Y191",
    "PI_603318",
    "PI_464929A",
    "Jin_Shan_Cha_Zhu_Shi_Dou",
    "Fen_Dou_No_65",
    "Bai_Mao_Dou",
    "PI_567364",
    "PI_393551",
    "PI_593983",
    "Hu_Pi_Dou"
  ],
  "reference_file": "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/completed_allele_line_info.txt",
  "genome_version": "Wm82.a2.v1",
  "gff_file": "/storage/htc/joshilab/yenc/projects/AlleleCatalog/data/just_genes_v2.gff",
  "gff_category": "gene",
  "gff_key": "Name",
  "output_folder": "/storage/htc/joshilab/yenc/projects/AlleleCatalog/output/",
  "memory": 100,
  "threads": 10
}
```

#### Run workflow with the Snakemake workflow management system

```
snakemake -j NUMBER_OF_JOBS --configfile CONFIGURATION_FILE --snakefile SNAKEMAKE_FILE

Mandatory Positional Argumants:
    NUMBER_OF_JOBS                          - the number of jobs
    CONFIGURATION_FILE                      - a configuration file
    SNAKEMAKE_FILE                          - the snakyVC.smk file that sit inside this repository 
```

## Examples

These are a few basic examples which show you how to use the AlleleCatalog:

```
cd /path/to/AlleleCatalog

snakemake -pj 10 --configfile inputs.json --snakefile AlleleCatalog.smk
```

```
cd /path/to/AlleleCatalog

snakemake --cluster "sbatch --account=xulab --cpus-per-task=10 --time=0-02:00 \
--partition=Lewis,BioCompute,hpc5,General --mem=64G" \
--jobs 25 --latency-wait 30 \
--configfile /storage/htc/joshilab/yenc/projects/AlleleCatalog/inputs.json \
--snakefile /storage/htc/joshilab/yenc/projects/AlleleCatalog/AlleleCatalog.smk
```