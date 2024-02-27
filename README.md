# AlleleCatalog

<!-- badges: start -->
<!-- badges: end -->

The AlleleCatalog pipeline is built for generating Allele Catalog datasets using next-generation sequencing (NGS) based genetic data and metadata.

## Requirements

In order to run the AlleleCatalog pipeline, users need to install Miniconda and prepare the Miniconda environment in their computing systems.

The required software, programming languages, and packages include:

```
bwa>=0.7.17
gatk4>=4.4.0.0
samtools>=1.6
htslib>=1.3
python>=3.12
snakemake>=8.4
scipy>=1.12
numpy>=1.26
pandas>=2.2
Beagle>=5.2
SnpEff>=4.3
```

Miniconda can be downloaded from [https://docs.anaconda.com/free/miniconda/](https://docs.anaconda.com/free/miniconda/).

For example, if users plan to install Miniconda3 Linux 64-bit, the wget tool can be used to download the Miniconda.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

To install Miniconda in a server or cluster, users can use the command below.

Please remember to replace the _<installation_shell_script>_ with the actual Miniconda installation shell script. In our case, it is **Miniconda3-latest-Linux-x86_64.sh**.

Please also remember to replace the _<desired_new_directory>_ with an actual directory absolute path.

```
chmod 777 -R <installation_shell_script>
./<installation_shell_script> -b -u -p <desired_new_directory>
rm -rf <installation_shell_script>
```

After installing Miniconda, initialization of Miniconda for bash shell can be done using the command below.

Please also remember to replace the _<desired_new_directory>_ with an actual directory absolute path.

```
<desired_new_directory>/bin/conda init bash
```

Installation of the Miniconda is required, and Miniconda environment needs to be activated every time before running the AlleleCatalog pipeline.

Write a Conda configuration file (.condarc) before creating a Conda environment:

```
nano ~/.condarc
```

Put the following text into the Conda configuration file (make sure you change _envs_dirs_ and _pkgs_dirs_) then save the file.

Please make sure not use tab in this yaml file, use 4 spaces instead.

Please make sure to replace _/new/path/to/_ with an actual directory absolute path.

```
envs_dirs:
    - /new/path/to/miniconda/envs
pkgs_dirs:
    - /new/path/to/miniconda/pkgs
channels:
    - conda-forge
    - bioconda
    - defaults
```

Create a Conda environment by specifying all required packages (option 1).

Please make sure to replace the _<conda_environment_name>_ with an environment name of your choice.

```
conda create -n <conda_environment_name> bioconda::gatk4 bioconda::samtools bioconda::bcftools bioconda::htslib \
bioconda::bedtools bioconda::bwa bioconda::snakemake bioconda::snakemake-executor-plugin-cluster-generic \
conda-forge::numpy conda-forge::pandas conda-forge::scipy
```

Create a Conda environment by using a yaml environment file (option 2).

Please make sure to replace the _<conda_environment_name>_ with an environment name of your choice.

```
conda create --name <conda_environment_name> --file AlleleCatalog-environment.yml
```

Create a Conda environment by using a explicit specification file (option 3).

Please make sure to replace the _<conda_environment_name>_ with an environment name of your choice.

```
conda create --name <conda_environment_name> --file AlleleCatalog-spec-file.txt
```

Activate Conda environment using conda activate command.

This step is required every time before running AlleleCatalog pipeline.

Please make sure to replace the _<conda_environment_name>_ with an environment name of your choice.

```
conda activate <conda_environment_name>
```

## Installation

You can install the AlleleCatalog from [Github](https://github.com/yenon118/AlleleCatalog.git) with:

```
git clone https://github.com/yenon118/AlleleCatalog.git
```

Create a "tools" folder within the AlleleCatalog pipeline clone.

```
cd AlleleCatalog
mkdir -p tools
cd tools
```

The Beagle imputation tool can be downloaded from [https://faculty.washington.edu/browning/beagle/beagle.html](https://faculty.washington.edu/browning/beagle/beagle.html) and placed inside the "tools" folder.

The SnpEff functional effect prediction tool can be downloaded from [https://pcingola.github.io/SnpEff/download/](https://pcingola.github.io/SnpEff/download/) and placed in the "tools" folder.

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

Below are some fundamental examples illustrating the usage of the AlleleCatalog pipeline.

Please adjust _/path/to/_ to an actual directory absolute path.

**Examples of running without an executor.**

```
cd /path/to/AlleleCatalog

snakemake -pj 3 --configfile inputs.json --snakefile AlleleCatalog.smk
```

**Examples of running with an executor.**

Snakemake version >= 8.0.0.

```
cd /path/to/AlleleCatalog

snakemake --executor cluster-generic \
--cluster-generic-submit-cmd "sbatch --account=xulab --time=0-02:00 \
--nodes=1 --ntasks=1 --cpus-per-task=3 \
--partition=Lewis,BioCompute,hpc5,General --mem=64G" \
--jobs 25 --latency-wait 60 \
--configfile lewis_slurm_inputs.json \
--snakefile AlleleCatalog.smk
```

Snakemake version < 8.0.0.

```
cd /path/to/AlleleCatalog

snakemake --cluster "sbatch --account=xulab --time=0-02:00 \
--nodes=1 --ntasks=1 --cpus-per-task=3 \
--partition=Lewis,BioCompute,hpc5,General --mem=64G" \
--jobs 25 --latency-wait 60 \
--configfile lewis_slurm_inputs.json \
--snakefile AlleleCatalog.smk
```

## Citation

Chan YO, Dietz N, Zeng S, Wang J, Flint-Garcia S, Salazar-Vidal MN, Škrabišová M, Bilyeu K, Joshi T: **The Allele Catalog Tool: a web-based interactive tool for allele discovery and analysis.** BMC Genomics 2023, 24(1):107.
