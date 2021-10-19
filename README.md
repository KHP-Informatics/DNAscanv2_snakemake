![alt text](https://github.com/hevmarriott/DNAscanv2/blob/master/DNAscan_logo.001.jpeg)

# DNAscanv2 (Snakemake Version)
```diff
+TO BE NOTED: We are always working to improve DNAscan so any bug reports, suggestions and general feedback would be highly welcome. 
```
## Table of Contents
1. [Introduction](#introduction)
2. [Documentation](#documentation)
    * [Obtaining](#obtaining)
    * [Installation](#installation)
    * [Dependencies](#dependencies)
    * [ALSgeneScanner](#alsgenescanner)
    * [Output](#output)
    * [Usage](#usage)
3. [Core Contributors](#core-contributors)
4. [Contributing](#contributing)
5. [Licence](#licence)

## Introduction 
DNAscan2 is a fast and efficient bioinformatics pipeline that allows for the analysis of DNA Next Generation Sequencing data. It is capable of screening for single nucleotide variants, small indels, structural variants, transposable elements, known and non-reference repeat expansions, short tandem repeats and viral (or any other organism's) genetic material. The identified variants can then be annotated using multiple databases including refGene, Clinvar, ExAC, dbSNP, CADD and gnomAD/1000 genomes variant frequency estimation. Information about the annotated variants can be obtained via the generation of several results reports. Sequencing and alignment reports can also be generated. More information about the architecture of DNAscan2 with all analysis options can be found on the main [DNAscan2 repository](https://github.com/hevmarriott/DNAscanv2/blob/master/README.md#introduction).

As the command line implementation of DNAscan2 does not allow for parallel/multi-sample processing and is not optimised for high performance computing (HPC) system and cluster execution, DNAscan2 has been redesigned as a snakemake workflow which can perform each step in a parallel manner on tens to hundreds of samples. This is beneficial for performance of larger scale genomic studies or for laboratories who have multiple sampels that they want to process. 

## Documentation
The general documentation for DNAscan2 is available by clicking on the following [link](https://github.com/hevmarriott/DNAscanv2/#documentation). 

### Obtaining
Please download this repository to your system:

```bash
git clone https://github.com/hevmarriott/DNAscan2_snakemake
```

### Installation
Before running the workflow, Snakemake and its associated packages have to be installed on your system. The best way to do this is via the Miniconda3 package. If you do not have this installed, run the following commands on your system to download this and the Conda package manager:

```bash 
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh -b -p /path/to/Miniconda3/installation/directory
```
In your $HOME/.bashrc file, add the following line:

```bash 
export PATH=/path/to/Miniconda3/installation/directory/Miniconda3/bin:$PATH 
```

Once Miniconda3 is installed on your system, create a standalone Conda environment containing Snakemake. The name (-n flag) of the environment should be relevant to its purpose:

```bash
conda create --name Snakemake snakemake
```

Then activate the environment:

```bash
conda activate Snakemake
```

NOTE: This environment should always be activated prior to running DNAscan2.

To submit jobs on the selected cluster, there are two methods (cluster execution profile and command line cluster configuration). The earlier one is currently recommended by the developers of Snakemake as cluster configuration is deprecated but can still be used. To create a job execution profile, install cookiecutter into the activated snakemake environment, which ensures that the job profile can be created:

```bash
conda install -c conda-forge cookiecutter
```

Navigate to the [Snakemake profile templates page](https://github.com/Snakemake-Profiles), choose your HPC/execution system and follow the instructions for installation into the DNAscanv2_snakemake directory. The results will be a named profile directory with the structure *HPC_system.account_name*, which you will need to specify when running DNAscan2. 

Alternatively, you can specify a cluster_config.json file on the command line. Instructions on how to do this is available [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration-deprecated).

To demonstrate the differences in execution between the cluster execution methods, consider a scenario where SNV/indel calling and filtering has to be performed on 10 samples on a SLURM cluster:

Method 1: Snakemake job execution profile
``` bash
snakemake -j 10 "output/results/{sample_1,sample_2...}/{sample_1,sample_2...}_sorted_filtered.vcf.gz" --profile slurm.<ACCOUNT_NAME>
``` 

Method 2: cluster_config.json
```bash 
snakemake -j 10 "output/results/{sample_1,sample_2...}/{sample_1,sample_2...}_sorted_filtered.vcf.gz" --cluster-config cluster_config.json --cluster "sbatch -p <PARTITION_NAME> <OTHER_SLURM_OPTIONS>" 
```

The full usage instructions are explained below. 

### Usage
The basic and advanced options of this workflow are the same as in the DNAscan2 command line implementation, which you can read in more detail [here](https://github.com/hevmarriott/DNAscanv2#usage). 

All of these options, along with paths to the reference genome/indexes and other tools unavailable via Conda can be specified from the config.yaml file. When defining directories, the paths have to end in '/'. Furthermore, all options that require boolean input (in the 'DNAscan2 command line' and 'DNAscan2 options' sections) must be populated with either 'true' or 'false', and the custom options for HISAT2, BWA, MELT and AnnotSV have to be kept as 'None' if you do not wish to provide the workflow with these values. 



### ALSgeneScanner 



### Output
The output tree of the DNAscan2 snakemake workflow is as follows:



### Dependencies
All of the necessary binary dependencies for each step in the DNAscan2 workflow are installed and deployed via the 'conda:' directive, using environment files located in the envs/ directory. 

There are some dependencies that can only be downloaded using direct installation, listed below with versions and download links:
* Strelka (SNV/indel variants) = [v2.9.10](https://github.com/Illumina/strelka/releases/tag/v2.9.10)
* Manta (structural variants) = [v1.6.0](https://github.com/Illumina/manta/releases/tag/v1.6.0)
* MELT (transposable elements) = [v2.2.2](https://melt.igs.umaryland.edu/downloads.php)
* ExpansionHunter Denovo (genome-wide short tandem repeats) = [v0.9.0](https://github.com/Illumina/ExpansionHunterDenovo/releases/tag/v0.9.0)
* ANNOVAR (SNV/indel + expansion/short tandem repeat annotation) = [latest](https://www.openbioinformatics.org/annovar/annovar_download_form.php)
* AnnotSV (structural variant/transposable element annotation) = [v3.0.8](https://github.com/lgmgeo/AnnotSV/releases/tag/v3.0.8)
* knotAnnotSV (annotated structural variant/transposable element report) = [v1.0.0](https://github.com/mobidic/knotAnnotSV/releases/tag/v1.0.0)

NOTE: if wanting to perform annotation with ANNOVAR and MELT, a manual registration step is required prior to download. 

Furthermore, this workflow requires that you provide the reference genome and associated HISAT2 and BWA indexes in the config.yaml file. Further details on how to to download and index the reference genome and obtain microbe databases (for custom microbe analysis) is available at the main [DNAscan2 repository](https://github.com/hevmarriott/DNAscanv2#how-to-download-the-reference-genome). 

For SNV/indel, expansion and STR annotation, you need to manually install these Annovar databases - refGene,dbnsfp33a,clinvar_20210501,intervar_20180118,avsnp147, exac03,1000g_aug_2015 and gnomad211_genome, according to the following example for hg19:

```bash
cd <ANNOVAR_DIR>
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene /path/to/annovar/humandb/
```
This should be repeated for every database listed above. The above list with their operations are already defined in the config.yaml file - if you do not want to use all of these listed databases, just delete their entry from the ANNOVAR_OPERATIONS and ANNOVAR_PROTOCOLS fields. 

## Core Contributors
- [Heather Marriott](heather.marriott@kcl.ac.uk), UK

For a full list of contributors see [LINK](https://github.com/hevmarriott/DNAscan2/CONTRIBUTORS.md)

## Contributing

Here’s how we suggest you go about proposing a change to this project:

1. [Fork this project][fork] to your account.
2. [Create a branch][branch] for the change you intend to make.
3. Make your changes to your fork.
4. [Send a pull request][pr] from your fork’s branch to our `master` branch.

Using the web-based interface to make changes is fine too, and will help you
by automatically forking the project and prompting to send a pull request too.

[fork]: https://help.github.com/articles/fork-a-repo/
[branch]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository
[pr]: https://help.github.com/articles/using-pull-requests/

## Licence
- [MIT](./LICENSE.txt)

