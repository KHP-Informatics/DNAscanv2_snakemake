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

Once Miniconda3 is installed on your system, create a Conda environment to download Snakemake to. The name (-n flag) of the environment should be relevant to its purpose:

```bash
conda create --name Snakemake snakemake
```

Then activate the environment:

```bash
conda activate Snakemake
```

NOTE: This environment should always be activated prior to running DNAscan2.

### Usage


### ALSgeneScanner 



### Output


### Dependencies
All of the necessary binary dependencies for each step in the DNAscan2 workflow are installed and deployed via the 'conda:' directive, using environment files located in the envs/ directory. 

There are some dependencies that work better if they are downloaded. 

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

