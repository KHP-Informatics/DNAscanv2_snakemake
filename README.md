![alt text](https://github.com/hevmarriott/DNAscanv2/blob/master/DNAscan_logo.001.jpeg)

# DNAscanv2 (Snakemake Version)
```diff
+TO BE NOTED: We are always working to improve DNAscan so any bug reports, suggestions and general feedback would be highly welcome. 
```
## Table of Contents
1. [Introduction](#introduction)
2. [Documentation](#documentation)
    * [Obtaining](#obtaining)
    * [Dependencies](#dependencies)
    * [ALSgeneScanner](#alsgenescanner)
    * [Output](#output)
    * [Usage](#usage)
3. [Core Contributors](#core-contributors)
4. [Contributing](#contributing)
5. [Licence](#licence)

## Introduction 
DNAscan2 is a fast and efficient bioinformatics pipeline that allows for the analysis of DNA Next Generation Sequencing data. It is capable of screening for single nucleotide variants, small indels, structural variants, transposable elements, known and non-reference repeat expansions, short tandem repeats and viral (or any other organism's) genetic material. The identified variants can then be annotated using multiple databases including refGene, Clinvar, ExAC, dbSNP, CADD and gnomAD/1000 genomes variant frequency estimation. Information about the annotated variants can be obtained via the generation of several results reports. Sequencing and alignment reports can also be generated. More information about the architecture of DNAscan2 with all analysis options can be found on the main [DNAscan2 repository](https://github.com/hevmarriott/DNAscanv2/blob/master/README.md#introduction).

However, the command line implementation of DNAscan2 does not allow for parallel/multi-sample processing, which for , can become a bottleneck for 

## Documentation
The general documentation for DNAscan2 is available by clicking on the following [link](https://github.com/hevmarriott/DNAscanv2/#documentation). 

NOTE: 

### Obtaining

### Usage

### ALSgeneScanner 

### Output

### Dependencies

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

