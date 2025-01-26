# QuasiFlow2

[![](https://img.shields.io/badge/nextflow-24.10.3-yellowgreen)](https://www.nextflow.io)
[![](https://img.shields.io/badge/uses-conda-yellowgreen)](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
[![](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## Acknowledgments

This project is based on work created by Alfred Ssekagiri https://github.com/AlfredUg.  
You can find the original repository here: https://github.com/AlfredUg/QuasiFlow.

## Introduction

QuasiFlow2 is a nextflow pipeline for reproducible analysis of NGS-based HIVDR testing data. It is meant to be used in conda environment. The pipeline takes raw sequence reads in FASTQ format as input, performs quality control, maps reads to a reference genome, and performs variant calling. It is based on the DLS1 pipeline of the original QuasiFlow https://github.com/AlfredUg/QuasiFlow (it works on nextflow 22.12.0 or lower, which supports DSL1). Quasiflow2 is made in the DSL2 version of the nextflow. It is tested with nextflow version 24.10.3, but it should work with newer ones, as well.  

## Installation

QuasiFlow2 requires **nextflow** (version 24.10.3 or higher), and it works under **conda** environment. It was tested within **miniforge**, but it should also work in other conda environments. Instruction for miniforge installation can be found here: https://github.com/conda-forge/miniforge. 

After the conda installation, the pipeline repository must be cloned into a desired directory. It is recommended to download it somewhere the most easily accessible, for example, the `$HOME` directory:

```bash
cd ~
git clone https://github.com/nedjoni/QuasiFlow2
```
It's necessary to create a **quasiflow2** working environment. To create it and install the required packages, you may use the command:

```bash
conda env create -f ~/QuasiFlow2/environment.yml -y
```

Activate the **nextflow** environment and confirm that installation was successful by printing out the help message:

```bash
conda activate quasiflow2
nextflow run ~/QuasiFlow2 --help
```

## Usage

The pipeline takes as input paired-end Illumina data in `fastq` or `fastq.gz` format. Data for the exercise can be downloaded and tested from the European Nucleotide Archive (ENA) using the `wget` command. This is paired-end data from a single sample of bioProject PRJDB3502.

```bash
wget -P ~/QuasiFlow2/fastq ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR030/DRR030218/DRR030218_1.fastq.gz 
wget -P ~/QuasiFlow2/fastq ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR030/DRR030218/DRR030218_2.fastq.gz 
```
For cleaner directory organization it is recommended to create a directory where all pipeline working directories will be created. Example:

```bash
mkdir ~/QuasiTest
cd ~/QuasiTest
```
Run QuasiFlow2 on a test dataset with default parameters under the `conda` profile. If you already installed the pipeline using the procedure explained in the installation, activate the **nextflow** environment run command:

```bash
nextflow run ~/QuasiFlow2
```

If you have fastq or fastq.gz reads in a specific location, point to the path as follows:
```bash
nextflow run ~/QuasiFlow2 --reads "path_to_the_reads_directory/*_{R1,R2,1,2}*.fastq{,.gz}" 
```
You have to provide adequate naming for the Illumina reads. Usually, they are named in the format "SampleName_S1_L001_R1_001.fastq.gz" or "SampleName_1.fastq". Nay other naming format would require adaptation of the reads path. Pipeline accepts both fastq/fastq.gz files. More than one pair of reads is supported.

If you want results in a specific location, point to the path as follows:
```bash
nextflow run path/to/QuasiFlow --outdir "path_to_the_results_directory" 
```

#### Pipeline 

**Outputs Quality control**

* `raw_reads_multiqc_report.html`: Aggregated quality control data and visualizations - one file for the entire dataset

**Variants and drug resistance outputs**

* `consensus*.fasta`: `FASTA` files of consensus sequences - one per sample
* `consensus*.json`: `JSON` files of detailed HIV drug resistance analysis - one per sample
* `dr_report*.csv`: `CSV` files of drug resistance mutations at different mutational frequencies - one per sample
* `filtered*.fastq`: `FASTQ` files of drug resistance mutations at different mutational frequencies - one per sample
* `mutation_report*.aavf`: `AAVF` files of amino acid variant calls - one per sample
* `hivdr*.html`: `HTML` Final drug resistance report - one per sample

**Pipeline information output**

* `QuasiFlow_DAG.html`: Graphical representation of the pipeline's processes/operators and channels between them.
* `QuasiFlow_report.html`: Overall start and completion time, CPU and memory usage.
* `QuasiFlow_timeline.html`: Timeline for all the processes executed in the pipeline.


## Parameters

### Input parameters

(Optional parameters)

* `--reads`: Path to input data (must be surrounded with quotes, default is the folder "fastq")

### HyDRA parameters

(Optional parameters)

* `--reporting_threshold`: Minimum mutation frequency percent to report.

* `--consensus_pct`: The minimum percentage of a base needs to be incorporated into the consensus sequence.

* `--length_cutoff`: Reads that fall short of the specified length will be filtered out.

* `--score_cutoff`: Reads that have a median or mean quality score (depending on the score type specified) less than the score cutoff value will be filtered out.

* `--min_variant_qual`: Minimum quality for the variant to be considered later on in the pipeline.

* `--min_dp`: Minimum required read depth for the variant to be considered later on in the pipeline.

* `--min_ac`: The minimum required allele count for the variant to be considered later on in the pipeline.

* `--min_freq`: The minimum required frequency for a mutation to be considered in the drug resistance report.


### Output parameters

(Optional parameters)

* `--outdir`: Directory path where results will be saved (default - results).

### Other arguments

(Optional parameters)

* `--overwrite`: Set to true to overwrite previous reports (default - false).

## Dependencies.

Below is the list of tools that are used in the QuasiFlow pipeline. These tools are readily available and may be installed using `conda` via `bioconda` channel.

+ [fastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc)
+ [MultiQC](https://multiqc.info/)
+ [Trim-galore](https://github.com/FelixKrueger/TrimGalore)
+ [Quasitools](https://phac-nml.github.io/quasitools/)

## Conclusion

This pipeline is made solely for more platform/universal use. It is made to be used with the current, DSL2 version of the nextflow. It is also simplified, so only a few steps are necessary for its installation and use.
It is made and tested under the Ubuntu WSL subsystem for Windows, but it should work under regular Ubuntu or other similar OS. 
In the original QuasiFlow (https://github.com/AlfredUg/QuasiFlow) there were additional reports available as well. I have tried to add it here as well, but it didn't work as intended. The problem lies with current databases that update periodically, so even if reports are created by the internal tools, they would not be entirely up to date. So, I've decided to exclude that option from the code, and you can use online tools for that purpose:

+ [HIVdb Program: Sequence Analysis](https://hivdb.stanford.edu/hivdb/by-sequences/) Uses fasta files. Very nice reports and HIVdb is actively maintained
+ [HIVdb Program: Sequence Reads (NGS) Analysis](https://hivdb.stanford.edu/hivdb/by-reads/) Uses aavf files. Variant data reports.
+ [HIV-GRADE](https://www.hiv-grade.de/grade_new/) Uses fasta files. Good report for drug resistance mutations. It has an option to include results from ANRS and HIVDB as well.

## Troubleshooting

Report any issue at https://github.com/nedjoni/QuasiFlow2/issues. I may be slower to respond because this project won't be maintained actively.

## License

QuasiFlow2 is licensed under GNU GPL v3.


## Special thanks

Special thanks go to Katja W. for introducing me to the original QuasiFlow pipeline.