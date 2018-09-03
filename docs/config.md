
## Pipeline configuration

Various config files are provided with the pipeline. Some contain pipeline-specific parameters for use with Nextflow in order to increase the flexibility and portability of the pipeline, and others contain parameters required for specific processes.

| Path                                                                                                                   | Used by       | Description                                                                                                                                                                                                                                                                                                                                                             |
| -------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`nextflow.config`](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/nextflow.config)                           | `nextflow`                                                                                                             | Config file used for the definition of Nextflow profiles and other standard parameters. Dont need to edit this file.                                                                                                                                           |
|                                                                                                                        |                                                                                                                        |                                                                                                                                                                                                                                                                |
| [`conf/genomes.config`](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/genomes.config)                   | `nextflow`                                                                                                             | You can customise this file if you prefer to store the genome parameters required by the pipeline in a single file. Alternatively, you can provide each of the parameters at the command-line (See [Reference genome](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/genome.md) section). |
| [`conf/base.config `](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/base.config)                        | `nextflow`                                                                                                             | Config file defining compute resources used by each process. If applicable, you can tailor the parameters in this file to suit your compute system.                                                                                                            |
| [`conf/conda.config`](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/conda.config)                       | `nextflow`                                                                                                             | Config file specifying parameters required to build Conda environment. If applicable, you will need to change the `beforeScript` parameter in this file to ensure the relevant Conda distribution is available on your system before each process is executed. |
| [`conf/babs_modules.config`](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/babs_modules.config)         | `nextflow`                                                                                                             | Config file specifying environment modules that can be used at The Francis Crick Institute to run the pipeline. You can customise this file to use specific versions of software for each process in the pipeline or you can use it as a template to build your own for use outside The Francis Crick Institute.                                          |
| [`conf/fastq_screen.conf.txt`](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/fastq_screen.conf.txt)     | [`FastQ Screen`](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html#configuration) | If you want the pipeline to run `fastq_screen` this file will need to be customised to reflect the paths to pre-created bowtie2 indices for the contaminant screen.                                                                                            |
| [`conf/multiqc_config.yaml`](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/multiqc_config.yaml)         | [`MultiQC`](http://multiqc.info/docs/#configuring-multiqc)                                                             | Config file for generating customised MultiQC report for BABS-ATACSeqPE pipeline. You wont need to edit this file unless you want to change aspects of the QC that are reported by MultiQC.                                                               |
| [`conf/bamtools_filter_pe.json`](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/bamtools_filter_pe.json) | [`BamTools`](https://insidedna.me/tool_page_assets/pdf_manual/bamtools.pdf)                                            | Config file used by BamTools for read alignment filtering. You wont need to customise this file unless you want to amend specific aspects of the read filtering criteria.                                                                                      |

## Software configuration

Depending on where and how you would like to run the pipeline, Nextflow needs to be configured to use the [required software](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/config.md#software-requirements).

### Conda

This is the easiest, most hassle-free and portable way of running the pipeline. A custom [Conda environment file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/environment.yaml) is provided with the pipeline for those that wish to run the pipeline on another system without having to painstakingly install all of the software and associated dependencies. This will require an internet connection on the command-line, and installation of [Anaconda or Miniconda](https://conda.io/docs/user-guide/install/index.html). Nextflow will rather amazingly create a temporary Conda environment by downloading and installing all of the required software before execution of the pipeline. **NOTE: This could take up to 45 minutes.**  

A [Conda config file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/conda.config) can be found in the `conf/` directory. By default, the pipeline will be executed using the `slurm` job submission system. You will need an account to use the HPC cluster on CAMP if you want to run this at The Francis Crick Institute. If you prefer to run the pipeline locally in serial just replace `executor = 'slurm'` to `executor = 'local'` in the Conda config file. However, the pipeline may take a very long time to complete and this is only recommended for testing purposes. Before the submission of each Nextflow process the `conda` command will need to be available on the command-line to activate the Conda environment. If you are not running the pipeline at The Francis Crick Institute you will need to edit `beforeScript = 'module purge && ml Anaconda2/5.1.0'` in the config file to load/use Conda.

By default, Conda will download and compile the packages for the environment in the user `HOME` directory, however, this has a space limitation of 5GB on CAMP which wont be large enough. To overcome this, you can create a file called `.condarc` in your `HOME` directory that instructs Conda to install the environment and associated packages in your `working` directory. Just copy the text below, place it in the `.condarc` file, and change the paths to reflect those that are most appropriate for you.  

```bash
envs_dirs:
  - /camp/stp/babs/working/patelh/conda/envs/
pkgs_dirs:
  - /camp/stp/babs/working/patelh/conda/pkgs/
```

It is also possible to create a permanent copy of the Conda environment for recurrent use and Nextflow can be configured accordingly, however this wont be outlined here.

**When running the pipeline just specify `-profile conda` in order to use this configuration.**

### Environment modules

If you are running this pipeline at The Francis Crick Institute all of the required software can be loaded with the environment module system on CAMP. A [module config file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/conf/babs_modules.config) specifying the software modules to load has been been created for this purpose. If you want to change the versions of software just edit this file before running the pipeline, however, this may require some testing to make sure the Nextflow processes requiring the software can still be run with the same command-line parameters.  

By default, the pipeline will be executed using the `slurm` job submission system. You will need an account to use the HPC cluster on CAMP if you want to run this at The Francis Crick Institute. If you prefer to run the pipeline locally in serial just replace `executor = 'slurm'` to `executor = 'local'` in the module config file. However, the pipeline may take a very long time to complete and this is only recommended for testing purposes.

**When running the pipeline just specify **`-profile babs_modules`** in order to use this configuration.**

### Standard

If you really need to run the pipeline in serial you can specify `-profile standard` when running the pipeline. This is the most basic profile and will require all the software to be available at the command-line. The processes will be submitted in serial on the logged in environment.

## Software requirements

The software below is required to run the pipeline:

|                                                                                  |                                                                       |                                                                  |
|----------------------------------------------------------------------------------|-----------------------------------------------------------------------|------------------------------------------------------------------|
| [nextflow](https://www.nextflow.io/)                                             | [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html)       | [Kent_tools](http://hgdownload.soe.ucsc.edu/admin/exe/)          |
| [R](https://www.r-project.org/)                                                  | [BWA](https://sourceforge.net/projects/bio-bwa/files/)                | [MACS2](https://github.com/taoliu/MACS)                          |
| [Python](https://www.python.org/downloads/)                                      | [picard](https://broadinstitute.github.io/picard/)                    | [HOMER](http://homer.ucsd.edu/homer/download.html)               |
| [Java](https://java.com/en/download/)                                            | [SAMtools](https://sourceforge.net/projects/samtools/files/samtools/) | [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)        |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)             | [BamTools](https://github.com/pezmaster31/bamtools)                   | [Pysam](http://pysam.readthedocs.io/en/latest/installation.html) |
| [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) | [BEDTools](https://github.com/arq5x/bedtools2/)                       | [MultiQC](http://multiqc.info/)                                  |

### R libraries

The following R libraries may need to be installed if unavailable. You can test this be loading the `R` module (if required), typing `R` at the command prompt and attempting to load the packages below e.g. `> library(optparse)` and so on. The pipeline assumes the correct R library path is set in order find the installed packages. If not, you can set this in the `.Rprofile` file in the user home directory or add a line which extends the `R` [libPaths](https://stat.ethz.ch/R-manual/R-devel/library/base/html/libPaths.html) in the executable R scripts in the `bin/` directory.

|                                                                         |                                                                                 |
|-------------------------------------------------------------------------|---------------------------------------------------------------------------------|
| [optparse](https://cran.r-project.org/web/packages/optparse/index.html) | [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html) |
| [ggplot2](https://ggplot2.tidyverse.org/)                               | [lattice](https://cran.r-project.org/web/packages/lattice/index.html)           |
| [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html) | [UpSetR](https://cran.r-project.org/web/packages/UpSetR/README.html)            |
| [scales](https://cran.r-project.org/web/packages/scales/index.html)     | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)       |
| [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html) | [vsn](https://bioconductor.org/packages/release/bioc/html/vsn.html)             |

### Linux utilities

Standard linux tools including `cut`, `awk`, `sort`, `mv`, `touch`, `echo`, `mkdir`, `paste`, `cp`, `ln`, `grep` are also used throughout the pipeline.
