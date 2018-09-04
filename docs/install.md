
## Install Nextflow

To run Nextflow it needs to be executable on the command-line. On an environment module system such as the one at The Francis Crick Institute this can be achieved by running the following command:

```bash
module load nextflow/0.30.0
```

**You need Nextflow version >= 0.30.0 to run this pipeline.**

See [nextflow.io][url_nextflow] for further information on how to install Nextflow.

## Obtaining the pipeline

There are various ways in which you can obtain the pipeline itself, however we recommend that you obtain a local copy of the pipeline by running the following [git][url_git] command in a directory where you want to perform the analysis:

```bash
git clone https://github.com/crickbabs/BABS-ATACSeqPE
```

The Nextflow pipeline and associated config and executable files will appear in the `BABS-ATACSeqPE/` directory.

```bash
cd BABS-ATACSeqPE
```
## Install the custom MultiQC module

With a `environment.yml` file, Nextflow is able to create an entire conda environment and execute all the processes inside it. This is probably the most portable way to combine Nextflow and Conda. However, this solution will not work if you want to use a modified version of a software, or just a software not available in Conda's repositories.

Here we want to use custom MultiQC modules, so we will use a conda environment in which those multiqc modules have been installed. This environment is already available for BABS members; the other users will have to create it.

### If you are from BABS

You don't need to do anything. The pipeline is using a `rnaseq_pipeline` conda  environment available in our shared space : `/camp/stp/babs/working/software/anaconda/envs`. Just be sure that anaconda is properly configured.

	$ readlink $HOME/.conda
	/camp/stp/babs/working/$USER/.conda

	$ readlink $HOME/.condarc
	/camp/stp/babs/working/$USER/.condarc

	$ cat /camp/stp/babs/working/$USER/.condarc
	envs_dirs:
	 - /camp/stp/babs/working/software/anaconda/envs
	pkgs_dirs:
	 - /camp/stp/babs/working/software/anaconda/pkgs

	$ cat /camp/stp/babs/working/$USER/.conda/environments.txt
	...
	/camp/stp/babs/working/software/anaconda/envs/rnaseq_pipeline
	...

### If you are not from BABS

You need to create an anaconda environment and install multiqc, other packages, and the multiqc plugins inside this environment. Here is the procedure.

First, you need to load Anaconda:

	$ module load Anaconda2/5.1.0

Then, you create a new environment named `rnaseq_pipeline`:

	$ conda create --yes --name rnaseq_pipeline python=3.6

After that you can install somes packages that the pipeline requires:

	$ conda install --yes --channel bioconda --name rnaseq_pipeline multiqc=1.5
	$ conda install --yes --channel anaconda --name rnaseq_pipeline openblas=0.2.20
	$ conda install --yes --channel bioconda --name rnaseq_pipeline htseq=0.9.1
	$ conda install --yes --channel anaconda --name rnaseq_pipeline pandas=0.23.3

Install the multiqc plugins:

	$ conda install --yes --name rnaseq_pipeline multiqc/multiqc_plugins-1.0-py36_2.tar.bz2

Finally the variable `ANACONDA_ENV` needs to be changed in the `main.nf` file. This variable has to be set to the path of the anaconda environment you just created.


You can now configure the pipeline to run on a Linux system of your choice. See [Pipeline configuration][url_doc_config] and [Reference genome][url_doc_genom] sections.

[url_nextflow]: http://www.nextflow.io
[url_git]: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
[url_doc_config]: https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/config.md)
[url_doc_genome]: https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/genome.md

