
## Install Nextflow

To run Nextflow it needs to be executable on the command-line. On an environment module system such as the one at The Francis Crick Institute this can be achieved by running the following command:

```bash
module load nextflow/0.30.2
```

**You need Nextflow version >= 0.30.2 to run this pipeline.**

See [nextflow.io](https://www.nextflow.io/) for further information on how to install Nextflow.

## Obtaining the pipeline

There are various ways in which you can obtain the pipeline itself, however we recommend that you obtain a local copy of the pipeline by running the following [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) command in a directory where you want to perform the analysis:

```bash
git clone https://github.com/crickbabs/BABS-ATACSeqPE
```

The Nextflow pipeline and associated config and executable files will appear in the `BABS-ATACSeqPE/` directory.

```bash
cd BABS-ATACSeqPE
```

You can now configure the pipeline to run on a Linux system of your choice. See [Pipeline configuration](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/config.md) and [Reference genome](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/genome.md) sections.
