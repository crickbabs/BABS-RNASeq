
## Run pipeline

When you have [installed][url_doc_install] the pipeline, [created a design file][url_doc_design] and [configured][url_doc_config] the pipeline you can run it with the command line.

### With a parameter file

If you want to run the workflow with a parameter file the command is:

	$ nextflow run main.nf -params-file examples/params/single_end.yml

, or:

	$ nextflow run main.nf -params-file examples/params/paired_end.yml

### With command line arguments

The native Nextflow's command line arguments are preceded by `-`, for example `-params-file`.
However, the `nextflow` command accepts any parameter if it is preceded by `--`.
You can also specify parameters at the command-line to override those provided in the config files. A list of these can be obtained with:

	$ nextflow run main.nf --help

For example, the `genome_version` parameters can be passed to `nextflow` like this:

	$ nextflow run main.nf --genome_version GRCh38

Thus, in this way, all the parameters of the pipeline can be passed via a shell script:

	$ nextflow run main.nf --genome_version GRCh38 --genome_release 86 --strandedness none --design examples/designs/single_end.csv --conda /camp/stp/babs/working/software/anaconda/envs/rnaseq_pipeline

, or:

	$ nextflow run main.nf --genome_version GRCh38 --genome_release 86 --strandedness none --design examples/designs/paired_end.csv --conda /camp/stp/babs/working/software/anaconda/envs/rnaseq_pipeline

[url_doc_install]: https://github.com/crickbabs/BABS-RNASeq/blob/master/docs/install.md
[url_doc_design]: https://github.com/crickbabs/BABS-RNASeq/blob/master/docs/design.md
[url_doc_config]: https://github.com/crickbabs/BABS-RNASeq/blob/master/docs/config.md

