
## Run pipeline

When you have [installed](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/install.md) and [configured](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/config.md) the pipeline, prepared the [reference genome](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/genome.md), and [created a design file](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/design.md) you can run it with the command below in the top-level directory of the pipeline:

```bash
nextflow run main.nf --design <DESIGN_FILE_PATH> --genome <GENOME_NAME> -profile <PROFILE_NAME>
```

You can also specify parameters at the command-line to override those provided in the config files. A list of these can be obtained with:

```bash
nextflow run main.nf --help
```

A template script ([`run_pipeline.sh`](https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/run_pipeline.sh)) outlining the steps and commands required to run the pipeline is also available to customise.
