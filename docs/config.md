
## Specifying Analysis Parameters

All the parameters except `modules` are mandatory.

|Parameter|Description|Required by|
|---|---|---|
|`genome_version`|the genome version of the organism in which the experiment has been done (`36` or `38` most of time). This parameter is required to build the path of the STAR genome index|`STAR`|
|`genome_release`|the genome release of the organism in which the experiment has been done (`86` or `89` most of time). This parameter is required to build the path of the STAR genome index|`STAR`|
|`design`|the path of the `csv` [design file][url_doc_design]|every processes because the sample column values will be used to name the output files, the experimental condition columns values will be used by `DESeq2`|
|`strandedness`|the strandedness of the reads (`none`, `forward`, `reverse`), but the pipeline determines this by itself anyway|`STAR`|
|`conda`|the absolute path to the conda environment|`STAR`|
|`modules`|the modules you want to use if you want a different version for the default||

There are two ways to specify the analysis parameters : with a parameter `yml` file or via command line arguments (see [Running the pipeline][url_doc_usage]).

[url_doc_usage]: https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/usage.md
[url_doc_design]: https://github.com/crickbabs/BABS-ATACSeqPE/blob/master/docs/design.md

