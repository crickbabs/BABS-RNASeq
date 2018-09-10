
## Creating a sample design file

The design file has to be a `csv` (not `tsv`).
The file must contain a column named **sample** which associates a name of your choice to each sample.
The file must also contain the location of the fastq file for each sample.
These columns must be named **file** for a single-end experiment, and **file1** and **file2** for a paired-end experiment.
All the other columns contain the different conditions of the experiment (treatment, time...).
For example, the design file for the single end test data contains:

	sample,file,treatment
	control1,fastq/single_end/sample1.fastq.gz,untreated
	control2,fastq/single_end/sample2.fastq.gz,untreated
	treatment2,fastq/single_end/sample5.fastq.gz,treated

And the one for the paired end test data contains:

	sample,file1,file2,treatment
	control2,fastq/paired_end/sample2_r1.fastq.gz,fastq/paired_end/sample2_r2.fastq.gz,untreated
	control3,fastq/paired_end/sample3_r1.fastq.gz,fastq/paired_end/sample3_r2.fastq.gz,untreated
	treatment3,fastq/paired_end/sample6_r1.fastq.gz,fastq/paired_end/sample6_r2.fastq.gz,treated

You can now [configure the pipeline][url_doc_config].

[url_doc_config]: https://github.com/crickbabs/BABS-RNASeq/blob/master/docs/config.md

