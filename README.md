
# ![BABS-ATACSeqPE][logo]

## Introduction

A [Nextflow][nextflow] pipeline for processing RNASeq sequencing data.

The pipeline was written by [The Bioinformatics & Biostatistics Group][url_babs] at [The Francis Crick Institute][url_crick], London.

## Pipeline summary

1. Raw read QC ([`FastQC`]([url_fastqc]))
2. Adapter trimming ([`cutadapt`][url_cutadapt])
3. Alignment and quantification ([`RSEM`][url_rsem], [`STAR`][url_star])
4. Sorting and indexing ([`SAMtools`][url_samtools])
5. Quality control metrics:
	* [`picard`][url_picard]:
		* Groups ([`AddOrReplaceReadGroups`][url_picard_group])
		* Duplicates ([`MarkDuplicates`][url_picard_duplicate])
		* Library complexity ([`EstimateLibraryComplexity`][url_picard_complexity])
		* Various metrics ([`CollectRnaSeqMetrics`][url_picard_rnaseqmetrics], [`CollectMultipleMetrics`][url_picard_multimetrics])
	* [`RSeQC`][url_rseqc]:
		* Samples quality ([`infer_experiment.py`][url_rseqc_infer_experiment], [`read_distribution.py`][url_rseqc_read_distribution], [`tin.py`][url_rseqc_tin])
		* Alternative splicing ([`junction_annotation.py`][url_rseqc_junction_annotation], [`junction_saturation.py`][url_rseqc_junction_saturation])
		* Mismatch ([`mismatch_profile.py`][url_rseqc_mismatch_profile])
	* [`RNA-SeQC`][url_rnaseqc]
6. Preparation for statistical analysis:
	* Create a count matrix ([`R`][url_r], [`SummarizedExperiment`][url_summarized_experiment])
	* Perform a principal component analysis ([`R`][url_r], [`DESeq2`][url_deseq2])
8. Collect and present a report ([`MultiQC`][url_multiqc])

## Documentation

The documentation for the pipeline can be found in the `docs/` directory:

1. [Installation](docs/install.md)
2. [Pipeline configuration](docs/config.md)
3. [Reference genome](docs/genome.md)
4. [Design file](docs/design.md)
5. [Running the pipeline](docs/usage.md)
6. [Output and interpretation of results](docs/output.md)
7. [Troubleshooting](docs/troubleshooting.md)

## Pipeline DAG

# ![BABS-RNASeqPE directed acyclic graph][dag]

## Credits

The pipeline was written by the [The Bioinformatics & Biostatistics Group][url_babs] at [The Francis Crick Institute][url_crick], London.

The pipeline was developed by [Gavin Kelly](mailto:gavin.kelly@crick.ac.uk), [Harshil Patel](mailto:harshil.patel@crick.ac.uk), [Nourdine Bah](mailto:nourdine.bah@crick.ac.uk) and [Philip East](mailto:philip.east@crick.ac.uk).

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

[url_babs]: https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/
[url_fastqc]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
[url_crick]: https://www.crick.ac.uk/
[url_nextflow]: http://www.nextflow.io
[url_nextflow_tuto]: http://www.nextflow.io/docs/latest/getstarted.html#get-started
[url_picard]: https://broadinstitute.github.io/picard/index.html
[url_picard_complexity]: https://broadinstitute.github.io/picard/command-line-overview.html#EstimateLibraryComplexity
[url_picard_duplicate]: https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
[url_picard_group]: https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups
[url_picard_multimetrics]: https://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics
[url_picard_rnaseqmetrics]: https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics
[url_cutadapt]: https://cutadapt.readthedocs.io/en/stable/
[url_star]: https://github.com/alexdobin/STAR
[url_rsem]: https://github.com/deweylab/RSEM
[url_rsem_calculate_expression]: http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html
[url_r]: https://www.r-project.org/
[url_deseq2]: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[url_summarized_experiment]: https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html
[url_samtools]: http://www.htslib.org/doc/samtools.html
[url_fastq_screen]: https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
[url_rseqc]: http://rseqc.sourceforge.net/
[url_rseqc_infer_experiment]: http://rseqc.sourceforge.net/#infer-experiment-py
[url_rseqc_junction_annotation]: http://rseqc.sourceforge.net/#junction-annotation-py
[url_rseqc_junstion_saturation]: http://rseqc.sourceforge.net/#junction-saturation-py
[url_rseqc_mismatch_profile]: http://rseqc.sourceforge.net/#mismatch-profile-py
[url_rseqc_read_distribution]: http://rseqc.sourceforge.net/#read-distribution-py
[url_rseqc_tin]: http://rseqc.sourceforge.net/#tin-py
[url_rnaseqc]: http://archive.broadinstitute.org/cancer/cga/rna-seqc
[url_multiqc]: http://multiqc.info/

[logo]: https://raw.githubusercontent.com/crickbabs/BABS-RNASeq/master/docs/images/BABS-RNASeq_logo.png
[dag]: https://raw.githubusercontent.com/crickbabs/BABS-RNASeq/master/docs/images/dag/dag.png

