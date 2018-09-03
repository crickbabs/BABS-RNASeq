
# ![BABS-ATACSeqPE](https://raw.githubusercontent.com/crickbabs/BABS-RNASeq/master/docs/images/BABS-RNASeq_logo.png)

## Introduction

A [Nextflow](https://www.nextflow.io/) pipeline for processing RNASeq sequencing data.

The pipeline was written by [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`cutadapt`](http://cutadapt.readthedocs.io/en/stable/installation.html))
3. Alignment and quantification ([`RSEM`](https://deweylab.github.io/RSEM/), [`STAR`](https://github.com/alexdobin/STAR))
4. Various metrics with:
	* ([`picard`][url_picard])
		* Library complexity ([`EstimateLibraryComplexity`][url_picard_complexity])
    * ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
5. Filtering to remove:
    * reads mapping to mitochondrial DNA ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
    * reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that arent marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads containing > 3 mismatches in either read of the pair ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    * reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    * reads that are soft-clipped ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    * reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html))
    * reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html))
    * reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html))
6. Merge alignments at replicate and sample level ([`picard`](https://broadinstitute.github.io/picard/))
    * Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    * Remove duplicate reads ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * Create normalised bigWig files scaled to 1 million mapped read pairs ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`wigToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    * Call broad peaks ([`MACS2`](https://github.com/taoliu/MACS))
    * Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
    * Merge peaks across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
    * Count reads in merged peaks from replicate-level alignments ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
    * Differential binding analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
7. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
8. Collect and present QC at the raw read, alignment and peak-level ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

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

# ![BABS-RNASeqPE directed acyclic graph](https://raw.githubusercontent.com/crickbabs/BABS-ATACSeqPE/master/docs/images/BABS-RNASeqPE_dag.png)

## Credits

The pipeline was written by the [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

The pipeline was developed by [Gavin Kelly](mailto:gavin.kelly@crick.ac.uk), [Harshil Patel](mailto:harshil.patel@crick.ac.uk), [Nourdine Bah](mailto:nourdine.bah@crick.ac.uk) and [Philip East](mailto:philip.east@crick.ac.uk).

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

[url_babs]: https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/
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
[url_deseq2]: https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
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

[dag]: png/dag.png
[pca]: png/single_end_pca.png
[onto]: png/single_end_ontology.png

