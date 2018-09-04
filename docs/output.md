
## Results and interpretation

The following directories will be created in the current directory after the pipeline has finished.

### `quality_control/`

|Directory|Description|
|---|---|
|`quality_control/fastqc/raw/`|[FastQC][url_fastqc] results **before** adapter trimming|
|`quality_control/fastqc/cutadapt/`|[FastQC][url_fastqc] results **after** adapter trimming|
|`quality_control/picard/`|[Picard][picard] metrics|
|`quality_control/rnaseqc/`|[RNA-SeQC][url_rnaseqc] metrics|
|`quality_control/rnaseqc/`|[RSeQC][url_rseqc] package metrics|
|`quality_control/multiqc_report.html/`|[MultiQC][url_multiqc] report of all the metrics|

### `results/`

|Directory|Description|
|---|---|
|`results/analysis/`|count matrix and PCA results|
|`results/cutadapt/`|trimmed FASTQ files from [cutadapt][url_cutadapt]|
|`results/rsem/`|results of the alignment from [RSEM][url_rsem] and [STAR][url_star]|

[url_picard]: https://broadinstitute.github.io/picard/index.html
[url_cutadapt]: https://cutadapt.readthedocs.io/en/stable
[url_star]: https://github.com/alexdobin/STAR
[url_rsem]: https://github.com/deweylab/RSEM
[url_rseqc]: http://rseqc.sourceforge.net
[url_multiqc]: http://multiqc.info

