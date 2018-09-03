
## Creating a sample design file

You will need to create a design file with information about the samples in your experiment before running the pipeline. It has to be a comma-separated file with 5 columns, and a header row as shown in the example below:

```bash
sample,replicate,run,fastq_1,fastq_2
control,1,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
control,2,1,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
treatment,1,1,AEG588A3_S5_L003_R1_001.fastq.gz,AEG588A3_S5_L003_R2_001.fastq.gz
treatment,2,1,AEG588A4_S3_L002_R1_001.fastq.gz,AEG588A4_S3_L002_R2_001.fastq.gz
treatment,2,2,AEG588A5_S4_L002_R1_001.fastq.gz,AEG588A5_S4_L002_R2_001.fastq.gz
```

| Column      | Description                                                                                                                               |
|-------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| `sample`    | Group identifier for sample.                                                                                                              |
| `replicate` | Integer representing replicate number.                                                                                                    |
| `run`       | Integer representing the number of times the same library has been sequenced. This will be used later for merging at the replicate-level. |
| `fastq_1`   | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz".                                             |
| `fastq_2`   | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz".                                             |
