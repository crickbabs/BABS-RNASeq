#!/usr/bin/env nextflow

/*
 * gavin.kelly@crick.ac.uk
 * harshil.patel@crick.ac
 * nourdine.bah@crick.ac.uk
 * philip.east@crick.ac.uk
 */

import java.nio.file.Paths

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            FUNCTIONS                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

/* ~~~~~~~~~~~~~~~~~~~ */
def absolute_path(path) {

	def f = new File(path)

	if ( ! f.isAbsolute() ) {
		return Paths.get(workflow.projectDir.toString(), path).toString()
	} else {
		return path
	}
}

/* ~~~~~~~~~~~~~~~~~~~~~ */
def check_file_path(path) {

	// a java file object
	def f = new File(path)
	def abs_path = ""
	
	// get absolute path
	if ( ! f.isAbsolute() ) {
		abs_path = Paths.get(workflow.launchDir.toString(), path).toString()
	} else {
		abs_path = path
	}

	// is the file exists ?
	def file = new File(abs_path)
	assert file.exists() : "Error " + file + " does not exist."

	return file.getAbsolutePath()
}

/* ~~~~~~~~~~~~~~~~~~~~~ */
def is_single_end(design) {

	// get header
	def header = ""
	new File(design).withReader{ header = it.readLine() }

	// check columns names
	def columns = header.split(",")
	def file_col = columns.findIndexOf{ it == "file" } >= 0 ? true : false 
	def file1_col = columns.findIndexOf{ it == "file1" } >= 0 ? true : false 
	def file2_col = columns.findIndexOf{ it == "file2" } >= 0 ? true : false

	// determine single or paired end
	if ( file1_col & file2_col & !file_col ) {
		return false
	} else if ( file_col & !file1_col & !file2_col ) {
		return true
	} else {
		println "Error the design file is not well formatted."
		System.exit(1)
	}
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
def get_rough_read_length(read_length) {
	
	// three babs star indices have been built with these read length parameters
	def starIndexReadLengths = [50, 75, 100]
	
	// take the index with the closest read length to the experiment's
	def diffs = []
	starIndexReadLengths.each() { length ->
		diff = (length - read_length.toInteger()).abs()
		diffs.add(diff)
	}
	def index = diffs.findIndexValues() { i -> i == diffs.min() }[0]
	def rough_read_length = starIndexReadLengths[index.toInteger()]

	return rough_read_length
}

/* ~~~~~~~~~~~~~~~~ */
def create_directory(path) {
	def dir = new File(path)
	if ( ! dir.exists() ) {
		dir.mkdirs()
	}
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          INPUT CHANNELS                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// parameters
design_filepath = check_file_path(params.design)
single_end = is_single_end(design_filepath)

if (single_end) {

	samples = Channel
					.fromPath(design_filepath)
					.splitCsv(header: true)
					.map{
						row ->
						[ name: row.sample, file: check_file_path(row.file) ] }
	
	Channel
		.fromPath(design_filepath)
		.splitCsv(header: true)
		.map{ row -> file(row.file) }
		.into{ fastq_files; fastq_files_fastqc; fastq_files_metadata }

} else {

	samples = Channel
					.fromPath(design_filepath)
					.splitCsv(header: true)
					.map{
						row ->
						[
							name: row.sample,
							file1: check_file_path(row.file1),
							file2: check_file_path(row.file2)
						]
					}

	Channel
		.fromPath(design_filepath)
		.splitCsv(header: true)
		.map{ row -> [row.file1, row.file2] }
		.flatten()
		.map{ file(it) }
		.into{ fastq_files; fastq_files_fastqc; fastq_files_metadata }
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                           READ LENGTH                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process read_length {

	tag { fastq }

	input:
		file(fastq) from fastq_files

	output:
		stdout into lengths

	script:
		template "read_length.sh"
}

lengths
	.map{ it.toInteger() }
	.into{ lengths_sum; lengths_count }

lengths_sum
	.sum()
	.set{ lengths_sum }

lengths_count
	.count()
	.set{ lengths_count }

lengths_sum
	.concat( lengths_count )
	.collect()
	.map{ it[0] / it[1] }
	.into{ read_length_ch; read_length_report }

read_length =
	read_length_ch
		.collect()
		.get()
		.getAt(0)

rough_read_length = get_rough_read_length(read_length)



///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         OTHER PARAMETERS                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// genome
version = params.genome_version.toString()
release = params.genome_release.toString()
genome = params.genomes[version][release]
fasta = genome.fasta

// annotation
gtf = genome.gtf
bed = genome.bed
refflat = genome.refflat
rrna_list = genome.rrna_list
rrna_interval_list = genome.rrna_interval_list
rnaseqc_gtf = genome.rnaseqc_gtf
index = genome.rsem.rsem_star[ rough_read_length + "bp" ].index

// directories
r_script_dir = absolute_path("scripts/r")
qc_dir = params.qc_dir
results_dir = params.results_dir

// conf
multiqc_conf = absolute_path( "multiqc/conf.yml" )

// anaconda 
conda_env = params.conda


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            WORKFLOW                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

////////////////
process fastqc {

	tag { name }

	input:
		file(fastq) from fastq_files_fastqc

	output:
		file("*fastqc*") into fastqc
	
	script:
	
		name = fastq.toString().replaceFirst(".fastq.gz", "")

		template "fastqc.sh"
}


//////////////////
process md5_fastq {

	input:
		file(fastq) from fastq_files_metadata

	output:
		file("md5.txt") into md5_raw

	"""
	md5sum $fastq > md5.txt
	"""
}

//////////////////
process cutadapt {

	tag { name }

	input:
		val sample from samples
	
	output:
		set val(name), file('*.log') into cutadapt_log_ungrouped
		set val(name), file('*.cutadapt.fastq.gz') \
			into trimmed_ungrouped_fastq 

	script:

		name = sample["name"]
		
		if (single_end) {

			fastq = sample["file"]
			output = name + ".cutadapt.fastq.gz"
			logfile = name + ".log"

			template "cutadapt/single_end.sh"

		} else {

			fastq1 = sample["file1"]
			fastq2 = sample["file2"]
			output1 = name + "_1.cutadapt.fastq.gz"
			output2 = name + "_2.cutadapt.fastq.gz"
			logfile = name + ".log"

			template "cutadapt/paired_end.sh"
		}
}



grouped_trimmed = Channel.create()
grouped_trimmed = trimmed_ungrouped_fastq.groupTuple()
grouped_log = Channel.create()
grouped_log = cutadapt_log_ungrouped.groupTuple()

process unify_lanes {
	input:
		set val(name), file(group_of_fastq) from grouped_trimmed
	output:
		set val(name), file("${name}.fastq.gz") into cutadapt_fastq, cutadapt_rsem

		
	"""
	cat ${group_of_fastq} > ${name}.fastq.gz
	"""
}

process unify_cutlog {
	input:
		set val(name), file(group_of_log) from grouped_log
	output:
		set val(name), file("${name}.log") into cutadapt_log
		
	"""
	cat ${group_of_log} > ${name}.log
	"""
}



/////////////////////////
process fastqc_cutadapt {

	tag { name }

	input:
		set val(name), file(fastq) from cutadapt_fastqc

	output:
		file "*fastqc*" into fastqc_cutadapt
	
	script:
	
		name = fastq.toString().replaceFirst(".fastq.gz", "")

		template "fastqc.sh"
}

//////////////
process rsem {

	tag { name }

	input:
		set val(name), file(fastq) from cutadapt_rsem

	output:
		set val(name), file("*.transcript.bam") into rsem_transcript
		set val(name), file("*.stat") into rsem_stat
		set val(name), file("*.STAR.genome.bam") into rsem_genome, rsem_genome_metadata
		set val(name), file("*.results") \
			into \
				rsem_results,
				rsem_results_analysis,
				rsem_results_metadata

	script:
		strandedness = params.strandedness
		cpus = task.cpus
		if (single_end) {

			template "rsem/star/single_end.sh"

		} else {

			template "rsem/star/paired_end.sh"
		}
}

process rsem_md5 {

	tag { name }

	input:
		set val(name), file(results) from rsem_results_metadata

	output:
		file("md5.txt") into md5_processed

	"""
	md5sum $results > md5.txt
	"""
}

////////////////////
process sort_index {

	tag { name }

	input:
		set val(name), file(bam) from rsem_genome

	output:
		set val(name), file("*.bam"), file("*.bai") \
			into \
				rsem_genome_indexed,
				rsem_genome_indexed_multiqc

	script:

		filename = name + ".sorted.bam"
		cpus = task.cpus

		template "samtools/sort_index.sh"
}

////////////////////
process insert_size {

	tag { name }

	input:
		set val(name), file(bam) from rsem_genome_metadata

	output:
		file("insert_size.txt") into insert_size

	"""
	samtools stats $bam  | grep "insert size " | sed --expression "s/SN/$name/" | awk -v FS='\t' -v OFS='\t' 'NR % 2 == 1 { o=\$0 ; next } { print o , \$3 }' > insert_size.txt
	"""
}

////////////////////
process metadata {

 	input:
	file 'raw.txt' from md5_raw.collectFile(name: 'raw.txt')
 	file 'proc.txt' from md5_processed.collectFile(name: 'proc.txt')
	file 'insert.txt' from insert_size.collectFile(name: 'insert.txt')

	output: file 'metadata.txt'

	"""
	echo "genome: " > metadata.txt
	echo "  version: " $version >> metadata.txt
	echo "  release: " $release >> metadata.txt
	echo "  genome: " $genome >> metadata.txt
        echo "  fasta: " $fasta >> metadata.txt
	echo "annotation: " >> metadata.txt
	echo "  gtf: " $gtf >> metadata.txt
	echo "  bed: " $bed >> metadata.txt
	echo "  refflat: " $refflat >> metadata.txt
	echo "  rrna_list: " $rrna_list >> metadata.txt
	echo "  rrna_interval_list: " $rrna_interval_list >> metadata.txt
	echo "  rnaseqc_gtf: " $rnaseqc_gtf >> metadata.txt
	echo "  index: " $index  >> metadata.txt
	awk 'BEGIN{print "raw_md5: "}{print "  - name: "\$2"\\n    md5: "\$1}' raw.txt >> metadata.txt
	awk 'BEGIN{print "processed_md5: "}{print "  - name: "\$2"\\n    md5: "\$1}' proc.txt >> metadata.txt
	awk -F '\t' 'BEGIN{print "insert: "}{print "  - name: "\$1"\\n    average: "\$3"\\n    sd: "\$4}' insert.txt >> metadata.txt
 	"""

 }



///////////////
process group {

	tag { name }

	input:
		set val(name), file(bam), file(bai) from rsem_genome_indexed

	output:
		set val(name), file("*.bam") into group

	script:

		tmp_dirname = "tmp"
		filename = name + ".rg.bam"

		template "picard/group.sh"

}

////////////////////
process duplicates {

	tag { name }

	input:
		set val(name), file(bam) from group

	output:
		set val(name), file("*.marked_duplicates") into duplicates_stats
		set val(name), file("*.bam") \
			into \
				duplicates_index,
				duplicates_complexity,
				duplicates_multimetrics,
				duplicates_rnaseqmetrics,
				duplicates_infer_experiment,
				duplicates_junction_annotation,
				duplicates_junction_saturation,
				duplicates_mismatch_profile,
				duplicates_read_distribution

	script:

		tmp_dirname = "tmp"
		filename = name + ".dupmarked.bam"
		metrics_filename = name + ".marked_duplicates"

		template "picard/duplicates.sh"
}

///////////////
process index {

	// custom label
	tag { name }

	input:
		set val(name), file(bam) from duplicates_index

	output:
		set val(name), file(bam), file("*.bai") \
			into \
				duplicates_bai,
				duplicates_bai_rnaseqc,
				duplicates_bai_multiqc

	script:

		template "samtools/index.sh"
}

////////////////////
process complexity {

	tag { name }

	input:
		set val(name), file(bam) from duplicates_complexity

	output:
		set val(name), file("*.complexity") into complexity

	script:

		tmp_dirname = "tmp"
		metrics_filename = name + ".complexity"

		template "picard/complexity.sh"
}

///////////////////////
process rnaseqmetrics {

	tag { name }

	input:
		set val(name), file(bam) from duplicates_rnaseqmetrics

	output:
		set val(name), file("*.rnaseqmetrics") into rnaseqmetrics

	script:

		tmp_dirname = "tmp"
		metrics_filename = name + ".rnaseqmetrics"

		template "picard/rnaseqmetrics.sh"
}

//////////////////////
process multimetrics {

	tag { name }

	input:
		set val(name), file(bam) from duplicates_multimetrics

	output:
		set val(name), file("*.pdf") into multimetrics_pdf
		set val(name), file("*_metrics") into multimetrics_metrics

	script:

		tmp_dirname = "tmp"
		metrics_filename = name + ".multimetrics"

		template "picard/multimetrics.sh"
}

//////////////////////////
process infer_experiment {

	tag { name }

	input:
		set val(name), file(bam) from duplicates_infer_experiment
	
	output:
		set val(name), file("*.infer_experiment*") into infer_experiment
	
	script:

		metrics_filename = name + ".infer_experiment"

		template "rseqc/infer_experiment.sh"
}

/////////////////////////////
process junction_annotation {

	tag { name }

	input:
		set val(name), file(bam) from duplicates_junction_annotation

	output:
		set val(name), file("*.junction_annotation*") into junction_annotation
	
	script:

		metrics_filename = name + ".junction_annotation"

		template "rseqc/junction_annotation.sh"
}

/////////////////////////////
process junction_saturation {

	tag { name }

	input:
		set val(name), file(bam) from duplicates_junction_saturation

	output:
		set val(name), file("*.junction_saturation*") into junction_saturation

	script:

		metrics_filename = name + ".junction_saturation"

		template "rseqc/junction_saturation.sh"
}

//////////////////////////
process mismatch_profile {

	tag { name }

	input:
		set val(name), file(bam) from duplicates_mismatch_profile

	output:
		set val(name), file("*.mismatch_profile*") into mismatch_profile

	script:

		metrics_filename = name + ".mismatch_profile"

		template "rseqc/mismatch_profile.sh"
}




///////////////////////////
process read_distribution {

	tag { name }

	input:
		set val(name), file(bam) from duplicates_read_distribution

	output:
		set val(name), file("*.read_distribution*") into read_distribution
	
	script:

		metrics_filename = name + ".read_distribution"

		template "rseqc/read_distribution.sh"
}

/////////////////////////////////////
process transcript_integrity_number {

	tag { name }

	input:
		set val(name), file(bam), file(bai) from duplicates_bai

	output:
		set val(name), file("*.{xls,txt}") into transcript_integrity_number

	script:

		template "rseqc/transcript_integrity_number.sh"
}

/////////////////
process rnaseqc {

	tag { name }

	input:
		set val(name), file(bam), file(bai) from duplicates_bai_rnaseqc

	output:
		set val(name), file("*rnaseqc*") into rnaseqc

	script:

		if (single_end) {

			metrics_filename = name + ".rnaseqc"

			template "rnaseqc/single_end.sh"

		} else {

			metrics_filename = name + ".rnaseqc"

			template "rnaseqc/paired_end.sh"
		}
}

////////////////////
process get_counts {

	tag { name }

	input:
		set val(name), file(results) from rsem_results_analysis

	output:
		set val(name), file("*.tsv") into counts
	
	script:

		filename = name + ".genes.results"
		counts = name + ".tsv"

		template "get_counts.sh"
		
}

//////////////////////
process build_matrix {

	input:
		file("*") from counts.map{x->x[1]}.collect()
	
	output:
		file("*.rda") into rda_matrix
		file("*.csv") into matrix
	
	shell:

		output_name = "matrix"

		"""
		Rscript ${r_script_dir}/matrix.r \
			--results-dir . \
			--design-file ${design_filepath} \
			--gtf-file ${gtf} \
			--output-name ${output_name}
		"""
}

////////////////////
process create_dds {

	input:
		file(rda) from rda_matrix
	
	output:
		file("*.rda") into (dds_vst, dds_rlog)
	
	script:
		
		rda_in = rda
		rda_out = "dds.rda"

		"""
		Rscript ${r_script_dir}/dds.r \
			--rda-se ${rda_in} \
			--formula "~ 1" \
			--ncores 32 \
			--rda-out ${rda_out}
		"""
}

////////////////////
process create_vst {

	input:
		file(rda) from dds_vst
	
	output:
		file("*.rda") into vst
	
	script:
		
		rda_in = rda
		rda_out = "vst.rda"

		"""
		Rscript ${r_script_dir}/vst.r ${rda_in} ${rda_out}
		"""
}

/////////////////////
process create_rlog {

	input:
		file(rda) from dds_rlog
	
	output:
		file("*.rda") into rlog
	
	script:
		
		rda_in = rda
		rda_out = "rlog.rda"

		"""
		Rscript ${r_script_dir}/rlog.r ${rda_in} ${rda_out}
		"""
}

////////////////
process do_pca {

	input:
		file(rda) from rlog

	output:
		file("*.rda") into pca
	
	script:

		rda_in = rda
		rda_out = "pca.rda"

		"""
		Rscript ${r_script_dir}/pca.r ${rda_in} ${rda_out}
		"""
}

/////////////////////
process pca_to_json {

	input:
		file(rda) from pca

	output:
		file("*.json") into pca_json
	
	script:

		json = "rnaseq_pca.json"

		"""
		Rscript ${r_script_dir}/pca_to_json.r \
			--rda ${rda} \
			--design ${design_filepath} \
			> ${json}
		"""
}

/////////////////
process multiqc {

	conda conda_env

	input:

		// cutadapt
		file("*") from cutadapt_log.map{x->x[1]}.collect()

		// fastqc
		file("*") from fastqc.collect()
		file("*") from fastqc_cutadapt.collect()

		// rsem
		file("*") from rsem_genome_indexed_multiqc.map{x->x[1]}.collect()
		file("*") from rsem_transcript.map{x->x[1]}.collect()
		file("*") from rsem_results.map{x->x[1]}.collect()
		file("*") from rsem_stat.map{x->x[1]}.collect()

		// picard
		file("*") from duplicates_stats.map{x->x[1]}.collect()
		file("*") from duplicates_bai_multiqc.map{x->x[1]}.collect()
		file("*") from complexity.map{x->x[1]}.collect()
		file("*") from rnaseqmetrics.map{x->x[1]}.collect()
		file("*") from multimetrics_metrics.map{x->x[1]}.collect()

		// rseqc
		file("*") from infer_experiment.map{x->x[1]}.collect()
		file("*") from junction_annotation.map{x->x[1]}.collect()
		file("*") from junction_saturation.map{x->x[1]}.collect()
		file("*") from mismatch_profile.map{x->x[1]}.collect()
		file("*") from read_distribution.map{x->x[1]}.collect()
		file("*") from transcript_integrity_number.map{x->x[1]}.collect()

		// rnaseqc
		file("*") from rnaseqc.map{x->x[1]}.collect()

		// pca
		file("*") from pca_json.collect()

	output:
		file "multiqc_data" into multiqc_data
		file "multiqc_report.html" into multiqc_report

	script:

		multiqc_template = "custom"

		template "multiqc.sh"
}

workflow.onComplete {
}

