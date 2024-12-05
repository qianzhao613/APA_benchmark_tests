#!/usr/bin/env nextflow

// set dsl version 
nextflow.enable.dsl=2

// include modules for workflow
include { 
	apa_annotation; 
	overlap_genes;
	nanopore_filter;
	identification; 
	nanopore_normalization; 
	apa_normalization;
	quantification; 
	differential_expression 
	} from "$launchDir/modules/apa_benchmark_process"

// define benchmark workflow
workflow {
	// Define inputs
	def sample_name = params.sample_name
	def apa_method = params.apa_method
	def gtf_file = file(params.gtf_file)
	def benchmark_rscript = file(params.benchmark_rscript)
	apa_output = channel.fromPath("${launchDir}/data/${sample_name}_illumina/*_${apa_method}_*", type: "dir", checkIfExists: true)
	barcode_file = channel.fromPath("${launchDir}/data/${sample_name}_illumina/${sample_name}_barcodes.tsv*", type: "file", checkIfExists: true)

	// Emit win_sizes to channel
	if (params.win_size instanceof String && params.win_size.startsWith('[') && params.win_size.endsWith(']')) { //
		def parsed_win_size = params.win_size.replaceAll(/\[|\]/, '').split(',').collect { it.trim().toInteger() }
		def (start, end, step) = parsed_win_size
		win_sizes = (start..end).step(step.toInteger()).toList()
	} else {
	// Treat win_size as a single number
		win_sizes = [params.win_size]
	}
	win_size_channel = channel.from(win_sizes)

	// annotation 
	//apa_annotation(benchmark_rscript, gtf_file, apa_output, barcode_file, sample_name, apa_method, params.core_num, params.pipe_folder)
	apa_annotation = channel.fromPath("${params.annotation_folder}/${sample_name}/${apa_method}_peak_annotations_by_SCAPE.qs", type: "file", checkIfExists: true)
	method_list = channel.of("Sierra", "scAPA", "polyApipe", "scAPAtrap", "SCAPTURE", "MAAPER","SCAPE")

	def overlap_genes_file = file("${params.annotation_folder}/${sample_name}/${sample_name}_overlap_gene_list.qs")

	if (overlap_genes_file.exists()) {
		// If file exists, create a channel from the file
		overlap_genes_channel = channel.fromPath(overlap_genes_file)
		overlap_genes_channel.view { "File exists, importing: $it" }
	} else {
		// If file does not exist, run the process
		overlap_genes_channel = overlap_genes(benchmark_rscript, method_list, sample_name, params.core_num, params.pipe_folder)
		overlap_genes_channel.view { "File does not exist, running process to generate: $it" }
	}

	// normalization
	// Check if the file exists
	def nanopore_normalization_file = file("${params.nanopore_normalization_folder}/${sample_name}_nanopore_normalized_isomatrix.qs")
	def apa_count_file = file("${params.apa_normalization_folder}/${sample_name}_${apa_method}_count_matrix.qs")
	def apa_normalization_file = file("${params.apa_normalization_folder}/${sample_name}_${apa_method}_normalized_by_ranger_matrix.qs")


	if (nanopore_normalization_file.exists()) {
		// If file exists, create a channel from the file
		nanopore_normalization_channel = Channel.fromPath(nanopore_normalization_file)
		nanopore_normalization_channel.view { "File exists, importing: $it" }
	} else {
		// If file does not exist, run the process
		nanopore_normalization_channel = nanopore_normalization(benchmark_rscript, sample_name, params.core_num, params.pipe_folder)
		nanopore_normalization_channel.view { "File does not exist, running process to generate: $it" }
	}
	
	if (apa_count_file.exists() && apa_normalization_file.exists()) {
		// If file exists, create a channel from the file
		apa_count_channel = channel.fromPath(apa_count_file)
		apa_normalization_channel = channel.fromPath(apa_normalization_file)

		apa_count_channel.view { "File1: $it" }
		apa_normalization_channel.view { "File2: $it" }
	} else {
		// If file does not exist, run the process
		apa_matrix_channel = apa_normalization(benchmark_rscript, apa_output, apa_annotation, sample_name, apa_method, params.core_num, params.pipe_folder) 

		apa_count_channel = apa_matrix_channel.map { it[0] }
		apa_normalization_channel = apa_matrix_channel.map { it[1] }

		apa_matrix_channel.view { result -> 
    "File does not exist, running process to generate: File1: ${result[0]}, File2: ${result[1]}" }
	}
	
	// nanopore filter
	def nanopore_filter_file = file("${params.annotation_folder}/${sample_name}/${sample_name}_${params.cell_type}_normalized_nanopore_isoform_table.qs")

	if (nanopore_filter_file.exists()) {
		// If file exists, create a channel from the file
		nanopore_filter_channel = Channel.fromPath(nanopore_filter_file)
		nanopore_filter_channel.view { "File exists, importing: $it" }
	} else {
		// If file does not exist, run the process
		nanopore_filter_channel = nanopore_filter(benchmark_rscript, nanopore_normalization_channel, sample_name, params.core_num, params.cell_type, params.pipe_folder)
		nanopore_filter_channel.view { "File does not exist, running process to generate: $it" }
	}
	
	if (params.cv_cutoff instanceof String && params.cv_cutoff.startsWith('[') && params.cv_cutoff.endsWith(']')) { //
		def parsed_cv_cutoff = params.cv_cutoff.replaceAll(/\[|\]/, '').split(',').collect { it.trim().toFloat() }
		def (start, end, step) = parsed_cv_cutoff
		cv_cutoffs = (start..end).step(step.toInteger()).toList()
	} else {
	// Treat win_size as a single number
		cv_cutoffs = [params.cv_cutoff]
	}

	// Emit win_sizes to channel
	cv_cutoff_channel = channel.from(cv_cutoffs)
	win_cv = win_size_channel.flatten().combine(cv_cutoff_channel.flatten())

	// Benchmark 1: Identification
	identification(benchmark_rscript, apa_count_channel.first(), sample_name, apa_method, win_cv, overlap_genes_channel.first(), params.core_num, params.pipe_folder)

	// Benchmark 2: Quantification

	quantification(benchmark_rscript, nanopore_normalization_channel.first(), apa_normalization_channel.first(), sample_name, apa_method, win_cv, overlap_genes_channel.first(), params.cell_type, params.core_num, params.pipe_folder)

	// Benchmark 3: Differential expression
	differential_expression(benchmark_rscript, apa_count_channel.first(), sample_name, apa_method, win_size_channel, overlap_genes_channel.first(), params.core_num, params.pipe_folder)
}


workflow.onComplete {
    println "Pipeline successfully completed at: $workflow.complete"
}

workflow.onError {
    println "Pipeline stopped with following error: $workflow.errorMessage"
}

