#!/usr/bin/env nextflow

// set dsl version 
nextflow.enable.dsl=2

// include modules for workflow
include { 
	apa_annotation; 
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

	// annotation 
	apa_annotation(benchmark_rscript, gtf_file, apa_output, barcode_file, sample_name, apa_method, params.core_num)
	
}

