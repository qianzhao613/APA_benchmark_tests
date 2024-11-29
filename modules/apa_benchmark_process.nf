#!/usr/bin/env nextflow

// set dsl version 
nextflow.enable.dsl=2
// defines process to perform APA tools benchmark
process apa_annotation {

        publishDir("${params.annotation_folder}/${sample_name}", mode: "copy", overwrite: true)

        input:
				path prog_dir
				path gtf_file
				path apa_output_dir
				path barcode_file
				val sample_name
				val apa_method
				val core_num

        output:
				path "${apa_method}_peak_annotations_by_SCAPE.qs", emit: apa_annotation_file

        script:
				"""
				Rscript ${prog_dir}/annotate_site_by_SCAPE.R --method ${apa_method} --nthread ${core_num} --gtf ${gtf_file} --APAdir ${apa_output_dir} --barcodes ${barcode_file}
				"""
}

process nanopore_normalization {

        publishDir("${params.nanopore_normalization_folder}", mode: "copy", overwrite: true)

        input:
				path prog_dir
				val sample_name
				val core_num

        output:
				path "${sample_name}_nanopore_normalized_isomatrix.qs", emit: nanopore_normalization_file

        script:
				"""
				Rscript ${prog_dir}/nanopore_normalization.R --sample ${sample_name} --nthread ${core_num}
				"""
}

process apa_normalization {

        publishDir("${params.apa_normalization_folder}", mode: "copy", overwrite: true)

        input:
				path prog_dir
				path apa_output_dir
				path apa_annotation
				val sample_name
				val apa_method
				val core_num

        output:
				tuple path("${sample_name}_${apa_method}_count_matrix.qs"), path("${sample_name}_${apa_method}_normalized_by_ranger_matrix.qs"), emit: apa_matrix_file

        script:
				"""
				Rscript ${prog_dir}/apa_normalization_v2.R --sample ${sample_name} --method ${apa_method} --APAannotation ${apa_annotation} --APAdir ${apa_output_dir} --nthread ${core_num}
				"""
}

process nanopore_filter {

        publishDir("${params.annotation_folder}/${sample_name}", mode: "copy", overwrite: true)

        input:
				path prog_dir
				path nanopore_normalization_file
				val sample_name
				val core_num
				val cell_type

        output:
				path "${sample_name}_${cell_type}_normalized_nanopore_isoform_table.qs", emit: nanopore_filter_file

        script:
				"""
				Rscript ${prog_dir}/nanopore_site_filter_v2.R --sample ${sample_name} --ground_truth ${nanopore_normalization_file} --nthread ${core_num} --celltype ${cell_type}
				"""
}

process overlap_genes {

        publishDir("${params.annotation_folder}/${sample_name}", mode: "copy", overwrite: true)

        input:
				path prog_dir
				val method_list
				val sample_name
				val core_num
				path annotation_file

        output:
				path "${sample_name}_overlap_gene_list.qs", emit: overlap_genes

        script:
				"""
				Rscript ${prog_dir}/apa_overlap_genes_v1.R --sample ${sample_name} --method_list ${method_list} --APAannotation ${params.annotation_folder}/${sample_name} --nthread ${core_num}
				"""
}

process identification {

        publishDir("${params.identification_folder}", mode: "copy", overwrite: true)

        input:
				path prog_dir
				path apa_output_dir
				val sample_name
				val apa_method
				tuple val(win_size), val(cv)
				path overlap_genes
				val core_num

        output:
				path "${apa_method}_${sample_name}_*_identification_*.qs"

        script:
				"""
				Rscript ${prog_dir}/benchmark_identification_v2.R --sample ${sample_name} --method ${apa_method} --apa_count ${apa_output_dir} --nthread ${core_num} --win ${win_size} --gene_list ${overlap_genes} --cv ${cv}
				"""
}



process quantification {

        publishDir("${params.quantification_results_folder}", mode: "copy", overwrite: true)

        input:
				path prog_dir
				path nanopore_normalization_data
				path apa_normalization_data
				val sample_name
				val apa_method
				tuple val(win_size), val(cv)
				path overlap_genes
				val cell_type
				val core_num

        output:
				path "${apa_method}_${sample_name}_*_quantification_*.qs"
				path "${apa_method}_${sample_name}_*_match_sites_matrix.qs"

        script:
				"""
				Rscript ${prog_dir}/benchmark_quantification_v2.R --sample ${sample_name} --method ${apa_method} --ground_truth ${nanopore_normalization_data} --APA_data ${apa_normalization_data} --nthread ${core_num} --nanocutoff 0 --apacutoff 0 --celltype ${cell_type} --win ${win_size} --gene_list ${overlap_genes} --cv ${cv}
				"""
}

process differential_expression {

        publishDir("${params.differential_results_folder}", mode: "copy", overwrite: true)

        input:
				path prog_dir
				path apa_output_dir
				val sample_name
				val apa_method
				val win_size
				path overlap_genes
				val core_num

        output:
				path "${apa_method}_${sample_name}_*_metrics.qs"
				path "${apa_method}_${sample_name}_*_match_sites_matrix.qs"
        script:
				"""
				Rscript ${prog_dir}/benchmark_differential_expression_isoform_v2.R --sample ${sample_name} --method ${apa_method} --apa_count ${apa_output_dir} --nthread ${core_num} --nanocutoff 0 --apacutoff 0 --win ${win_size} --gene_list ${overlap_genes}
				"""
}
