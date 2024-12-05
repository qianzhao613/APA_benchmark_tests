#!/usr/bin/env Rscript

# load library ------------------------------------------------------------
library(dplyr, quietly = T, warn.conflicts = F)
library(argparse, quietly = T, warn.conflicts = F)

# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Normalize APA matrix')

parser$add_argument('--sample', 
                    dest='sample_name',
                    action='store',
                    help='the sample name (eg. MOBV1)')
parser$add_argument('--method', 
                    dest='method',
                    action='store',
                    help='the APA identification method (eg. SCAPE)')
parser$add_argument('--APAannotation', 
                    dest='APAannotation',
                    action='store',
                    default='',
                    help='path to the APA annotation file (qs formart)')
parser$add_argument('--APAdir', 
                    dest='APAdir',
                    action='store',
                    default='',
                    help='path to APA output')
parser$add_argument('--nthreads', 
                    dest='nthreads',
                    action='store',
                    type='integer',
                    default=8,
                    help='the number of threads (default: 8)')
parser$add_argument('--pipelinedir', 
                    dest='pipelinedir',
                    action='store',
                    default='.',
                    help='the pipeline folder')
					
args <- parser$parse_args()

## input load once
core_num <- args[["nthreads"]]
sample_name <- args[["sample_name"]] 
apa_method <- args[["method"]] 
apa_mat_file <- args[["APAdir"]]
annot_file <- args[["APAannotation"]]
pipeline_dir <- args[["pipelinedir"]]

scale_factor <- 10000

# function ----------------------------------------------------------------
import::here(glue::glue("{pipeline_dir}/bin/rscript/benchmark_functions.R"), 
             generate_cell_list, apa_matrix_extract)
			 
# data --------------------------------------------------------------------

## cell metadata (load once)
cell_metadata_file <- glue::glue("{pipeline_dir}/data/{sample_name}_illumina/{sample_name}_cell_expression_annotated.qs")

cell_ids <- qs::qread(file = cell_metadata_file, nthreads = core_num) %>% 
  dplyr::mutate(cell_ids = gsub("-1", "", cell_ids)) %>% # need to REMOVE "-1" if not remove before
  dplyr::pull(cell_ids)

## short reads
illumina_matrix <- glue::glue("{pipeline_dir}/data/{sample_name}_illumina/{sample_name}_illumina_counts_sparseMartix.qs")
illumina_matrix <- qs::qread(illumina_matrix, nthreads = core_num)

colnames(illumina_matrix) <- gsub("-1", "", colnames(illumina_matrix))
illumina_matrix <- illumina_matrix[, cell_ids]


## apa expression matrix
feature_list <- c("3UTRs", "Exon", "Intron", "LastExon1Kb", "3UTRs_1kb", "3UTRs_2kb") # this is based on SCAPE annotation strategy
apa_annot <- qs::qread(annot_file, nthreads =  core_num) %>% dplyr::filter(Type %in% feature_list) 

apa_matrix <- apa_matrix_extract(expr_file = apa_mat_file, annot = apa_annot, method = apa_method, barcodes = cell_ids)

## save raw counts
output_file1 <- glue::glue("{sample_name}_{apa_method}_count_matrix.qs")
qs::qsave(apa_matrix, file = output_file1, nthreads = core_num)

# 2nd strategy
apa_matrix_annot <- apa_matrix %>% 
  dplyr::distinct(peakID, annot.symbol, .keep_all = TRUE) %>% 
  dplyr::rename(transcriptId = annot.tx_name, geneId = annot.symbol, ensemblId = annot.gene_id) %>% 
  dplyr::select(peakID, ensemblId, geneId, transcriptId) # 15032 -- 15027 peakID, ensemblId, geneId, transcriptId

apa_matrix <- apa_matrix %>%
  dplyr::distinct(peakID, annot.symbol, !!!syms(cell_ids), .keep_all = TRUE) %>% 
  dplyr::group_by(peakID) %>% 
  dplyr::summarise(across(all_of(cell_ids), sum)) %>% 
  tibble::column_to_rownames("peakID") 

apa_matrix_rownames <- rownames(apa_matrix)

apa_matrix <- apa_matrix %>% 
  Matrix::as.matrix() %>% 
  as(., "dgCMatrix")

if(isTRUE(all.equal(colnames(apa_matrix), colnames(illumina_matrix)))){

  col_sums <- Matrix::colSums(illumina_matrix)
  diag_col_sums_inv <- Matrix::Diagonal(x = 1 / col_sums)
  normalize_apa_matrix <- apa_matrix %*% diag_col_sums_inv
  normalize_apa_matrix <- normalize_apa_matrix * scale_factor
  colnames(normalize_apa_matrix) <- colnames(apa_matrix)
  
}

normalized_matrix_rowname <- rownames(normalize_apa_matrix)

normalize_apa_matrix <- normalize_apa_matrix %>% 
  as.matrix() %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(peakID = normalized_matrix_rowname) %>% 
  dplyr::left_join(apa_matrix_annot, by = "peakID")

## save normalized counts
output_file2 <- glue::glue("{sample_name}_{apa_method}_normalized_by_ranger_matrix.qs")
qs::qsave(normalize_apa_matrix, file = output_file2, nthreads = core_num)
