#!/usr/bin/env Rscript
library(dplyr, warn.conflicts = F)
library(DESeq2, warn.conflicts = F)
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts = F))
library(argparse, quietly = T, warn.conflicts = F)

# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Differential APA')

parser$add_argument('--sample', 
                    dest='sample_name',
                    action='store',
                    help='the sample name (eg. MOBV1)')
parser$add_argument('--method', 
                    dest='method',
                    action='store',
                    help='the APA identification method (eg. SCAPE)')
parser$add_argument('--apa_count', 
                    dest='apa_count_file',
                    action='store',
                    default='',
                    help='path to apa output')
parser$add_argument('--nthreads', 
                    dest='nthreads',
                    action='store',
                    type='integer',
                    default=8,
                    help='the number of threads (default: 8)')
parser$add_argument('--nanocutoff', 
                    dest='nanopore_cutoff',
                    action='store',
                    type='integer',
                    default=0,
                    help='the cutoff of nanopore (default: 0)')
parser$add_argument('--apacutoff', 
                    dest='apa_cutoff',
                    action='store',
                    type='integer',
                    default=0,
                    help='the cutoff of APA method (default: 0)')
parser$add_argument('--win', 
                    dest='win',
                    action='store',
                    type='integer',
                    default=0,
                    help='Users can choose window size.')
parser$add_argument('--gene_list', 
                    dest='gene_list',
                    action='store',
                    default="",
                    help='Users can choose gene list.')
parser$add_argument('--cv', 
                    dest='cv',
                    action='store',
                    default="none",
                    help='Users can choose CV cutoff.')
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
apa_count_file <- args[["apa_count_file"]]
nanopore_reads_cutoff <- args[["nanopore_cutoff"]] 
apa_reads_cutoff <- args[["apa_cutoff"]]
win_size <- args[["win"]]
gene_list <- args[["gene_list"]]
cv_cutoff <- args[["cv"]]
pipeline_dir <- args[["pipelinedir"]]

padj_threshold <- 0.05
single_tail_compare <- FALSE
multi_tail_compare <- FALSE 
search_top_gene <- FALSE 
apa_top_gene_num <- 0 

# functions ---------------------------------------------------------------
import::here(glue::glue("{pipeline_dir}/bin/rscript/benchmark_functions.R"),
             "create_pseudobulk", 
			 "two_cluster_deseq2", 
			 "bed2gr", 
			 "peak2bed", 
			 "matched_sites_expression", 
			 .character_only = TRUE)
			 
# annotation and metadata -------------------------------------------------
mm10_cr_gtf <- qs::qread(glue::glue("{pipeline_dir}/reference/mm10_cr_annot.qs"), nthreads = core_num)

## cell metadata
cell_metadata_file <- glue::glue("{pipeline_dir}/data/{sample_name}_illumina/{sample_name}_cell_expression_annotated.qs") 
cell_metadata <- qs::qread(file = cell_metadata_file, nthreads = core_num) %>% dplyr::mutate(cell_ids = gsub("-1", "", cell_ids)) # need to REMOVE "-1" if not remove before

high_percent_types <- cell_metadata %>% dplyr::group_by(seurat_clusters) %>% dplyr::summarise(n = dplyr::n()) %>% dplyr::arrange(n) %>% dplyr::slice_tail(n=2)
cell_type1 <- high_percent_types[["seurat_clusters"]][1] %>% as.character()
cell_type2 <- high_percent_types[["seurat_clusters"]][2] %>% as.character()


gene_list <- qs::qread(gene_list, nthreads = core_num)

# nanopore ----------------------------------------------------------------
## read nanopore data, convert into ground truth polyAsite matrix (time comsumming probably make an independent step) and filter all zero rows
nanopore_file <-  glue::glue("{pipeline_dir}/data/{sample_name}_illumina/{sample_name}_isomatrix.txt.gz")
nanopore_matrix <- vroom::vroom(nanopore_file, show_col_types = F) %>%
  dplyr::filter(!is.na(transcriptId)) %>%
  dplyr::mutate(transcriptId = gsub("\\.[0-9]{1,2}", "", transcriptId)) %>%
  dplyr::left_join(mm10_cr_gtf, by = c("transcriptId" = "transcript_id")) %>%
  dplyr::filter(!is.na(peakID)) %>%
  dplyr::filter(geneId %in% unique(gene_list)) %>%
  tidyr::unite(col = "peakID", c("geneId", "peakID"), sep = ":") %>% 
  dplyr::group_by(peakID) %>%
  dplyr::summarise(across(all_of(cell_metadata[["cell_ids"]]), sum)) %>%
  tibble::column_to_rownames("peakID")



if(cv_cutoff == "none"){
  
  print("No CV filter for ONT data!")
  
} else {
  isoform_table_file <- glue::glue("{pipeline_dir}/output/apa_annotation/{sample_name}/{sample_name}_{cell_type}_normalized_nanopore_isoform_table.qs")
  cv_cutoff <- as.numeric(cv_cutoff)
  transcript_ids <- qs::qread(file = isoform_table_file, nthreads = core_num) %>%
    dplyr::filter(cv < cv_cutoff) %>%
    dplyr::pull(transcript_id) %>%
    unique()

  nanopore_total <- nanopore_total %>% dplyr::filter(transcriptId %in% transcript_ids)
  print(glue::glue("CV filer {cv_cutoff} apply on ONT data!"))
}

# APA tools ---------------------------------------------------------------
## import apa data
apa_matrix <- qs::qread(apa_count_file, nthreads = core_num) %>% 
  dplyr::filter(annot.symbol %in% unique(gene_list)) %>%
  tidyr::unite(col = "peakID", c("annot.symbol", "peakID"), sep = ":") %>% 
  dplyr::select(all_of(c("peakID", cell_metadata[["cell_ids"]]))) %>% 
  dplyr::group_by(peakID) %>% 
  dplyr::summarise(dplyr::across(dplyr::all_of(cell_metadata[["cell_ids"]]), sum)) %>% 
  tibble::column_to_rownames("peakID")

# compare DE --------------------------------------------------------------


## find differential expressed sites by DESeq2 [pseudobulk]
suppressMessages({
  apa_de <- two_cluster_deseq2(apa_matrix, c(cell_type1, cell_type2), cell_metadata, count_filter = 0)
  nano_de <- two_cluster_deseq2(nanopore_matrix, c(cell_type1, cell_type2), cell_metadata, count_filter = 0)
})

rm(apa_matrix, nanopore_matrix, mm10_cr_gtf)
gc()

apa_total <- apa_de %>% dplyr::filter(padj < padj_threshold) %>% 
  tidyr::separate(col = "peakID", into = c("geneId", "seqnames", "sites", "strand"), sep = ":") %>% 
  tidyr::unite(col = "peakID", c("seqnames", "sites", "strand"), sep = ":") %>% 
  dplyr::mutate(total = log2FoldChange, transcriptId = peakID) 

nanopore_total <- nano_de %>% dplyr::filter(padj < padj_threshold) %>% 
  tidyr::separate(col = "peakID", into = c("geneId", "seqnames", "sites", "strand"), sep = ":") %>% 
  tidyr::unite(col = "peakID", c("seqnames", "sites", "strand"), sep = ":") %>% 
  dplyr::mutate(total = log2FoldChange, transcriptId = peakID) 

intersect_gene_name <- intersect(apa_total$geneId, nanopore_total$geneId)

## get matched site_level expression data
matched_expression <- parallel::mclapply(
  intersect_gene_name,
  matched_sites_expression,
  apa_data = apa_total,
  nanopore_data = nanopore_total,
  win_size =  win_size,
  mc.cores = core_num)

matched_expression <- matched_expression %>% 
  purrr::keep(purrr::negate(is.null)) %>% 
  purrr::reduce(rbind) %>% 
  dplyr::distinct() %>%
  dplyr::group_by(nanopore_sites) %>% # if one ground truth matches two or more apa sits
  dplyr::filter(apa_score == max(apa_score)) %>% # only keep the apa site with maxmium expression
  dplyr::ungroup() %>% 
  dplyr::group_by(apa_sites) %>% # if the rest sites contain the situation that two or more ground truth match one apa sites
  dplyr::filter(nanopore_score == max(nanopore_score)) # only keep the nanopore site with maxmium expression

# calculate r
r_pearson <- tryCatch({
  cor.test(matched_expression$nanopore_score, matched_expression$apa_score, method = "pearson")$estimate
}, error = function(e) {
  NA  # in case not enough finite observations
})

#r_spearman_all_gene <- cor.test(matched_expression$nanopore_score, matched_expression$apa_score, method = "spearman")$estimate
#r_spearman_per_gene <- cor.test(matched_expression$nanopore_rank, matched_expression$apa_rank, method = "pearson")$estimate

matched_expression <- matched_expression %>% tidyr::separate(col = "apa_sites", into = c("geneId", "apa_sites"), sep = "\\|")

differential_metrics <- c(
  win_size = win_size, 
  sample = sample_name,
  method = apa_method,
  true_pos = nrow(matched_expression), 
  predicted_site = nrow(apa_total), 
  groundtruth_sites = nrow(nanopore_total),
  precision = nrow(matched_expression)/nrow(apa_total),
  recall = nrow(matched_expression)/nrow(nanopore_total),
  r_pearson = r_pearson, 
  #r_spearman_all_gene = r_spearman_all_gene, 
  #r_spearman_per_gene = r_spearman_per_gene, 
  n_total_gene = length(intersect_gene_name)
)

# Initialize the base filename
differential_output_filename <- glue::glue("{apa_method}_{sample_name}_")

# Add specific components based on the conditions
if (isTRUE(search_top_gene)) {
  differential_output_filename <- glue::glue(
    "{differential_output_filename}top_{apa_top_gene_num}_nanopore_{nanopore_reads_cutoff}_differential_win{win_size}_cv{cv_cutoff}_metrics.qs"
  )
} else {
  # Initialize base with common components
  differential_output_filename <- glue::glue(
    "{differential_output_filename}filter_{apa_reads_cutoff}_nanopore_{nanopore_reads_cutoff}_differential_win{win_size}_cv{cv_cutoff}_metrics"
  )
  
  # Add specific tail components
  if (isTRUE(multi_tail_compare)) {
    differential_output_filename <- glue::glue("{differential_output_filename}_multi_tails.qs")
  } else if (isTRUE(single_tail_compare)) {
    differential_output_filename <- glue::glue("{differential_output_filename}_single_tail.qs")
  } else {
    # Add the default .qs extension if no tail comparison is true
    differential_output_filename <- glue::glue("{differential_output_filename}.qs")
  }
}

qs::qsave(matched_expression,
          file = gsub(".qs", "_match_sites_matrix.qs", differential_output_filename),
          nthreads = core_num)

qs::qsave(differential_metrics, 
          file = differential_output_filename,
          nthreads = core_num)
