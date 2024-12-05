#!/usr/bin/env Rscript

library(dplyr, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts = F))
library(argparse, quietly = T, warn.conflicts = F)

# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Quantification performance')

parser$add_argument('--sample', 
                    dest='sample_name',
                    action='store',
                    help='the sample name (eg. MOBV1)')
parser$add_argument('--method', 
                    dest='method',
                    action='store',
                    help='the APA identification method (eg. SCAPE)')
parser$add_argument('--ground_truth', 
                    dest='ground_truth',
                    action='store',
                    default='',
                    help='path to the nanopore normalization file (qs format)')
parser$add_argument('--APA_data', 
                    dest='APA_data',
                    action='store',
                    default='',
                    help='path to APA normalization data (qs format)')
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
parser$add_argument('--celltype', 
                    dest='cell_type',
                    action='store',
                    default="all",
                    help='Users can choose at which level to generate ground-truth data, such as pseudo (all), cell type (ct), or  single cell (sc).')
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
nanopore_reads_cutoff <- args[["nanopore_cutoff"]] 
apa_reads_cutoff <- args[["apa_cutoff"]]
cell_type <- args[["cell_type"]] # all/ct/sc
nanopore_normalization_file <- args[["ground_truth"]]
apa_normalization_file <- args[["APA_data"]]
win_size <- args[["win"]]
cv_cutoff <- args[["cv"]]
gene_list <- args[["gene_list"]]
pipeline_dir <- args[["pipelinedir"]]

single_tail_compare <- FALSE
multi_tail_compare <- FALSE #TRUE
search_top_gene <- FALSE #TRUE

apa_top_gene_num <- 0 # 3000, 6000, 9000

mm10_cr_gtf <- qs::qread(glue::glue("{pipeline_dir}/reference/mm10_cr_annot.qs"), nthreads = core_num)

# match_matrix_output_dir <- glue::glue("{output_dir}/match_sites_expression")
# 
# if(dir.exists(match_matrix_output_dir)){print("Dir exists!")} else {dir.create(match_matrix_output_dir, recursive = T)}

# functions ---------------------------------------------------------------
import::here(glue::glue("{pipeline_dir}/bin/rscript/benchmark_functions.R"), 
             apa_matrix_extract,
             generate_cell_list,
             peak2bed, 
             bed2gr,
             matched_sites_expression)

# cell metadata ---------------------------------------------------------

cell_metadata_file <- glue::glue("{pipeline_dir}/data/{sample_name}_illumina/{sample_name}_cell_expression_annotated.qs") 

cell_metadata <- qs::qread(file = cell_metadata_file, nthreads = core_num) %>% dplyr::mutate(cell_ids = gsub("-1", "", cell_ids)) # need to REMOVE "-1" if not remove before
cell_ids <- generate_cell_list(cell_metadata = cell_metadata, type = cell_type, core_num = core_num)
cell_ids <- cell_ids[[cell_type]]

# nanopore ----------------------------------------------------------------
## ground truth

# nanopore_normalization_file <- glue::glue("{quantification_dir}/ground_truth/{sample_name}_nanopore_normalized_isomatrix.qs")
normalize_nanopore_matrix <- qs::qread(nanopore_normalization_file, nthreads = core_num)

# sum up each rows, then filter, and finally add ref annotation
nanopore_total <- normalize_nanopore_matrix %>%
  dplyr::mutate(total = rowSums(pick(all_of(cell_ids)))) %>%
  dplyr::select(geneId, transcriptId, total) %>%
  dplyr::filter(total > nanopore_reads_cutoff) %>%
  dplyr::left_join(mm10_cr_gtf, by = c("transcriptId" = "transcript_id")) %>% 
  dplyr::filter(!is.na(peakID)) 

if(isTRUE(multi_tail_compare)){
  
  multi_tail_genes <- nanopore_total %>% 
    dplyr::group_by(geneId) %>% 
    dplyr::summarise(n=n()) %>% 
    dplyr::filter(n >= 2) %>% 
    dplyr::pull(geneId)
  
  nanopore_total <- nanopore_total %>% dplyr::filter(geneId %in% multi_tail_genes)
  
}

if(isTRUE(single_tail_compare)){
  
  single_tail_genes <- nanopore_total %>% 
    dplyr::group_by(geneId) %>% 
    dplyr::summarise(n=n()) %>% 
    dplyr::filter(n == 1) %>% 
    dplyr::pull(geneId)
  
  nanopore_total <- nanopore_total %>% dplyr::filter(geneId %in% single_tail_genes)
  
}

isoform_table_file <- glue::glue("{pipeline_dir}/output/apa_annotation/{sample_name}/{sample_name}_{cell_type}_normalized_nanopore_isoform_table.qs")

if(cv_cutoff == "none"){
  
  print("No CV filter for ONT data!")
  
} else {
  cv_cutoff <- as.numeric(cv_cutoff)
  transcript_ids <- qs::qread(file = isoform_table_file, nthreads = core_num) %>%
    dplyr::filter(cv < cv_cutoff) %>%
    dplyr::pull(transcript_id) %>%
    unique()
  
  nanopore_total <- nanopore_total %>% dplyr::filter(transcriptId %in% transcript_ids)
  print(glue::glue("CV filer {cv_cutoff} apply on ONT data!"))
}

# read APA data -----------------------------------------------------------
## apa expression matrix
normalize_apa_matrix <- qs::qread(apa_normalization_file, nthreads = core_num)

apa_total <- normalize_apa_matrix %>% 
  dplyr::mutate(total = rowSums(pick(all_of(cell_ids)))) %>% 
  dplyr::select(peakID, ensemblId, geneId, transcriptId, total) %>% 
  dplyr::distinct() %>%  # Sierra gene_name inconsistent with SCAPE annotation, so sometimes duplicate
  dplyr::filter(total > apa_reads_cutoff) #27737

### if need to use top genes instead of filtering
if(isTRUE(search_top_gene)){
  
  apa_gene_name <- apa_total %>% 
    dplyr::filter(!is.na(geneId)) %>% 
    dplyr::group_by(geneId) %>% 
    dplyr::summarise(across(total, sum)) %>% 
    dplyr::arrange(dplyr::desc(total)) %>% 
    dplyr::slice(1:apa_top_gene_num) %>% 
    dplyr::pull(geneId)
  
} else {
  
  apa_gene_name <- unique(apa_total$geneId)
  
}

# overlap genes -----------------------------------------------------------
## first, keep overlapped genes
### extract overlap gene list
# intersect_gene_name <- intersect(unique(nanopore_total$geneId), apa_gene_name) # 9409
gene_list <- qs::qread(gene_list, nthreads = core_num)
intersect_gene_name <- intersect(unique(nanopore_total$geneId), unique(gene_list)) 

### keep overlapping genes in matrix
nanopore_total <- nanopore_total %>% 
  dplyr::filter(geneId %in% intersect_gene_name) %>%  # 25080
  dplyr::filter(!is.na(peakID)) # assign peakID

apa_total <- apa_total %>% filter(geneId %in% intersect_gene_name) # 24777

## second: overlap sites (true positive) between ground truth sites and predicted sites
# win_size_vector <- seq(from = 0, to = 300, by = 10)

# get the number of overlap sites at the gene level and add counts data

#get matched site_level expression data
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

qs::qsave(matched_expression, 
          file = glue::glue("{apa_method}_{sample_name}_filter_{apa_reads_cutoff}_nanopore_{nanopore_reads_cutoff}_win{win_size}_cv{cv_cutoff}_match_sites_matrix.qs"),
          nthreads = core_num)
# #calculate r2
# r_squared <- summary(lm(nanopore_score ~ apa_score, data = matched_expression))$r.squared
# log_r_squared <- summary(lm(log(nanopore_score, base=2) ~ log(apa_score, base=2), data = matched_expression))$r.squared  # 0.4048
# return(c(win_size, r_squared, log_r_squared))

# calculate r
r_pearson <- cor.test(log(matched_expression$nanopore_score, base=2), log(matched_expression$apa_score, base=2), method = "pearson")$estimate
r_spearman_all_gene <- cor.test(matched_expression$nanopore_score, matched_expression$apa_score, method = "spearman")$estimate
r_spearman_per_gene <- cor.test(matched_expression$nanopore_rank, matched_expression$apa_rank, method = "pearson")$estimate


quantification_metrics <- c(
  win_size = win_size, 
  sample = sample_name,
  method = apa_method,
  r_pearson = r_pearson, 
  r_spearman_all_gene = r_spearman_all_gene, 
  r_spearman_per_gene = r_spearman_per_gene, 
  n_overlap_site = nrow(matched_expression), 
  n_overlap_gene = length(unique(matched_expression$geneId)),
  n_apa_sites = nrow(apa_total), 
  n_nanopore_sites = nrow(nanopore_total),
  overlap_gene = length(intersect_gene_name)
  )
# # final metrics table
# quantification_metrics <- quantification_metrics %>%
#   purrr::reduce(rbind) %>%
#   dplyr::as_tibble(.name_repair = 'unique') %>%
#   setNames(c("win_size", "r_pearson", "r_spearman_all_gene", "r_spearman_per_gene", "n_overlap_site", "n_overlap_gene")) %>%
#   dplyr::mutate(method = apa_method, overlap_gene = length(intersect_gene_name), n_apa_sites = nrow(apa_total), n_nanopore_sites = nrow(nanopore_total))

# Initialize the base filename
quantification_output_filename <- glue::glue("{apa_method}_{sample_name}_")

# Add specific components based on the conditions
if (isTRUE(search_top_gene)) {
  quantification_output_filename <- glue::glue(
    "{quantification_output_filename}top_{apa_top_gene_num}_nanopore_{nanopore_reads_cutoff}_quantification_win{win_size}_cv{cv_cutoff}_metrics.qs"
  )
  
} else {
  # Initialize base with common components
  quantification_output_filename <- glue::glue(
    "{quantification_output_filename}filter_{apa_reads_cutoff}_nanopore_{nanopore_reads_cutoff}_quantification_win{win_size}_cv{cv_cutoff}_metrics"
  )
  
  # Add specific tail components
  if (isTRUE(multi_tail_compare)) {
    quantification_output_filename <- glue::glue("{quantification_output_filename}_multi_tails.qs")
  } else if (isTRUE(single_tail_compare)) {
    quantification_output_filename <- glue::glue("{quantification_output_filename}_single_tail.qs")
  } else {
    # Add the default .qs extension if no tail comparison is true
    quantification_output_filename <- glue::glue("{quantification_output_filename}.qs")
  }
}

qs::qsave(quantification_metrics, 
          file = quantification_output_filename,
          nthreads = core_num)
