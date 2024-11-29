#!/usr/bin/env Rscript

library(dplyr, warn.conflicts = F)
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts = F))
library(argparse, warn.conflicts = F)

# functions ---------------------------------------------------------------
import::here("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/bin/rscript/benchmark_functions.R", 
             generate_cell_list,
             peak2bed, 
             bed2gr, 
             find_overlap_sites_per_gene)

# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Identification performance')

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

args <- parser$parse_args()

## input load once
core_num <- args[["nthreads"]]
sample_name <- args[["sample_name"]] 
apa_method <- args[["method"]] 
apa_count_file <- args[["apa_count_file"]]

nanopore_reads_cutoff <- args[["nanopore_cutoff"]] 
apa_reads_cutoff <- args[["apa_cutoff"]]
cell_type <- args[["cell_type"]] # all/ct/sc
win_size <- args[["win"]]
cv_cutoff <- args[["cv"]]
gene_list <- args[["gene_list"]]

single_tail_compare <- FALSE
multi_tail_compare <- FALSE 
search_top_gene <- FALSE 
apa_top_gene_num <- 0 

# pre-load data -----------------------------------------------------------
## cell metadata (load once)
cell_metadata_file <- glue::glue("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/data/{sample_name}_illumina/{sample_name}_cell_expression_annotated.qs") #"E18_cell_expression_annotated.qs"
cell_metadata <- qs::qread(file = cell_metadata_file, nthreads = core_num) %>% dplyr::mutate(cell_ids = gsub("-1", "", cell_ids)) # need to REMOVE "-1" if not remove before
cell_ids <- generate_cell_list(cell_metadata = cell_metadata, type = cell_type, core_num = core_num)
cell_ids <- cell_ids[[cell_type]]

## ground truth
mm10_cr_gtf <- qs::qread("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/reference/mm10_cr_annot.qs", nthreads = core_num)

nanopore_total <- glue::glue("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/data/{sample_name}_illumina/{sample_name}_isomatrix.txt.gz") %>% 
  vroom::vroom(., show_col_types = F) %>%
  dplyr::mutate(
    transcriptId = gsub("\\.[0-9]{1,2}", "", transcriptId),
    total = rowSums(pick(all_of(cell_ids)))) %>%
  dplyr::select(geneId, transcriptId, total) %>%
  dplyr::filter(total > nanopore_reads_cutoff) %>%
  dplyr::left_join(mm10_cr_gtf, by = c("transcriptId" = "transcript_id")) %>% 
  dplyr::filter(!is.na(peakID)) # assign peakID, not run the current version

if(isTRUE(multi_tail_compare)){
  
  multi_tail_genes <- nanopore_total %>% 
    dplyr::group_by(geneId) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::filter(n >= 2) %>% 
    pull(geneId)
  
  nanopore_total <- nanopore_total %>% 
    dplyr::filter(geneId %in% multi_tail_genes)
  
}

if(isTRUE(single_tail_compare)){
  
  single_tail_genes <- nanopore_total %>% 
    dplyr::group_by(geneId) %>% 
    dplyr::summarise(n=n()) %>% 
    dplyr::filter(n == 1) %>% pull(geneId)
  
  nanopore_total <- nanopore_total %>% 
    dplyr::filter(geneId %in% single_tail_genes)
  
}


if(cv_cutoff == "none"){
  
  print("No CV filter for ONT data!")
  
} else {
  isoform_table_file <- glue::glue("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/APA_benchmark_tests/output/apa_annotation/{sample_name}/{sample_name}_{cell_type}_normalized_nanopore_isoform_table.qs")

  cv_cutoff <- as.numeric(cv_cutoff)
  transcript_ids <- qs::qread(file = isoform_table_file, nthreads = core_num) %>%
    dplyr::filter(cv < cv_cutoff) %>%
    dplyr::pull(transcript_id) %>%
    unique()

  nanopore_total <- nanopore_total %>% dplyr::filter(transcriptId %in% transcript_ids)
  print(glue::glue("CV filer {cv_cutoff} apply on ONT data!"))
}

# read APA data -----------------------------------------------------------
## input change every time
# apa_dir <- glue::glue("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/data/{sample_name}_illumina") 

### apa expression matrix

counts <- qs::qread(apa_count_file, nthreads = core_num)

apa_total <- counts %>% 
  dplyr::mutate(total = rowSums(pick(all_of(cell_ids)))) %>% 
  dplyr::rename(transcriptId = annot.tx_name, geneId = annot.symbol, ensemblId= annot.gene_id) %>% 
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
#intersect_gene_name <- intersect(unique(nanopore_total$geneId), apa_gene_name) # 9409
gene_list <- qs::qread(gene_list, nthreads = core_num)
intersect_gene_name <- intersect(unique(nanopore_total$geneId), unique(gene_list)) 

### keep overlapping genes in matrix
nanopore_total <- nanopore_total %>% 
  dplyr::filter(geneId %in% intersect_gene_name) %>%  # 25080
  dplyr::filter(!is.na(peakID)) # assign peakID

apa_total <- apa_total %>% dplyr::filter(geneId %in% intersect_gene_name) # 24777

## second: overlap sites (true positive) between ground truth sites and predicted sites
# win_size_vector <- seq(from = 0, to = 300, by = 10)

# identify overlap sites
overlap_sum <- parallel::mclapply(
  intersect_gene_name, 
  find_overlap_sites_per_gene, # call the function
  apa_data = apa_total,
  nanopore_data = nanopore_total,
  win_size = win_size,
  mc.cores = core_num
) %>% 
  purrr::reduce(`+`) %>% # sum up the number of sites
  setNames(c("true_pos", "predicted_sites", "groundtruth_sites")) %>% 
  append(c(
    win_size = win_size,
    sample = sample_name,
    method = apa_method, 
    overlap_gene = length(intersect_gene_name), 
    n_apa_sites = nrow(apa_total), 
    n_nanopore_sites = nrow(nanopore_total)
  ))

# # calculate metrics and return [true_pos/pred = precision; true_pos/ground = recall]
# precision <- overlap_sum[['true_pos']]/overlap_sum[['predicted_sites']]
# recall <- overlap_sum[['true_pos']]/overlap_sum[['groundtruth_sites']]
# 
# # final metrics table
# identification_metrics <- identification_metrics %>% 
#   purrr::reduce(rbind) %>% 
#   dplyr::as_tibble(.name_repair = "unique_quiet") %>% 
#   setNames(c("win_size", "precision", "recall")) %>% 
#   dplyr::mutate(method = apa_method, overlap_gene = length(intersect_gene_name), n_apa_sites = nrow(apa_total), n_nanopore_sites = nrow(nanopore_total))


# Initialize the base filename
identification_output_filename <- glue::glue("{apa_method}_{sample_name}_")

# Add specific components based on the conditions
if (isTRUE(search_top_gene)) {
  identification_output_filename <- glue::glue(
    "{identification_output_filename}top_{apa_top_gene_num}_nanopore_{nanopore_reads_cutoff}_identification_win{win_size}_cv{cv_cutoff}_metrics.qs"
  )
} else {
  # Initialize base with common components
  identification_output_filename <- glue::glue(
    "{identification_output_filename}filter_{apa_reads_cutoff}_nanopore_{nanopore_reads_cutoff}_identification_win{win_size}_cv{cv_cutoff}_metrics"
  )
  
  # Add specific tail components
  if (isTRUE(multi_tail_compare)) {
    identification_output_filename <- glue::glue("{identification_output_filename}_multi_tails.qs")
  } else if (isTRUE(single_tail_compare)) {
    identification_output_filename <- glue::glue("{identification_output_filename}_single_tail.qs")
  } else {
    # Add the default .qs extension if no tail comparison is true
    identification_output_filename <- glue::glue("{identification_output_filename}.qs")
  }
}


qs::qsave(overlap_sum, 
          file = identification_output_filename,
          nthreads = core_num)

