#!/usr/bin/env Rscript

# load library ------------------------------------------------------------
library(dplyr, quietly = T, warn.conflicts = F)
library(argparse, quietly = T, warn.conflicts = F)
# functions ---------------------------------------------------------------

# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Normalize nanopore expression')

parser$add_argument('--sample', 
                    dest='sample_name',
                    action='store',
                    help='the sample name (eg. MOBV1)')

parser$add_argument('--nthreads', 
                    dest='nthreads',
                    action='store',
                    type='integer',
                    default=8,
                    help='the number of threads (default: 8)')

args <- parser$parse_args()

## input load once
core_num <- args[["nthreads"]]
sample_name <- args[["sample_name"]] #"MOBV1" "E18" MOBV1  args[1] 

# data --------------------------------------------------------------------

## cell metadata (load once)
cell_metadata_file <- glue::glue("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/data/{sample_name}_illumina/{sample_name}_cell_expression_annotated.qs") 
cell_ids <- qs::qread(file = cell_metadata_file, nthreads = core_num) %>% 
  dplyr::mutate(cell_ids = gsub("-1", "", cell_ids)) %>% # need to REMOVE "-1" if not remove before
  dplyr::pull(cell_ids)

## ground truth
nanopore_file <-  glue::glue("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/data/{sample_name}_illumina/{sample_name}_isomatrix.txt.gz") 
nanopore_matrix <- vroom::vroom(nanopore_file, show_col_types = F)
scale_factor <- 10000

normalize_nanopore_matrix <- nanopore_matrix %>% 
  dplyr::reframe(across(all_of(cell_ids), function(x){x/sum(x)*scale_factor})) %>% 
  dplyr::mutate(
    geneId = nanopore_matrix[["geneId"]],
    transcriptId = gsub("\\.[0-9]{1,2}", "", nanopore_matrix[["transcriptId"]])
  )

## save
# output_dir <- "~/nanopore_validation/quantification_validation/ground_truth/"
output_file <- glue::glue("{sample_name}_nanopore_normalized_isomatrix.qs")
qs::qsave(normalize_nanopore_matrix, file = output_file, nthreads = core_num)


