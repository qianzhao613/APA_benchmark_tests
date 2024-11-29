#!/usr/bin/env Rscript


# load library ------------------------------------------------------------
library(dplyr)
library(argparse, warn.conflicts = F)
# functions ---------------------------------------------------------------
# entropy function
calculate_entropy <- function(x) {
  
  p <- x / sum(x)
  p <- p[p > 0] # Remove zero probabilities to avoid log(0)
  
  return(-sum(p * log(p)))
}

# uncertainty function
gene_uncertainty_metrics <- function(express_matrix, alpha = 0.05){
  
  sample_size <- ncol(express_matrix)
  gene_mean <- rowMeans(express_matrix)
  gene_sum <- rowSums(express_matrix)
  gene_variance <- apply(express_matrix, 1, var) # optimization based on definition didn't work, probably due to not matrix format or float calculation
  gene_sd <- sqrt(gene_variance) # optimize from apply(matrix, 1, sd)
  gene_entropy <- apply(express_matrix, 1, calculate_entropy)
  gene_sparse <- rowSums(express_matrix == 0) / sample_size # optimize from apply(matrix, 1, function(x) sum(x==0)/length(x))

  gene_table <- data.frame(
    gene_sum = gene_sum,
    gene_mean =  gene_mean, 
    gene_variance = gene_variance, 
    gene_sd = gene_sd,
    sample_size = sample_size,
    gene_sparse = gene_sparse,
    gene_entropy = gene_entropy
  ) %>% 
    dplyr::filter(gene_sum > 0) %>% 
    tibble::rownames_to_column("gene_name") %>% 
    dplyr::mutate(
      lower_bound = gene_mean - (qt(p = alpha/2, df = (sample_size - 1), lower.tail = F) * (gene_sd / sqrt(sample_size))), # mean - marginal error
      upper_bound = gene_mean + (qt(p = alpha/2, df = (sample_size - 1), lower.tail = F) * (gene_sd / sqrt(sample_size))), # mean + marginal error
      cv = gene_sd/gene_mean
    )   # optimize from apply(matrix, 1, function(x) calculate_CI(mean(x), length(x), sd(x), 0.05))
  
  return(gene_table)
}

# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Quantification performance')

parser$add_argument('--sample', 
                    dest='sample_name',
                    action='store',
                    help='the sample name (eg. MOBV1)')
parser$add_argument('--ground_truth', 
                    dest='ground_truth',
                    action='store',
                    default='',
                    help='path to the nanopore normalization file (qs format)')
parser$add_argument('--nthreads', 
                    dest='nthreads',
                    action='store',
                    type='integer',
                    default=8,
                    help='the number of threads (default: 8)')
parser$add_argument('--celltype', 
                    dest='cell_type',
                    action='store',
                    default="all",
                    help='Users can choose at which level to generate ground-truth data, 
                    such as pseudo (all), cell type name.')

args <- parser$parse_args()

# inputs ------------------------------------------------------------------

## input load once

sample_name <- args[["sample_name"]] 
nanopore_normalization_file <- args[["ground_truth"]]
core_num <- args[["nthreads"]]
cell_type <- args[["cell_type"]] # all/ct/sc

normalization = TRUE

# method1 -----------------------------------------------------------------
cell_metadata_file <- glue::glue("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/data/{sample_name}_illumina/{sample_name}_cell_expression_annotated.qs") 
cell_metadata <- qs::qread(file = cell_metadata_file, nthreads = core_num) %>% dplyr::mutate(cell_ids = gsub("-1", "", cell_ids)) # need to REMOVE "-1" if not remove before

if(isTRUE(cell_type == "all")){
  
  cell_ids <- cell_metadata %>% dplyr::pull(cell_ids)
  
} else {
  
  cell_ids <- cell_metadata %>% 
    dplyr::filter(seurat_clusters == cell_type) %>% 
    dplyr::pull(cell_ids) # keep cells in the specific cell type
  
}

# raw or normalization matrix
if(isTRUE(normalization)){
  
  nanopore_matrix <- qs::qread(nanopore_normalization_file, nthreads = core_num) %>% 
    tidyr::unite(col = gene_name, c("geneId", "transcriptId"), sep = ":", remove = T) %>% 
    tibble::column_to_rownames("gene_name") %>% 
    dplyr::select(dplyr::all_of(cell_ids))

} else {
  # raw UMI coutns
  nanopore_file <-  glue::glue("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/data/{sample_name}_illumina/{sample_name}_isomatrix.txt.gz") 
  
  nanopore_matrix <- vroom::vroom(nanopore_file, show_col_types = F) %>% 
    tidyr::unite(col = gene_name, c("geneId", "transcriptId"), sep = ":", remove = T) %>% 
    tibble::column_to_rownames("gene_name") %>% 
    dplyr::select(dplyr::all_of(cell_ids))
}

gene_table <- gene_uncertainty_metrics(nanopore_matrix) %>% 
  dplyr::mutate(sample = sample_name, type = cell_type) %>% 
  tidyr::separate(col = gene_name, into = c("gene_name", "transcript_id"), sep = ":") %>% 
  dplyr::mutate(transcript_id = stringr::str_replace(transcript_id, pattern = "\\.[0-9]{1,}", replacement = ""))# 14s

qs::qsave(gene_table, file = glue::glue("{sample_name}_{cell_type}_normalized_nanopore_isoform_table.qs"), nthreads = core_num)


