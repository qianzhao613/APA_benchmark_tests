#!/usr/bin/env Rscript


# load library ------------------------------------------------------------
library(dplyr, warn.conflicts = F)
library(argparse, warn.conflicts = F) # 2.2.2
# import functions --------------------------------------------------------
import::here("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/bin/rscript/benchmark_functions.R", generate_cell_list)

# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Generate nanopore identification benchmark data')

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
parser$add_argument('--nanocutoff', 
                    dest='nanopore_cutoff',
                    action='store',
                    type='integer',
                    default=0,
                    help='the cutoff of nanopore (default: 0)')
parser$add_argument('--celltype', 
                    dest='cell_type',
                    action='store',
                    default="all",
                    help='Users can choose at which level to generate ground-truth data, such as pseudo (all), cell type (ct), or  single cell (sc).')

args <- parser$parse_args()


# inputs ------------------------------------------------------------------

# inputs from parameters
core_num <- args[["nthreads"]]
sample_name <- args[["sample_name"]] #"MOBV1" "E18" MOBV1  args[1] 
nanopore_reads_cutoff <- args[["nanopore_cutoff"]] # 0, 3, 5, 10, 15, 20, 25, 30
cell_type <- args[["cell_type"]] # all/ct/sc

# reference annotation [gene + peakID]
mm10_cr_gtf <- qs::qread("reference/mm10_cr_annot.qs", nthreads = core_num)

# cell metadata [cell_ids + counts_nanopore + counts_without_introns + seurat_clusters]
cell_metadata_file <- glue::glue("/mnt/mr01-home01/m57549qz/scratch/nextflow_test/data/{sample_name}_illumina/{sample_name}_cell_expression_annotated.qs")
cell_metadata <- qs::qread(file = cell_metadata_file, nthreads = core_num) %>% dplyr::mutate(cell_ids = gsub("-1", "", cell_ids)) # REMOVE "-1" if not remove before in case some outputs do not have "-1"

cell_ids <- generate_cell_list(cell_metadata = cell_metadata, type = cell_type, core_num = core_nums)

## read nanopore count matrix
nanopore_matrix <- vroom::vroom(glue::glue("~/nanopore_validation/{sample_name}_illumina/{sample_name}_isomatrix.txt"), show_col_types = F)


## calucalte total
nanopore_total_counts <- lapply(names(cell_ids), function(x){
  # extract cell ids
  selected_cell_ids <- cell_ids[[x]]
  
  # sum up each rows, then filter, and finally add ref annotation
  nanopore_total <- nanopore_matrix %>%
    dplyr::mutate(
      transcriptId = gsub("\\.[0-9]{1,2}", "", transcriptId),
      total = rowSums(pick(all_of(selected_cell_ids)))) %>%
    dplyr::select(geneId, transcriptId, total) %>%
    dplyr::filter(total > nanopore_reads_cutoff) %>%
    dplyr::left_join(mm10_cr_gtf, by = c("transcriptId" = "transcript_id")) %>% 
    dplyr::filter(!is.na(peakID)) # assign peakID, not run the current version
  
  # save
  x <- stringr::str_replace_all(x, "[/\\-+() ]", "_")
  nanopore_output_file <- glue::glue("{sample_name}_{x}_filter_{nanopore_reads_cutoff}_nanopore.qs")
  qs::qsave(nanopore_total, file = nanopore_output_file, nthreads = core_nums)
  
})

# output files ------------------------------------------------------------
# output_dir <- glue::glue("~/nanopore_validation/identification_validation/ground_truth/{sample_name}")
# 
# if(!file.exists(output_dir)){
#   
#   dir.create(output_dir, recursive = T)
#   
# }
