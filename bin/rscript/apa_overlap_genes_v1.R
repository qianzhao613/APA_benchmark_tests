#!/usr/bin/env Rscript


# load library ------------------------------------------------------------
library(dplyr)
library(argparse, warn.conflicts = F)
# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Gene overlap')

parser$add_argument('--sample', 
                    dest='sample_name',
                    action='store',
                    help='the sample name (eg. MOBV1)')
parser$add_argument('--method_list', 
                    dest='method_list',
                    nargs='+',
                    default='SCAPE',
                    action='store',
                    help='the APA identification methods (eg. SCAPE)')
parser$add_argument('--APAannotation', 
                    dest='APAannotation',
                    action='store',
                    default='',
                    help='path to the APA annotation file (qs formart)')
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
apa_method_vector <- args[["method_list"]] 
annot_file <- args[["APAannotation"]]
pipeline_dir <- args[["pipelinedir"]]

genes_list <- lapply(apa_method_vector, function(apa_method){
  
  annotation_file <- glue::glue("{annot_file}/{apa_method}_peak_annotations_by_SCAPE.qs")
  
  genes <- qs::qread(annotation_file, nthreads = core_num) %>% dplyr::pull(annot.symbol)
  
  return(genes)
})

intersect_gene_name <- Reduce(intersect, genes_list)

qs::qsave(intersect_gene_name, file = glue::glue("{sample_name}_overlap_gene_list.qs"), nthreads = core_num)
