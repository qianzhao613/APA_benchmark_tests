#!/usr/bin/env Rscript

library(dplyr, warn.conflicts = F)
library(argparse, warn.conflicts = F)

# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Identification performance')

parser$add_argument('--method', 
                    dest='method',
                    action='store',
                    help='the APA identification method (eg. SCAPE)')
parser$add_argument('--species', 
                    dest='species',
                    action='store',
                    default='Mm10',
                    help='the species used for annotation (eg. Mm10)')
parser$add_argument('--nthreads', 
                    dest='nthreads',
                    action='store',
                    type='integer',
                    default=8,
                    help='the number of threads (default: 8)')
parser$add_argument('--gtf', 
                    dest='gtf',
                    action='store',
                    default='genes.gtf',
                    help='the genome annotation file (gtf formart)')
parser$add_argument('--APAdir', 
                    dest='APAdir',
                    action='store',
                    default='illumina',
                    help='path to apa output')
parser$add_argument('--barcodes', 
                    dest='barcodes',
                    action='store',
                    default='genes.gtf',
                    help='the genome annotation file (gtf formart)')
parser$add_argument('--pipelinedir', 
                    dest='pipelinedir',
                    action='store',
                    default='.',
                    help='the pipeline folder')
args <- parser$parse_args()

## input load once
apa_method <- args[["method"]] 
species_name <- args[["species"]] 
gtf_file <- args[["gtf"]] 
apa_output_dir <- args[["APAdir"]] 
barcode_file <- args[["barcodes"]] 
core_num <- args[["nthreads"]]
pipeline_dir <- args[["pipelinedir"]]
# functions ---------------------------------------------------------------
import::here(glue::glue("{pipeline_dir}/bin/rscript/benchmark_functions.R"), 
             apa_annotation)
			 
# scape annotation -------------------------------------------------------------------
# # It will consume a lot of time if it is the first time to annotate.

apa_annot <- apa_annotation(
  apa_dir = apa_output_dir, 
  method = apa_method, 
  gtf_file = gtf_file, 
  barcode_file = barcode_file,
  species = species_name, 
  core_num = core_num)

qs::qsave(apa_annot, file = glue::glue("{apa_method}_peak_annotations_by_SCAPE.qs"), nthreads = 8)




