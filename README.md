# APA benchmark tests
This is a benchmark pipeline based on Nextflow, designed to evaluate the performance of APA (Alternative Polyadenylation) analysis tools for single-cell transcriptomics (or spatial transcriptomics) based on 3' sequencing.

## Overview

![pipeline overview](readme_figures/Figure1B_renew.png)


The pipeline is divided into two main parts:

1. **Data Preparation**: This part reads the outputs from APA analysis tools and ONT (Oxford Nanopore Technologies), converts them into a standard format, and normalizes the data.

2. **Performance Evaluation**: This part consists of three modules:
   - **Site Identification Comparison**: Compares the identified polyadenylation sites.
   - **Site Quantification Comparison**: Compares the quantification of polyadenylation sites.
   - **Differential Expression Analysis Comparison**: Compares the differential expression analysis of polyadenylation sites.
  
## Requirments

1. Install Nextflow (version 23.04.3.5875)
2. Install R (version >= 4.1)
3. Install some R packages:
   + GenomicRanges
   + DESeq2
   + argparse
   + qs
   + import

## Directory structure

![directory structure](readme_figures/directory_structure.png)

This pipeline consists of four directories: 

+ bin: R scripts 
+ modules: Nextflow processes 
+ Reference: annotation files (gtf and fasta) from Cell Ranger 
+ data: It stores the required sample data in the ${sample_name}_illumina folder, and the site annotation (in qs format) in the apa_annotation folder. The sample data includes: 
	- outputs from various APA analysis methods, 
	- Illumina data (from Cell Ranger), 
	- ONT data, 
	- metadata.
		
## Usage

The usage of the benchmark pipeline is as follows:

```
nextflow run main.nf --sample_name MOBV1 --apa_method Sierra --core_num 4 --win_size 50 --cv_cutoff "none"
```

## Citation
@article {Zhao2024.11.29.626111,
	author = {Zhao, Qian and Rattray, Magnus},
	title = {Guidelines for alternative polyadenylation identification tools using single-cell and spatial transcriptomics data},
	elocation-id = {2024.11.29.626111},
	year = {2024},
	doi = {10.1101/2024.11.29.626111},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Background Many popular single-cell and spatial transcriptomics platforms exhibit 3{\textquoteright} bias, making it challenging to resolve all transcripts but potentially more feasible to resolve alternative polyadenylation (APA) events. Despite the development of several tools for identifying APA events in scRNA-seq data, a neutral benchmark is lacking, complicating the choice for biologists.Results We categorized existing APA analysis tools into three main classes, with the alignment-based class being the largest and we further divided this category into four sub-types. We compared the performance of methods from each algorithmic subtype in terms of site identification, quantification, and differential expression analysis across four single-cell and spatial transcriptomic datasets, using matched nanopore data as ground truth. No single method showed absolute superiority in all comparisons. Therefore, we selected representative methods (Sierra, scAPAtrap, and SCAPE) to deeply analyze the impact of different algorithmic choices on performance. SCAPE which is based on the distance estimation demonstrated less sensitivity to changes in read length and sequencing depth. It identified the most sites and achieved high recall but does not account for the impact of alternative splicing on site identification, leading to a loss in precision. Sierra that fits a coverage distribution is sensitive to changes in sequencing depth and identifies relatively fewer sites, but it considers the impact of junction reads on site identification and this results in relatively high precision. scAPAtrap combines peak calling and soft clipping, both of which are sensitive to sequencing depth. Moreover, soft clipping is particularly sensitive to read length, with increased read length leading to more false positive sites. Quantification consistency was affected by Cell Ranger versions and parameters, influencing downstream analysis but having less effect on differential expression between cell types.Conclusions Each method has unique strengths. SCAPE is recommended for low-coverage data, scAPAtrap for moderate read lengths including intergenic sites, and Sierra for high-depth data with alternative splicing considerations. Filtering low-confidence sites, choosing appropriate mapping tools, and optimizing window size can improve performance.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2024/12/04/2024.11.29.626111},
	eprint = {https://www.biorxiv.org/content/early/2024/12/04/2024.11.29.626111.full.pdf},
	journal = {bioRxiv}
}