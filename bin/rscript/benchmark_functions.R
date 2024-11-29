# extract expression matrix
apa_matrix_extract <- function(expr_file, annot, method, barcodes){
  
  if(method == "Sierra"){
    
    ## file path
    sierra_mat_file <- file.path(expr_file, "Sierra_counts", "matrix.mtx.gz")
    sierra_sit <- file.path(expr_file, "Sierra_counts", "sitenames.tsv.gz")
    sierra_bcd <- file.path(expr_file, "Sierra_counts", "barcodes.tsv.gz")
    
    ## Read in expression matrix
    counts <- Matrix::readMM(sierra_mat_file)
    gene_ids <- readr::read_tsv(sierra_sit, col_names = FALSE, show_col_types = FALSE)$X1
    cell_ids <- readr::read_tsv(sierra_bcd, col_names = FALSE, show_col_types = FALSE)$X1
    
    ## Make the column names as the cell IDs and the row names as the gene IDs
    rownames(counts) <- gene_ids
    colnames(counts) <- gsub("-1", "", cell_ids)
    counts <- as.data.frame(as.matrix(counts))
    
    counts <- counts %>% 
      tibble::rownames_to_column("peakID") %>% 
      tidyr::separate(col = peakID, into = c("gene_name", "chr", "region", "strand"), sep = ":") %>% 
      tidyr::separate(col = region, into = c("start", "end"), sep = "-") %>% 
      dplyr::mutate(strand = dplyr::if_else(strand == 1, "+", "-"),
                    coord = dplyr::if_else(strand == "+", end, start)) %>% 
      tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":")
  }
  
  if(method == "SCAPE"){
    
    ## file path
    scape_mat_file <- file.path(expr_file, "pasite.csv.gz")
    
    ## rename and reformat
    counts <- vroom::vroom(scape_mat_file, show_col_types = F) %>% 
      tibble::column_to_rownames("...1") %>% 
      dplyr::rename_with(~gsub("-1", "", .x)) %>% 
      tibble::rownames_to_column("peakID") %>% 
      tidyr::separate(col = "peakID", into = c("chr", "coord", "score", "strand"), sep = ":") %>% 
      tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":") 
    
  }
  
  if(method == "scAPAtrap"){
    ## file path
    scapatrap_expr_file <- file.path(expr_file, "expma.qs")
    
    ## rename
    counts <- qs::qread(scapatrap_expr_file) %>% tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":")
    
  }
  
  if(method == "polyApipe"){
    
    ## file path
    polyapipe_mat_file <- list.files(expr_file, pattern = "counts.tab.gz", full.names = T)
    
    ## rename and reformat
    counts <- vroom::vroom(polyapipe_mat_file, show_col_types = FALSE) %>% 
      dplyr::filter(cell %in% barcodes) %>% 
      tidyr::separate(col = gene, into = c("chr", "coord", "strand"), sep = "_") %>%
      dplyr::mutate(strand = recode(strand, "r" = "-", "f" = "+")) %>% 
      tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":") %>% 
      tidyr::pivot_wider(names_from = cell, values_from = count, values_fill = 0) 
    
  }
  
  if(method == "SCAPTURE"){
    
    ## file path
    scapture_mat_file <- list.files(expr_file, pattern = "PASquant.KeepCell.UMIs.tsv.gz", full.names = T)
    scapture_bed_file <- list.files(expr_file, pattern = "PASquant.KeepPAS.bed", full.names = T)
    ## rename and reformat
    
    scapture_bed <- vroom::vroom(scapture_bed_file, col_names = F) %>% 
      dplyr::mutate(pa_site = dplyr::if_else(X6 == "+", X3, X2)) %>% 
      tidyr::unite(col = "peakID", c(X1, pa_site, X6), sep = ":") %>% 
      dplyr::rename(peak_name = X4) %>% 
      dplyr::select(peakID, peak_name)
    
    counts <- vroom::vroom(scapture_mat_file, show_col_types = FALSE) %>% 
      tibble::column_to_rownames("gene") %>% 
      dplyr::rename_with(~gsub("-1", "", .x)) %>% 
      tibble::rownames_to_column("peak_name") %>%
      dplyr::left_join(scapture_bed, by = "peak_name")

  }

  if(method == "scAPA"){
    
    ## file path
    peak_file <- list.files(expr_file, pattern = "Peaks.RDS", full.names = T, recursive = T)
    
    ## rename and reformat
    scAPA_bed <- peak_file %>% 
      readRDS(.) %>% 
      .@row.Data %>% 
      tidyr::separate_rows(Chr, Start, End, Strand, sep = ";") %>%
      dplyr::group_by(GeneID) %>%
      dplyr::summarise(
        Chr = unique(Chr),
        Start = dplyr::first(Start),
        End = dplyr::last(End),
        Length = dplyr::first(Length),
        Strand = unique(Strand)
      ) %>% 
      dplyr::mutate(
        Chr = stringr::str_c("chr", Chr),
        coord = dplyr::if_else(Strand == "+", End, Start)) %>% 
      tidyr::unite(col = "peakID", c("Chr", "coord", "Strand"), sep = ":") %>% 
      dplyr::select(GeneID, peakID)
    
    counts <- peak_file %>% 
      readRDS(.) %>% 
      .@cells.counts %>% 
      tibble::remove_rownames() %>% 
      tibble::column_to_rownames("Peak_ID") %>% 
      dplyr::rename_with(~gsub(".*_", "", .x)) %>% 
      dplyr::rename_with(~gsub("-1", "", .x)) %>% 
      tibble::rownames_to_column("GeneID") %>%
      dplyr::left_join(scAPA_bed, by = "GeneID")
    
  }

  if(method == "MAAPER"){
    peak_file <- list.files(expr_file, pattern = "pas.txt", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file) %>% dplyr::pull(pas)
    ## create a pseudo counts matrix
    counts <- matrix(
      sample(1:10, length(barcodes) * length(pas_name), replace = TRUE), 
      nrow = length(pas_name), 
      dimnames = list(pas_name, barcodes)
    ) %>% 
      dplyr::as_tibble(.) %>% 
      dplyr::mutate(peakID = pas_name)
    
  }

  if(method == "DarPars2"){
    
    file_list <- list.files(expr_file, pattern = "^Dapars2_result_temp.*\\.txt$", recursive = TRUE, full.names = TRUE)
    
    pas_name <- lapply(file_list, function(darpars_output){
      
      pa_site <- vroom::vroom(darpars_output) %>% 
        tidyr::separate(col = Gene, into = c("esembl_id", "symbol", "chr", "strand"), sep = "\\|") %>% 
        tidyr::separate(col = Loci, into = c("chr", "region"), sep = ":") %>% 
        tidyr::separate(col = region, into = c("start", "end"), sep = "-") %>% 
        dplyr::mutate(
          distal_sites = as.integer(dplyr::if_else(strand == "+", end, start)),
          Predicted_Proximal_APA = as.integer(Predicted_Proximal_APA)) %>% 
        tidyr::pivot_longer(col = c("Predicted_Proximal_APA", "distal_sites"), names_to = "type", values_to = "coord") %>% 
        tidyr::unite(col = "pas_name", c("chr", "coord", "strand"), sep = ":") %>% 
        dplyr::distinct(pas_name)
      
      return(pa_site)
    }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::pull(pas_name)
    
    ## create a pseudo counts matrix
    counts <- matrix(
      sample(1:10, length(barcodes) * length(pas_name), replace = TRUE), 
      nrow = length(pas_name), 
      dimnames = list(pas_name, barcodes)
    ) %>% 
      dplyr::as_tibble(.) %>% 
      dplyr::mutate(peakID = pas_name)
  }

  ## add annotation
  counts <- counts %>% 
    dplyr::left_join(annot, by = c("peakID" = "pa_site")) %>% 
    dplyr::filter(!is.na(seqnames))
  
  
  return(counts)
}

# annotate pa sites by SCAPE
apa_annotation <- function(apa_dir, method=c("SCAPE", "Sierra", "scAPAtrap"), gtf_file, barcode_file, species="Mm10", core_num=8){
  require(magrittr)
  require(SCAPE)
  require(GenomicRanges)
  Sys.setenv(VROOM_CONNECTION_SIZE = 500072)
  
  if(method == "SCAPE"){
    peak_file <- list.files(apa_dir, pattern = "pasite.csv.gz", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file, col_select = 1, col_names = c("pas_site"))
    pas_name <- pas_name %>% 
      dplyr::filter(!is.na(pas_site)) %>% 
      tidyr::separate(col = "pas_site", into = c("chr", "coord", "score", "strand"), sep = ":") %>% 
      tidyr::unite(col = "pas_site", c("chr", "coord", "strand"), sep = ":") %>% 
      dplyr::pull(pas_site)
  }
  
  if(method == "Sierra"){
    peak_file <- list.files(apa_dir, pattern = "sitenames.tsv.gz", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file, col_names = c("gene_name", "chr", "region", "strand"))
    pas_name <- pas_name %>% 
      dplyr::mutate(strand = dplyr::if_else(strand=="1", "+", "-")) %>% 
      tidyr::separate(col="region", into=c("start", "end"), sep = "-") %>% 
      dplyr::mutate(sitename = dplyr::if_else(strand=="+", end, start)) %>% 
      tidyr::unite(col = "peakID", c("chr", "sitename", "strand"), sep = ":") %>% 
      dplyr::pull(peakID)
  }
  
  if(method == "scAPAtrap"){
    peak_file <- list.files(apa_dir, pattern = "expma.qs", full.names = T, recursive = T)
    expma <- qs::qread(peak_file, nthreads=8)
    pas_name <- expma %>% 
      dplyr::select(1:6) %>% 
      tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":") %>% 
      dplyr::pull(peakID)
  }
  
  if(method == "MAAPER"){
    peak_file <- list.files(apa_dir, pattern = "pas.txt", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file) %>% dplyr::pull(pas)
    
  }
  
  if(method == "SCAPTURE"){
    peak_file <- list.files(apa_dir, pattern = "PASquant.KeepPAS.bed", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file, col_names = F)
    pas_name <- pas_name %>% 
      dplyr::mutate(pa_site = dplyr::if_else(X6 == "+", X3, X2)) %>% 
      dplyr::mutate(pas_name = stringr::str_c(X1, pa_site, X6, sep = ":")) %>% 
      dplyr::pull(pas_name)
  }
  
  if(method == "polyApipe"){
    peak_file <- list.files(apa_dir, pattern = "counts.tab.gz", full.names = T, recursive = T)
    barcodes <- readr::read_tsv(barcode_file, col_names = FALSE, show_col_types = FALSE)$X1 %>% gsub("-1", "", .)
    
    pas_name <- vroom::vroom(peak_file) %>% 
      dplyr::filter(cell %in% barcodes) %>% 
      dplyr::group_by(gene) %>% 
      dplyr::summarise(across(count, sum)) %>% 
      dplyr::filter(count > 0) %>% 
      tidyr::separate(col = gene, into = c("chr", "coord", "strand"), sep = "_") %>% 
      dplyr::mutate(strand = recode(strand, "r" = "-", "f" = "+")) %>% 
      tidyr::unite(col = "pas_name", c("chr", "coord", "strand"), sep = ":") %>% 
      dplyr::pull(pas_name)
  }
  
  if(method == "scAPA"){
    require(scAPA)
    peak_file <- list.files(apa_dir, pattern = "Peaks.RDS", full.names = T, recursive = T)
    
    pas_name <- peak_file %>% 
      readRDS(.) %>% 
      .@row.Data %>%
      tidyr::separate_rows(Chr, Start, End, Strand, sep = ";") %>%
      dplyr::group_by(GeneID) %>%
      dplyr::summarise(
        Chr = unique(Chr),
        Start = dplyr::first(Start),
        End = dplyr::last(End),
        Length = dplyr::first(Length),
        Strand = unique(Strand)
      ) %>% 
      dplyr::mutate(
        Chr = stringr::str_c("chr", Chr),
        coord = dplyr::if_else(Strand == "+", End, Start)) %>% 
      tidyr::unite(col = "pas_name", c("Chr", "coord", "Strand"), sep = ":") %>% 
      dplyr::pull(pas_name)
  }

  if(method == "DarPars2"){
    
    file_list <- list.files(apa_dir, pattern = "^Dapars2_result_temp.*\\.txt$", recursive = TRUE, full.names = TRUE)
    
    pas_name <- lapply(file_list, function(darpars_output){
      pa_site <- vroom::vroom(darpars_output) %>% 
        tidyr::separate(col = Gene, into = c("esembl_id", "symbol", "chr", "strand"), sep = "\\|") %>% 
        tidyr::separate(col = Loci, into = c("chr", "region"), sep = ":") %>% 
        tidyr::separate(col = region, into = c("start", "end"), sep = "-") %>% 
        dplyr::mutate(
          distal_sites = as.integer(dplyr::if_else(strand == "+", end, start)),
          Predicted_Proximal_APA = as.integer(Predicted_Proximal_APA)) %>% 
        tidyr::pivot_longer(col = c("Predicted_Proximal_APA", "distal_sites"), names_to = "type", values_to = "coord") %>% 
        tidyr::unite(col = "pas_name", c("chr", "coord", "strand"), sep = ":") %>% 
        dplyr::distinct(pas_name)

      return(pa_site)
    }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::pull(pas_name)
    
  }

  annot_info <- AnnotationSite(pas_name,
                               gtf_file,
                               species, #'Mm10'
                               cores = core_num) 
  
  return(annot_info)
}

# generate cell list
generate_cell_list <- function(cell_metadata, type = "all", core_num = 1L){
  cell_type_annot <- unique(as.character(cell_metadata$seurat_clusters))
  
  if(type == "sc"){
    cids <- cell_metadata$cell_ids
    names(cids) <- cell_metadata$cell_ids
  }
  
  if(type == "ct"){
    cids <- parallel::mclapply(
      cell_type_annot,
      FUN = function(x) {
        cell_ids <- cell_metadata %>%
          dplyr::filter(seurat_clusters == x) %>%
          dplyr::pull("cell_ids")
        return(cell_ids)},
      mc.cores = core_num)
    
    names(cids) <- cell_type_annot
  }
  
  if(type == "all"){
    cids <- list()
    cids[["all"]] <- cell_metadata$cell_ids
  }
  
  return(cids)
}

# convert peak ID to bed format dataframe
peak2bed <- function(peak_list){
  
  require(dplyr)
  
  peak_bed <- peak_list %>% 
    tidyr::separate(col = peakID, into = c("seqnames", "sites", "strand"), sep = ":") %>% 
    tidyr::unite(c("transcriptId", "geneId"), col = "name" , sep = "|") %>% 
    dplyr::mutate(start = dplyr::if_else(strand == "+", as.integer(sites), as.integer(sites) - 1),
                  end = dplyr::if_else(strand == "+", as.integer(sites) + 1, as.integer(sites))) %>% 
    dplyr::rename(score = total) %>% 
    dplyr::select(seqnames, start, end, name, score, strand)
  
  return(peak_bed)
}

# convert bed format dataframe to gr object
bed2gr <- function(peak_bed){
  
  require(GenomicRanges)
  
  gr_obj <- GRanges(
    seqnames = peak_bed$seqnames, 
    ranges = IRanges(start = peak_bed$start, end = peak_bed$end),
    names = peak_bed$name,
    strand = peak_bed$strand,
    score = peak_bed$score)
  
  return(gr_obj)
}

# find overlap sites between nanopore and apa data
find_overlap_sites_per_gene <- function(apa_data, nanopore_data, gene_name, win_size=0){
  
  require(dplyr)
  require(GenomicRanges)
  
  # convert apa dataframe into gr object
  apa_gr <- apa_data %>% 
    dplyr::filter(geneId == gene_name) %>% 
    peak2bed() %>% 
    bed2gr()
  
  # merge transcripts with the same end for nanopore data
  nanopore_gr <- nanopore_data %>% 
    dplyr::filter(geneId == gene_name) %>% #Aagab Tmem203 Rpl9
    peak2bed() %>% 
    dplyr::group_by(seqnames, start, end, strand) %>% 
    dplyr::summarise(name=stringr::str_c(name,collapse = ";"), score=sum(score)) %>% 
    bed2gr()
  
  # find overlap sites
  overlap_sites <- GenomicRanges::findOverlaps(
    query = apa_gr, 
    subject = nanopore_gr, 
    maxgap = win_size, 
    ignore.strand = FALSE)
  
  # true positive; predicted sites (apa length= tp + fp); ground truth (nanopore length=tp+fn)
  return(c(length(unique(overlap_sites@from)), overlap_sites@nLnode, overlap_sites@nRnode)) #return(overlap_sites)
}

# find overlap sites between nanopore and apa data and merge table
matched_sites_expression <- function(apa_data, nanopore_data, gene_name, win_size=0){
  # merge transcripts with the same end for nanopore data
  nanopore_gr <- nanopore_data %>% 
    dplyr::filter(geneId == gene_name) %>% #Aagab Tmem203 Rpl9
    peak2bed() %>% 
    dplyr::group_by(seqnames, start, end, strand) %>% 
    dplyr::summarise(name=stringr::str_c(name,collapse = ";"), score=sum(score)) %>% 
    bed2gr()
  
  # convert apa dataframe into gr object
  apa_gr <- apa_data %>% 
    dplyr::filter(geneId == gene_name) %>% 
    peak2bed() %>% 
    bed2gr()
  
  # find overlap sites
  overlap_sites <- GenomicRanges::findOverlaps(
    query = apa_gr, 
    subject = nanopore_gr, 
    maxgap = win_size, 
    ignore.strand = FALSE)
  
  # combine apa and nanopore dataframe
  if(length(overlap_sites@from) > 0){
    
    apa_df <- apa_gr[overlap_sites@from, ] %>% 
      as.data.frame() %>% 
      tidyr::unite(col="apa_sites", c("seqnames", "start", "end", "strand", "names"), sep = ":") %>% 
      dplyr::select(apa_sites, score) %>% 
      dplyr::rename(apa_score = score) %>% 
      dplyr::mutate(apa_rank = dplyr::dense_rank(apa_score))
    
    nanopore_df <- nanopore_gr[overlap_sites@to, ] %>% 
      as.data.frame() %>% 
      tidyr::unite(col="nanopore_sites", c("seqnames", "start", "end", "strand", "names"), sep = ":") %>% 
      dplyr::select(nanopore_sites, score) %>% 
      dplyr::rename(nanopore_score = score) %>% 
      dplyr::mutate(nanopore_rank = dplyr::dense_rank(nanopore_score))
    
    merge_df <- cbind(apa_df, nanopore_df)
    
    return(merge_df)
  } else {
    return(NULL)
  }
}

create_pseudobulk <- function(filtered_matrix, cell_type, cell_id, num_batches = 3, random_seed = 123){
  
  set.seed(random_seed)
  
  # ramdom generate batch
  batch_assignments <- sample(rep(1:num_batches, length.out = length(cell_id)))
  
  pseudobulk_counts <- lapply(1:num_batches, function(batch){
    rowSums(filtered_matrix[, batch_assignments == batch, drop = FALSE])
  }) %>% 
    setNames(paste0(cell_type, "_batch", 1:num_batches)) %>% 
    dplyr::bind_cols()
  
  return(pseudobulk_counts)
}

two_cluster_deseq2 <- function(matrix, cell_types, cell_metadata, num_batches = 3, random_seed = 123, count_filter = 0){
  
  require(dplyr)
  require(DESeq2)
  # create pseudobulk counts
  pseudobulk_counts <- lapply(cell_types, function(cell_type){
    
    cell_ids <- cell_metadata %>% dplyr::filter(seurat_clusters == cell_type) %>% dplyr::pull(cell_ids)
    
    pseudobulk_per_type <- matrix %>% 
      dplyr::select(dplyr::all_of(cell_ids)) %>% 
      create_pseudobulk(., cell_type, cell_ids, num_batches = num_batches, random_seed = random_seed)
    
    return(pseudobulk_per_type)
    
  }) %>% 
    dplyr::bind_cols(.name_repair = "unique") %>% 
    as.matrix()
  
  mode(pseudobulk_counts) <- "numeric"
  rownames(pseudobulk_counts) <- rownames(matrix)
  pseudobulk_counts <- pseudobulk_counts[rowSums(pseudobulk_counts) > count_filter,]
  
  # design matrix
  condition <- factor(rep(c(cell_type1, cell_type2), each = num_batches))
  
  # create DESeq object
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = pseudobulk_counts, 
    colData = data.frame(condition = condition), 
    design = ~ condition)
  
  dds <- DESeq2::DESeq(dds)
  
  res <- DESeq2::results(dds, contrast = c("condition", cell_type2, cell_type1)) %>% # type2 vs type1
    as.data.frame() %>% 
    tibble::rownames_to_column("peakID") 
  
  rm(dds); gc()
  # print(head(res[order(res$padj), ], 10))
  # summary(res)
  # write.csv(as.data.frame(res), file = "pseudobulk_differential_expression_genes.csv", row.names = TRUE)
  return(res)
}

