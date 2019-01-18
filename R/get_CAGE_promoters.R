get_CAGE_promoters <- function (tc_gr, promoters, keep_all_promoters = FALSE, upstream=500, downstream=500,  protein_coding_only=TRUE) {
  #' Create promoters centered on dominant TSS CAGE tag clusters with annotations from Biomart promoters.
  #'
  #' @param tc_gr GenomicRanges object containing tag clusters extracted from CAGEset
  #' @param promoters GenomicRanges object containing promoters extracted from Biomart
  #' @param keep_all_promoters logical: in case there are multiple TC assigned to a promoter,
  #'                  do you want to keep all of them or just the one with highest expression
  #' @param upstream upstream region from TSS to be included in promoter region
  #' @param downstream downstream region from TSS to be included in promoter region
  #' @param protein_coding_only logical: will returned promoters be promoter of
  #'              only promoter coding genes
  #'
  #' This function creates promoter ranges from CAGE tag clusters.
  #'
  #' @return a promoter regions centered on dominant TSS extracted from CAGE tag clusters
  #' @importFrom GenomicRanges resize findOverlaps
  #' @import dplyr
  #'
  #' @keywords promoters
  #' @export

  tc_prom_overlap <- findOverlaps(tc_gr, promoters)
  tc_prom_overlap_df <- as.data.frame(tc_prom_overlap)
  tc_prom_overlap_df$tpm <- tc_gr$tpm[tc_prom_overlap_df$queryHits]
  tc_prom_overlap_df$ensembl_gene_id <-
    promoters$ensembl_gene_id[tc_prom_overlap_df$subjectHits]

  #retain only one transcript per hit
  tc_prom_overlap_df <-
    tc_prom_overlap_df[!duplicated(tc_prom_overlap_df[,c("queryHits", "ensembl_gene_id")]), ]


  # if multiple TCs hit the same gene, retain only the TC with most tags: in case keep_all_promoters=FALSE
  tc_prom_overlap_df <- tc_prom_overlap_df %>% arrange(ensembl_gene_id, desc(tpm))
  if (keep_all_promoters == FALSE){
    tc_prom_overlap_df <- tc_prom_overlap_df[!duplicated(tc_prom_overlap_df[,c("ensembl_gene_id")]), ]
  }

  #### assign Ensembl id to TCs
  mcols(tc_gr)$ensembl_gene_id <- "NA"
  mcols(tc_gr)[tc_prom_overlap_df$queryHits, ]$ensembl_gene_id <-
    mcols(promoters)[tc_prom_overlap_df$subjectHits, ]$ensembl_gene_id

  genes_df <- mcols(promoters)[,c("ensembl_gene_id", "external_gene_name", "gene_biotype")]
  genes_df <- genes_df[!duplicated(genes_df), ]
  genes_with_tcs <- merge(genes_df, as.data.frame(tc_gr)[tc_prom_overlap_df$queryHits, ],
                          by = c("ensembl_gene_id"),
                          all= keep_all_promoters)

  ## filter out non-protein coding genes if param = TRUE
  if (protein_coding_only) genes_with_tcs <- genes_with_tcs[genes_with_tcs$gene_biotype == "protein_coding", ]
  genes_with_tcs <- genes_with_tcs[, -14:-15] # drop width quantile information

  cage_promoters <- makeGRangesFromDataFrame(genes_with_tcs[, -5:-6],
                                             start.field = "dominant_ctss",
                                             end.field = "dominant_ctss",
                                             keep.extra.columns = TRUE)
  cage_promoters <- promoters(cage_promoters, upstream, downstream)
  cage_promoters <- unique(cage_promoters)
  return(cage_promoters)
}
