get_biomart_promoters <- function (mart, upstream=500, downstream=500, output="GR", protein_coding_only=TRUE) {
  
  #' Extract promoter regions from biomaRt object
  #'
  #' @param mart object of class Mart created using the useMart from BiomaRt package
  #' @param upstream upstream region from TSS to be included in promoter region
  #' @param downstream downstream region from TSS to be included in promoter region
  #' @param output character defining output type. 
  #'               use "DF" and "GR" or "data.frame" and "GenomicRanges"
  #' @param protein_coding_only logical: will returned promoters be promoter of 
  #'              only promoter coding genes 
  #' 
  #' This function extracts promoter regions from a specified mart object. 
  #' Each promoter will be annotated with gene_biotype, tpm expression of a 
  #' promoter and interquantile width of a promoter.
  #' Width of extracted promoters is defined by upstream and downstream parameters. 
  #' Generated promoter regions can be exported as GenomicRanges object or, 
  #' as a data frame. 
  #' 
  #' @return a promoter regions of specified width as a data frame or GenomicRanges object 
  #' @importFrom biomaRt getBM
  #' @importFrom GenomicRanges makeGRangesFromDataFrame promoters
  #' @importFrom GenomeInfoDb seqlevelsStyle
  #' 
  #' @keywords promoters
  #' @export
  
  
  biomart_genes <- getBM(attributes = c("ensembl_gene_id",
                                      "ensembl_transcript_id",
                                      "external_gene_name",
                                      "chromosome_name",
                                      "transcript_start",
                                      "transcript_end",
                                      "strand",
                                      "gene_biotype"),
                       mart = mart)
  
  if (protein_coding_only) biomart_genes <- biomart_genes[biomart_genes$gene_biotype == "protein_coding",]
  
  biomart_genes$strand <- ifelse(biomart_genes$strand == 1, "+", 
                                 ifelse(biomart_genes$strand == -1, "-", "*"))
  biomart_genes_gr <- makeGRangesFromDataFrame(biomart_genes,
                                             start.field = "transcript_start",
                                             end.field = "transcript_end",
                                             keep.extra.columns = TRUE)
  seqlevelsStyle(biomart_genes_gr) <- "UCSC"
  
  #make promoters using coordinates upstream, downstream
  biomart_promoters <- promoters(biomart_genes_gr,
                               upstream = upstream,
                               downstream = downstream)
  biomart_promoters <- unique(biomart_promoters)
  
  if (tolower(output) %in% c("gr", "genomicranges")){
    return(biomart_promoters)} else if (tolower(output) %in% c("df", "data.frame", "dataframe")){
      promoters_df <- as.data.frame(biomart_promoters)
      names(promoters_df)[1] <- "chr"
      return(promoters_df)} else 
        stop("Output has to be either GenomicRanges promoter object or data.frame!")
}

