create_promoter_matrix <- function(all_hits, optimal_threshold, promoters){
  #'
  #'
  #' @param all_hits
  #' @param optimal_threshold
  #' @param promoters
  #' @param
  #' @param
  #' @param
  #'
  #'
  #'
  #' @return promoter_matrix object that for each promoter annotates
  #' @importFrom
  #' @importFrom
  #'
  #' @keywords
  #' @export

  ## winsorize all hits to quantiles
  qnt <- quantile(all_hits$absScore, probs=c(.15, .85), na.rm = T)
  H <- IQR(all_hits$absScore, na.rm = T) * 1.5 # intarquantile range
  all_hits$absScore[all_hits$absScore < (qnt[1] - H)] <- (qnt[1] - H)
  all_hits$absScore[all_hits$absScore > (qnt[2] + H)] <- (qnt[2] + H)


  # now that we have capped all values, will just normalise by max
  norm_score <- all_hits
  norm_score$normScore <- normalise(norm_score$absScore)
  norm_score$loc_prob <- prom_freq_norm(norm_score$start)[norm_score$start]
  norm_score$seqnames <- as.numeric(as.character(norm_score$seqnames))

  ### combine scores
  norm_score$comb <- norm_score$normScore * norm_score$loc_prob
  norm_score <- norm_score[order(norm_score$comb, decreasing = TRUE), ]

  ## threshold all hits
  norm_score <- norm_score[norm_score$comb > optimal_threshold, ]

  promoter_matrix <- mcols(promoters)[, c("external_gene_name", "interquantile_width")] #, "tpm.dominant_ctss", "tpm"
  promoter_matrix$TF <- 0
  promoter_matrix$TF[norm_score$seqnames] <- 1
  colnames(promoter_matrix)[ncol(promoter_matrix)] <- as.character(all_hits$TF[1])

  return(promoter_matrix)
}
