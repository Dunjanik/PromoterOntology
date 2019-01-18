filter_hits <- function(all_hits, threshold){
  #'
  #' A function that filters out false positive sequence scan hits once threshold is provided.
  #'
  #' @param all_hits A dataframe containint all hits reported by sequence scan
  #' @param threshold A threshold of combined score of location and sequence similarity that is used to filter out false positive scan hits.
  #'
  #'
  #' @return a dataframe with hits that are considered to be true positive sequence scan hits
  #'
  #' @keywords
  #' @export

  ## winsorize all hits to quantiles
  qnt <- quantile(all_hits$absScore, probs=c(.15, .85), na.rm = T)
  H <- IQR(all_hits$absScore, na.rm = T) * 1.5 # intarquantile range
  all_hits$absScore[all_hits$absScore < (qnt[1] - H)] <- (qnt[1] - H)
  all_hits$absScore[all_hits$absScore > (qnt[2] + H)] <- (qnt[2] + H)


  # now that we have capped all values, will just normalise by max
  all_hits$normScore <- normalise(all_hits$absScore)
  all_hits$loc_prob <- prom_freq_norm(all_hits$start)[all_hits$start]


  ### combine scores
  all_hits$comb <- all_hits$normScore * all_hits$loc_prob
  all_hits <- all_hits[order(all_hits$comb, decreasing = TRUE), ]

  ## lets have only one-best hit per promoter
  all_hits <- all_hits[!duplicated(all_hits$seqnames), ]

  all_hits <- all_hits[all_hits$comb >= threshold, ]

  return(all_hits)
}


