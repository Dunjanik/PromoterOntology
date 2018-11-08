scan_promoters <- function(seq, pwm, min_score, save_intermediate_files=FALSE, label=NA){
  #' Scan promoters with PWM and return hits in a tidy format
  #'
  #' @param seq A DNAStringSet that will be scanned with PWM
  #' @param pwm PWMatrix or PWMatrixList object that sequences are going to be scanned against.
  #' @param min_score A minimum sequence identity needed to report sequence fragment as hit.
  #'            Can be given an character string in the format of "80%" or as a single absolute value between 0 and 1.
  #' @param save_intermediate_files logical: do you want to save a SiteSet object from sequence scan
  #' @param label if save_intermediate_files = TRUE, label will be used as a filename
  #'
  #' @return A data.frame that reports all hits from sequence seach in tidy format.
  #' @importFrom TFBSTools searchSeq
  #'
  #' @keywords
  #' @export
  message("Scanning Sequences")
  results <- searchSeq(pwm, subject = seq, seqname = label, min.score = min_score, strand = "+")

  if (save_intermediate_files) saveRDS(results, paste0("searchSeq_", label, ".RDS"))

  # lets calculate matscore for all thresholds :)
  all_hits <- lapply(results, function(x) as.data.frame(x))
  all_hits <- do.call(rbind, all_hits)
  all_hits <- all_hits[, c("seqnames", "start", "end", "absScore", "TF", "siteSeqs")]

  return(all_hits)
}
