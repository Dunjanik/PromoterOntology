find_threshold <- function(all_hits, constraints="default"){
  #'
  #' A function that creates normalises score and location position
  #' from sequence scann hits and suggests an optimal threshold for filtering
  #' false positive hits.
  #' To filter false positives we have implemented location constraints of core promoter TFs
  #' To use default constraints of these PWMs, use "default" parameter for constraints parameter.
  #' User can also submit custom constraints too.
  #'
  #' @param all_hits A dataframe containint all hits reported by sequence scan
  #' @param constraints location constraints where PWM can be found. User can use default
  #' defined location constraints, or define custom constraints in form of a numeric vector with start and end coordinate.
  #'
  #'
  #' @return suggested threshold value for filtering out false positive hits
  #' @importFrom inflection findiplist
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


  ### combine scores
  norm_score$comb <- norm_score$normScore * norm_score$loc_prob
  norm_score <- norm_score[order(norm_score$comb, decreasing = TRUE), ]

  ## lets have only one-best hit per promoter
  norm_score <- norm_score[!duplicated(norm_score$seqnames), ]

  if (constraints == "default"){
    #data("tf_constraints")
    n <- which(tf_constraints$TF == unique(norm_score$TF))
    if (length(n) == 0){
      message("No position constraints choosen.")
      cons <- unique(norm_score$start)
    } else {
      cons <- as.numeric(tf_constraints[n, c(2:3)])
      cons <- cons + 500
      cons <- seq(from=cons[1], to=cons[2])
    }
  } else {
    if (!is.numeric(constraints) & length(constraints) != 2){
      stop("Constraints should be a vector of two coordinates!")
    } else {
      cons <- constraints + 500
      cons <- seq(from=cons[1], to=cons[2])
    }
  }

  norm_score <- norm_score[norm_score$start %in% cons, ]

  ## after filtering out hits outside defined regions define optimal threshold for combined score

  diff_x <- diff(norm_score$comb)

  nr <- length(diff_x)
  w <- 5
  i <- 1:nr
  iw <- embed(i,w)[, w:1]
  dfw <- apply(iw, 2, function(x) diff_x[x])

  min <- rowSums(dfw)[which.min(rowSums(dfw))]
  whi <- which(rowSums(dfw) == min)

  infl <- c(FALSE, diff(abs(diff(norm_score$comb))>0)!=0)

  A <- findiplist(norm_score$comb, 1:nrow(norm_score), 0, doparallel=TRUE)
  thr <- c(norm_score$comb[max(whi) + 1], A[1, 3])

  optimal_threshold <- thr[which.min(thr)]

  return(optimal_threshold)
}



normalise <- function(x){
  min <- min(x)
  max <- max(x)
  norm <- (x - min)/(max - min)
  return(norm)
}



prom_freq_norm <- function(x){
  ta <- table(x)
  n_ta <- as.numeric(names(ta))
  missin <- (1:max(x))[!(1:max(x) %in% as.numeric(names(ta)))]
  y <- rep(0, length(missin))
  ta <- c(ta, y)
  names(ta) <- c(n_ta, missin)
  ta <- ta[order(as.numeric(names(ta)))]

  ta <- normalise(ta)
  return(ta)
}
