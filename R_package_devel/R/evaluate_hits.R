evaluate_hits <- function(all_hits ){


  ## winsorize all hits to quantiles
  qnt <- quantile(all_hits$absScore, probs=c(.15, .85), na.rm = T)
  caps <-  quantile(all_hits$absScore, probs=c(.05, .95), na.rm = T)
  H <- IQR(all_hits$absScore, na.rm = T) * 1.5 # intarquantile range
  all_hits$absScore[all_hits$absScore < (qnt[1] - H)] <- caps[1]
  all_hits$absScore[all_hits$absScore > (qnt[2] + H)] <- caps[2]


  # now that we have capped all values, will just normalise by max
  norm_score <- all_hits
  norm_score$normScore <- normalise(norm_score$absScore)
  norm_score$loc_prob <- prom_freq_norm(norm_score$start)[norm_score$start]


  ### combine scores
  norm_score$comb <- norm_score$normScore * norm_score$loc_prob
  norm_score <- norm_score[order(norm_score$comb, decreasing = TRUE), ]

  x <- ecdf(norm_score$comb)


  A <- findiplist(x(norm_score$comb), 1:nrow(norm_score), 0,doparallel=TRUE)


}

all_hits <- pwm_prom[[10]]
findiplist(norm_score$comb, 1:nrow(norm_score), 0, doparallel=FALSE)
plot(ecdf(norm_score$comb))

x <- norm_score$comb
y <- 1:nrow(norm_score)
lo <- loess(x~y, span = 0.5)
plot(y, x)
xl <- seq(min(y),max(y), (max(y) - min(y))/1000)
lines(xl, predict(lo,xl), col='red', lwd=2)

fit <- lm(y ~ x)
summary(fit)


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
# Evaluate hits:
#     plot all Hits
#     plot max hit per promoter
#     normalise score and position
#     plot binned hits heatmap
#
