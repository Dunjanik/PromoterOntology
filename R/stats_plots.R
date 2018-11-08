stats_plots <- function(all_hits){

  stat <- all_hits %>% group_by(start) %>% summarize(mean=mean(absScore), median=median(absScore), sd=sd(absScore))
  max_prom <- all_hits %>% group_by(seqnames) %>% summarise(loc=start[which.max(absScore)])


  bin_size <- lapply(pwmMatches, function(x) c(round(ceiling(x / 10) * 1:9), x))
  bin_prob <- mapply(function(x, y) data.frame(sort(x$prob)[y], stringsAsFactors = FALSE), x=norm_score, y=bin_size)
  names(bin_prob) <- TF_names
  bin_score <- mapply(function(x, y) data.frame(sort(x$score)[y], stringsAsFactors = FALSE), x=norm_score, y=bin_size)
  names(bin_score) <- TF_names


  ## plot binned heatmap
  for(i in 1:length(norm_score)){
    norm_score[[i]]$bin.prob <- findInterval(norm_score[[i]]$prob, bin_prob[[i]])
    norm_score[[i]]$bin.score <- findInterval(norm_score[[i]]$score, bin_score[[i]])
    # now i am getting cases with max value as sixth cathegory so i will slightly increase the last baundry

    norm_score[[i]]$bin.score[norm_score[[i]]$bin.score == 10] <- 9
    norm_score[[i]]$bin.prob[norm_score[[i]]$bin.prob == 10] <- 9

    norm_score[[i]]$class <- paste0(norm_score[[i]]$bin.prob, "_", norm_score[[i]]$bin.score)

  }


  classes <- lapply(norm_score, function(x) table(x$class))

  my_palette <- colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(n = 299)

  ref.clas5 <- c(paste0(rep(0:4, each=5), "_", rep(0:4, times=5)))
  ref.clas10 <- c(paste0(rep(0:9, each=10), "_", rep(0:9, times=10)))

  pdf(paste0("PWMheatmap5_new", stageID, "_11.pdf"))
  for(i in 1:length(classes)){
    mat <- rep(0, times=100)
    names(mat) <- ref.clas10
    mat[names(mat) %in% names(classes[[i]])] <- classes[[i]]
    mat <- matrix(mat, nrow =10, byrow = T)
    heatmap.2(mat, Colv = "NA", Rowv = FALSE,
              dendrogram="none", col=my_palette, trace="none",
              main = TF_names[i], xlab = "score bin", ylab="position bin")
  }

  dev.off()


  pdf(paste0("scoreProbCombHistograms_", stageID, "_11.pdf"))
  for(i in 1:length(norm_score)){
    hist(norm_score[[i]]$comb, breaks=500, main=TF_names[i])
    hist(norm_score[[i]]$comb, breaks=500, xlim=c(0, 0.2), main=TF_names[i])
    plot(ecdf(norm_score[[i]]$comb))
  }
  dev.off()

}
