plot_multiple_POs <- function(list_POs, use_cage=TRUE){
  #'
  #' Plot mutliple PO results in one barchart
  #'
  #' This function allows to plot results of promoter ontology for multiple groups of genes in one plot, making comparison between different gene groups easy. Function returns ggplot2 object, so in the case additional customisation of plot is needed, extra parameters can be added additionaly with "+" sign.
  #'
  #' @param list_POs a list object where each sublist is a result of promoter_onto function
  #' @param use_cage logical specifying whether to plot
  #'
  #'
  #' @return plot object to which you can add extra ggplot parameters
  #' @import tidyverse
  #'
  #'

  n <- length(list_POs)

  ## keep only one test from over-under pair
  list_POs <- lapply(list_POs, function(x){
    x[[1]] <- x[[1]] %>%
      group_by(TF) %>%
      top_n(-1, pval);
    x
  })

  tf_df <- lapply(list_POs, function(x){
    data.frame(feature = as.character(x[[1]]$TF),
               pval = p.adjust(x[[1]]$pval, "BH"),
               class = as.character(x[[1]]$class))
  })
  tf_df <- do.call(rbind, tf_df)
  tf_df$group <- rep(names(list_POs), times=sapply(list_POs, function(x) nrow(x[[1]])))

  if(use_cage == TRUE){
    cage_df <- lapply(list_POs, function(x){
      data.frame(feature = x[[2]]$feature,
                 pval = p.adjust(x[[2]]$pval, "BH"),
                 class = x[[2]]$class)
    })
    cage_df <- do.call(rbind, cage_df)
    cage_df$group <- rep(names(list_POs), times=sapply(list_POs, function(x) nrow(x[[2]])))
    tf_df <- rbind(tf_df, cage_df)

  }

  tf_df$pval <- log10(tf_df$pval)
  ## lets see if this is going to help with plotting
  tf_df$pval <- tf_df$pval * ifelse(tf_df$class == "under", 1, -1)

  plot <- ggplot(tf_df, aes(x=feature, y=pval, class=group, fill=group))+
    geom_bar(stat="identity", position=position_dodge())+
    theme_minimal()+
    geom_hline(yintercept = log10(0.05), colour="red", linetype=3)+
    geom_hline(yintercept = -log10(0.05), colour="red", linetype=3)+
    scale_fill_manual(values=wes_palette(n=(n+1), name="Darjeeling1", type="continuous")[2:(n+1)])+
    ylab("-log10(pval)")+
    labs(subtitle="BH multiple testing corrected")

  return(plot)

}
