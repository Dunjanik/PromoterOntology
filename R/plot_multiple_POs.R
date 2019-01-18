plot_multiple_POs <- function(list_POs, use_cage=TRUE, tunable_plot = TRUE, title = "Promoter structer for these gene groups:"){
  #'
  #' Plot mutliple PO results in one barchart
  #'
  #' This function allows to plot results of promoter ontology for multiple groups of genes in one plot, making comparison between different gene groups easy. Function returns ggplot2 object, so in the case additional customisation of plot is needed, extra parameters can be added additionaly with "+" sign.
  #'
  #' @param list_POs a list object where each sublist is a result of promoter_onto function
  #' @param use_cage logical specifying whether to plot
  #' @param tunable_plot logical, specifiying if all results will be presented in a single plot or the outplut will be a modular ggplot object that can later be faceted...
  #' @param title title of a plot used only in a non-tunable version
  #'
  #' @return plot object to which you can add extra ggplot parameters
  #' @import tidyverse
  #'
  #' @export

  n <- length(list_POs)

  if(tunable_plot == TRUE){
    ## keep only one test from over-under pair
    list_POs <- lapply(list_POs, function(x){
      x[[1]] <- x[[1]] %>%
        group_by(TF) %>%
        top_n(-1, pval);
      x
    })
  }

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

  ## reorder features on the graph based on significance
  order <- tf_df %>% group_by(feature) %>% summarise(cum = sum(abs(pval))) %>%
    arrange(cum) %>% dplyr::select(feature)
  order <- order$feature
  tf_df$feature <- factor(tf_df$feature, levels = order)
  tf_df$pval <- tf_df$pval * ifelse(tf_df$class == "under", 1, -1)
  ## coords for annotation text
  coords <- c(min(tf_df$pval), max(tf_df$pval))

  if(tunable_plot == TRUE){
    ## number of features to serve as coordinates for hline
    l <- length(unique(tf_df$feature))
    l <- 1:(l-1)
    l <- l + 0.5

    pl <- ggplot(tf_df, aes(x=feature, y=pval, class=group, fill=group))+
      geom_bar(width=0.75, stat="identity", position=position_dodge())+
      theme_minimal()+ coord_flip() +
      theme(axis.text.y = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))+
      geom_hline(yintercept = log10(0.05), colour="red", linetype=3)+
      geom_hline(yintercept = -log10(0.05), colour="red", linetype=3)+
      scale_fill_manual(values=wes_palette(n=(n+1), name="Darjeeling1", type="continuous")[2:(n+1)])+
      geom_vline(xintercept = l, colour="grey", linetype=3)+
      ylab("-log10(pval)")+
      labs(subtitle="BH multiple testing corrected")+
      annotate("text", y = coords[1] * 0.88, x = 1, label = "italic(Depleted)", parse = TRUE) +
      annotate("text", y = coords[2] * 0.88, x = 1, label = "italic(Enriched)", parse = TRUE)
    return(pl)
  }else{
    tf_df$pval <- abs(tf_df$pval)
    plot_under <- ggplot(filter(tf_df, class == "under"), aes(x=feature, y=pval, class=group, fill=group))+
      geom_bar(width=0.75, stat="identity", position=position_dodge())+
      theme_minimal() + coord_flip()+scale_y_reverse()+
      scale_fill_manual(values=wes_palette(n=(n+1), name="Darjeeling1", type="continuous")[2:(n+1)])+
      ylab("-log10(pval)")+ guides(fill=FALSE)+theme(axis.text.y = element_text(size=12))+
      annotate("text", y = abs(coords[1]) * 0.88 , x = 1, label = "italic(Depleted)", parse = TRUE) +
      labs(subtitle = "Depleted features", x = "Promoter features")+
      geom_hline(yintercept = -log10(0.05), colour="red", linetype=3)

    plot_over <- ggplot(filter(tf_df, class == "over"), aes(x=feature, y=pval, class=group, fill=group))+
      geom_bar(width=0.75, stat="identity", position=position_dodge())+
      theme_minimal() + coord_flip()+
      scale_fill_manual(values=wes_palette(n=(n+1), name="Darjeeling1", type="continuous")[2:(n+1)])+
      annotate("text", y = coords[2] * 0.88, x = 1, label = "italic(Enriched)", parse = TRUE)+
      ylab("-log10(pval)")+
      labs(subtitle = "Enriched features", x = "Promoter features")+
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank())+
      geom_hline(yintercept = -log10(0.05), colour="red", linetype=3)

    grid.arrange(plot_under, plot_over, nrow = 1,
                 top=title)
  }

}
