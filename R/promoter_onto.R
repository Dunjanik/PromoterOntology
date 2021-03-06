promoter_onto <- function(promoter_matrix, sample, proms, features = c(1:11, 16:17), save_plot=TRUE, label, all=FALSE){
  #' Overenrichment analysis of promoter features for a group of genes
  #'
  #'
  #' Promoter_onto performs oveerenrichment analysis of promoter features for genes in sample compared to the background. Promoter_matrix object contains annotations of promoter features for a background set. Overenrichment analysis also looks at CAGE features like tpm expression and promoter interquantile width, for that reason promoter input is needed too. By default this function plots overenrichment plots as a pdf to the file in current directory.
  #'
  #' @param promoter_matrix promoter matrix that for each promoter (row), describes its core promoter features (columns)
  #' @param sample gene names for which we will calculate overrepresentation of features
  #' @param proms genomic ranges object containing promoters of all active genes centered on dominant TSS
  #' @param features which features to test, provided as vector of column indices of promoter_matrix matrix
  #' @param save_plot should results be saved, or just presented in console, if yes, file is going to be saved as PDF to working directory
  #' @param label filename for the plot
  #' @param all weather to report stats of all features (TRUE), or only significan ones (FALSE)
  #'
  #'
  #'
  #' @return a list with significance features of overrepresentation test
  #' @import tidyverse
  #' @import GenomicRanges
  #' @import ggthemes
  #' @import reshape2
  #' @import ggplot2
  #' @import wesanderson
  #'
  #' @keywords overrepresentation
  #' @export
  #'
  # function that plots all over/underrepresentation plots into a pdf
  # and returns a df with all features that are significantly overrepresented

  ### overrep of TF features
  overTF <- overrep_features(promoter_matrix, sample, features = features)
  overTF$TF <- rownames(overTF)
  overTF <- arrange(overTF, desc(pval))
  overTF$TF <- factor(overTF$TF, levels=unique(overTF$TF))
  overTF_plot <- ggplot(overTF, aes(x=geneRatio, y=TF, size=count, color=pval)) +
    geom_point(alpha=0.99) + ## illustrator catch
    scale_colour_gradientn(colours = colorRampPalette(c("brown", "red", "orange", "green", "blue"))(5),
                           values=c(0, 0.01,0.05,0.2,0.4, 1)) +
    ylab("TFs") + ggtitle("Feature overrepresentation in the sample") +
    theme_fivethirtyeight() + theme(plot.background = element_rect(fill = 'white'),
                                    panel.background = element_rect(fill = 'white'),
                                    legend.background = element_rect(fill = 'white'),
                                    axis.title = element_text(size = rel(1.4)),
                                    axis.text = element_text(size = rel(1.4)))

  ### underrep of TF features
  underTF <- overrep_features(promoter_matrix, sample, features = features, enrich = "under")
  underTF$TF <- rownames(underTF)
  underTF <- arrange(underTF, desc(pval))
  underTF$TF <- factor(underTF$TF, levels=unique(underTF$TF))
  underTF_plot <- ggplot(underTF, aes(x=geneRatio, y=TF, size=count, color=pval)) +
    geom_point(alpha=0.99) + ## illustrator catch
    scale_colour_gradientn(colours = colorRampPalette(c("brown", "red", "orange", "green", "blue"))(5),
                           values=c(0, 0.01,0.05,0.2,0.4, 1)) +
    ylab("TFs") + ggtitle("Feature underrepresentation in the sample") +
    theme_fivethirtyeight() +  theme(plot.background = element_rect(fill = 'white'),
                                     panel.background = element_rect(fill = 'white'),
                                     legend.background = element_rect(fill = 'white'),
                                     axis.title = element_text(size = rel(1.4)),
                                     axis.text = element_text(size = rel(1.4)))

  ### CAGE features
  feats_to_add <- mcols(proms)[, c("external_gene_name", "ensembl_gene_id", "nr_ctss",
                                   "tpm", "tpm.dominant_ctss", "interquantile_width")]
  colnames(feats_to_add)[1] <- "gene"

  ## make sure promoter matrix doesnt have "tpm", "iq" columns
  cols <- c("tpm", "tpm.dominant_ctss", "interquantile_width")
  p <- which(colnames(promoter_matrix) %in% cols)
  promoter_matrix <- promoter_matrix[, -p]
  promoter_matrix <- merge(promoter_matrix, feats_to_add, by="gene")
  rm(feats_to_add)

  # I have noticed that i have some suplicated entries. Will get rid of those
  promoter_matrix <- promoter_matrix[!duplicated(promoter_matrix$ensembl_gene_id), ]
  promoter_matrix <- as.data.frame(promoter_matrix)
  promoter_matrix <- promoter_matrix %>% group_by(gene) %>% arrange(gene, desc(tpm)) %>% top_n(1, wt = "tpm")

  # CAGE features are not norm distributed
  n <- ncol(promoter_matrix)
  to_plot <- promoter_matrix[, c("gene", "tpm", "tpm.dominant_ctss", "interquantile_width")]
  to_plot$sample <- "background"
  to_plot$sample[to_plot$gene %in% sample] <- "sample"
  to_plot <- as.data.frame(to_plot)
  to_plot_m <- melt(to_plot, id.vars=c("gene", "sample"))

  cage_plot <- ggplot(to_plot_m, aes(variable, y=value, fill = sample)) +
    geom_violin() + scale_y_log10() + geom_boxplot(aes(fill = sample), width = 0.1, position = position_dodge(width = 0.9), outlier.color = NA) +
    ggtitle("Distributions of CAGE features for \nSample and the background")+
    theme_fivethirtyeight() + scale_fill_manual(values=wes_palette(n=2, name="Darjeeling1")[2:1]) +
    theme(plot.background = element_rect(fill = 'white'),
          panel.background = element_rect(fill = 'white'),
          legend.background = element_rect(fill = 'white'),
          axis.title = element_text(size = rel(1.3)),
          axis.text = element_text(size = rel(1.1)))

  cage_sig <- c(tpm=wilcox.test(filter(to_plot, sample=="sample")$tpm, filter(to_plot, sample=="background")$tpm)$p.val,
                tpm_dominant=wilcox.test(filter(to_plot, sample=="sample")$tpm.dominant_ctss, filter(to_plot, sample=="background")$tpm.dominant_ctss)$p.val,
                iq=wilcox.test(filter(to_plot, sample=="sample")$interquantile_width, filter(to_plot, sample=="background")$interquantile_width)$p.val)

  to_plot_list <- split(to_plot, f = to_plot$sample)
  to_plot_list <- lapply(to_plot_list, function(x){
    apply(x[, 2:4], 2, mean)
  })
  cage_fold_enrich <- to_plot_list$sample / to_plot_list$background

  cage_sig <- data.frame(feature = names(cage_sig),
                         fold_enrich = cage_fold_enrich,
                        pval = cage_sig)
  cage_sig$class <- ifelse(cage_sig$fold_enrich > 1, "over", "under")
    ###### all plots #######
  if(save_plot == TRUE){
    pdf(paste0(label, ".pdf"), useDingbats = FALSE)
    plot(overTF_plot)
    plot(underTF_plot)
    plot(cage_plot)
    dev.off()
  }else{
    plot(overTF_plot)
    plot(underTF_plot)
    plot(cage_plot)
  }

  ######### significant features #######
  signif <- rbind(cbind(overTF, class="over"),
                  cbind(underTF, class="under"))
  if(all==FALSE){
    signif <- signif[signif$pval < 0.1, ]
  }

  return(list(signif, cage_sig))

}
