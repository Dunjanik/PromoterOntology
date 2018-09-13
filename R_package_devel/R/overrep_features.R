overrep_features <- function(prom_matrix, sample,
                             features = c(2:14, 18:19), enrich = "over"){

  sample_matrix <- prom_matrix[prom_matrix$gene %in% sample, features]
  prom_matrix_s <- prom_matrix[, features]
  if (enrich == "over"){
    pval <- phyper(sum(sample_matrix[, 1]) - 1,   # hits in a sample
                   sum(prom_matrix_s[, 1]),             # hits in the universe
                   (nrow(prom_matrix_s) - sum(prom_matrix_s[, 1])), # universe size minus hits
                   nrow(sample_matrix),                  # sample size
                   lower.tail = FALSE)
    message("Calculating overenrichment")

    geneRatio <- sum(sample_matrix[, 1]) / nrow(sample_matrix)
    count <- sum(sample_matrix[, 1])

    for(i in 2:ncol(sample_matrix)){
      pval <- c(pval, phyper(sum(sample_matrix[, i]) - 1,
                             sum(prom_matrix_s[, i]),
                             (nrow(prom_matrix_s) - sum(prom_matrix_s[, i])),
                             nrow(sample_matrix),
                             lower.tail = FALSE) )

      geneRatio <- c(geneRatio,
                     sum(sample_matrix[, i]) / nrow(sample_matrix))

      count <- c(count,
                 sum(sample_matrix[, i]))

    }
  } else { if (enrich == "under"){
    pval <- phyper(sum(sample_matrix[, 1]),
                   sum(prom_matrix_s[, 1]),
                   (nrow(prom_matrix_s) - sum(prom_matrix_s[, 1])),
                   nrow(sample_matrix),
                   lower.tail = TRUE)
    message("Calculating underenrichment")

    geneRatio <- sum(sample_matrix[, 1]) / nrow(sample_matrix)
    count <- sum(sample_matrix[, 1])
    for(i in 2:ncol(sample_matrix)){
      pval <- c(pval, phyper(sum(sample_matrix[, i]),
                             sum(prom_matrix_s[, i]),
                             (nrow(prom_matrix_s) - sum(prom_matrix_s[, i])),
                             nrow(sample_matrix),
                             lower.tail = TRUE) )

      geneRatio <- c(geneRatio,
                     sum(sample_matrix[, i]) / nrow(sample_matrix))

      count <- c(count,
                 sum(sample_matrix[, i]))

    }
  } else stop("Type of enrichment can be either over or under!") }

  result <- cbind(pval, geneRatio, count)
  result <- as.data.frame(result)
  rownames(result) <- colnames(sample_matrix)

  return(result)

}
