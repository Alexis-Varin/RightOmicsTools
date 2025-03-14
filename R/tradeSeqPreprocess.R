#' @title Remove genes with very low expression for \pkg{tradeSeq} analysis
#'
#' @description The functions \code{\link[tradeSeq]{fitGAM}} and \code{\link[tradeSeq]{evaluateK}}, used to analyze gene expression along trajectories, usually require both considerable calculation power (access to computer clusters) and time (up to days) to complete. While one might assume that computing time is uniform across genes, it appears however that model fitting takes a few seconds or less for the vast majority of genes, and that only a small fraction disproportionately affects calculation time, with model fitting taking up to hours or even days for a single gene. Further examination shows that all these genes are among the lowest expressed, and that there is an inverse correlation between the number of cells with a non-zero count for a particular gene and model fitting time, which becomes exponentially longer as the number of expressing cells tends towards zero. Based on these observations, this function was built to considerably speed up computing time / lower computing power requirements while retaining as much information as possible. It filters out genes below a certain expression threshold (10 expressing cells by default), effectively removing genes whose expression is likely just noise (and thus the cause of lengthy and difficult model fitting), rather than a signal of biological significance. This is especially common when cells are subset from a larger object for trajectory inference, and the expression matrix contains many zero or near-zero expression genes, as these would be primarily expressed in excluded cells. Additionally, it offers a second method of filtering, if computing resources are still a limiting factor and further reduction in the number of genes is needed, by retaining only the top genes according to deviance, which provides a better information estimate than selecting variable genes, as the latter tends to still include a few genes with near-zero expressing cells. Finally, if \code{\link[tradeSeq]{fitGAM}} and \code{\link[tradeSeq]{evaluateK}} will be run with a condition, it also allows filtering based on the specific metadata that will be used.
#'
#' @param sds A \pkg{SingleCellExperiment} object containing the pseudotime values of one or several lineages, computed using \code{\link[slingshot]{slingshot}} (usually, the input object to \code{\link[tradeSeq]{fitGAM}}).
#' @param min.cells Numeric. The minimum number of cells with 1 count or more to keep a gene.
#' @param nb.genes Numeric. The maximum number of genes to keep according to deviance, as computed by \code{scry{devianceFeatureSelection}} internally. If \code{NULL}, all genes are kept.
#' @param condition Character or Factor. Either the name of a metadata present in \code{sds} (for example, 'treatment', 'disease', etc) to subset the object and separately remove lowly expressed genes, or the identities, as character or factor, of length equal to the number of cells in \code{sds}.
#' @param plot.genes Logical. If \code{TRUE}, a bar plot showing the number of genes removed according to the number of cells with at least 1 count or more is displayed.
#' @param output.data Logical. If \code{TRUE}, the function will return a \code{list} containing the filtered \pkg{SingleCellExperiment} object, a \code{data.frame} object with the number of cells with at least 1 count or more for each removed gene, and a \code{data.frame} object with the data used to build the bar plot.
#'
#' @return A filtered \pkg{SingleCellExperiment} object, or a \code{list} containing the filtered \pkg{SingleCellExperiment} object, a \code{data.frame} object with the number of cells with at least 1 count or more for each removed gene, and a \code{data.frame} object with the data used to build the bar plot.
#'
#' @import SingleCellExperiment
#' @import slingshot
#' @import scry
#' @import ggplot2
#' @import scales
#' @importFrom stats aggregate
#' @importFrom SummarizedExperiment assays
#' @export

tradeSeqPreprocess = function(sds,
                              min.cells = 10,
                              nb.genes = NULL,
                              condition = NULL,
                              plot.genes = TRUE,
                              output.data = FALSE) {

  nb.cells = NULL

  if (is.character(condition) || is.factor(condition)) {
    if (length(condition) == 1) {
      condition = data.frame(colData(sds)[ , condition, drop = FALSE])
    }
    else {
      condition = data.frame(condition = condition, row.names = rownames(colData(sds)))
    }
    condition.cells = list()
    for (i in 1:length(as.character(unique(condition[ , 1])))) {
      condition.cells[[i]] = rownames(condition[which(condition[ , 1] == as.character(unique(condition[ , 1]))[i]), , drop = FALSE])
    }
  }

  mat.to.subset = list(assays(sds)$counts[rowSums(assays(sds)$counts > 0) < min.cells, ])
  gc(verbose = FALSE)
  sum.mat = data.frame(nb.cells = rowSums(mat.to.subset[[1]] > 0), condition = rep("Global", nrow(mat.to.subset[[1]])))
  if(is.data.frame(condition)) {
    for (i in 2:length(as.character(unique(condition[ , 1])))) {
      mat.to.subset[[i]] = assays(sds)$counts[ , which(colnames(assays(sds)$counts) %in% condition.cells[[i]]), drop = FALSE][rowSums(assays(sds)$counts[ , which(colnames(assays(sds)$counts) %in% condition.cells[[i]]), drop = FALSE] > 0) < min.cells, , drop = FALSE]
      mat.to.subset[[i]] = mat.to.subset[[i]][which(!rownames(mat.to.subset[[i]]) %in% unique(unlist(lapply(mat.to.subset[1:(i-1)], rownames)))), , drop = FALSE]
      gc(verbose = FALSE)
      sum.mat = rbind(sum.mat, data.frame(nb.cells = rowSums(mat.to.subset[[i]] > 0), condition = rep(as.character(unique(condition[ , 1]))[i], nrow(mat.to.subset[[i]]))))
    }
  }
  genes.to.subset = unique(unlist(lapply(mat.to.subset, rownames)))
  sds2 = sds[setdiff(rownames(assays(sds)$counts), genes.to.subset), ]
  message("Removed ", length(genes.to.subset), " genes, ", nrow(assays(sds2)$counts), " remaining (", round(nrow(assays(sds2)$counts)/nrow(assays(sds)$counts)*100, 2), "%)")

  stats.on.mat = as.data.frame(table(sum.mat$nb.cells, sum.mat$condition))
  colnames(stats.on.mat) = c("nb.cells", "condition", "nb.genes")

  if (isTRUE(plot.genes)) {
    print(ggplot(data = stats.on.mat, aes(x = nb.cells, y = nb.genes, fill = condition)) +
            geom_bar(stat = "identity", col = "transparent", position = position_stack(reverse = T)) +
            scale_fill_manual(values = if (is.data.frame(condition)) if (isTRUE(grepl("Global", sum.mat$condition))) c("grey60",hue_pal()(length(unique(sum.mat$condition))-1)) else hue_pal()(length(unique(sum.mat$condition))) else "grey60") +
            labs(x = paste0("Number of cells with 1 count or more per gene"), y = "Number of genes removed", fill = "Identity") +
            geom_text(data = aggregate(nb.genes ~ nb.cells, data = stats.on.mat, sum), aes(y = nb.genes, label = nb.genes), vjust = -0.5) +
            annotate("text", x = max(as.numeric(aggregate(nb.genes ~ nb.cells, data = stats.on.mat, sum)$nb.cells)), y = max(aggregate(nb.genes ~ nb.cells, data = stats.on.mat, sum)$nb.genes), label = paste0("n = ", length(genes.to.subset), " genes removed"), vjust = -0.5, hjust = 1) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
            theme_classic() +
            theme(legend.position = ifelse(is.data.frame(condition), "right", "none")))
  }

  if (is.numeric(nb.genes)) {
    if (nb.genes > nrow(assays(sds2)$counts)) {
      stop("The number of genes to keep is greater than the number of genes remaining")
    }
    sds2 = sds2[names(rowData(devianceFeatureSelection(sds2, nkeep = nb.genes))$binomial_deviance), ]
    message("Kept top ", nrow(assays(sds2)$counts), " genes according to deviance")
    gc(verbose = FALSE)
  }
  if (!is.data.frame(condition)) {
    stats.on.mat = stats.on.mat[, -2, drop = FALSE]
    sum.mat = sum.mat[, -2, drop = FALSE]
  }

  if (isTRUE(output.data)) {
    return(list("slingshot object" = sds2, "Removed genes and number of cells" = sum.mat, "Barplot data" = stats.on.mat))
  }
  else {
    return(sds2)
  }
}
