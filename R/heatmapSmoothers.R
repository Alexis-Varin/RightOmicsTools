#' @title Heatmap of average fitted GAM smoothers along pseudotime
#'
#' @description This function generates a heatmap to visualize gene fitted GAM smoothers, computed using \code{\link[tradeSeq]{fitGAM}}, along pseudotime.
#'
#' @param sds A \pkg{SingleCellExperiment} object containing the pseudotime values of one or several lineages, computed using \code{\link[slingshot]{slingshot}} (usually, the input object to \code{\link[tradeSeq]{fitGAM}}).
#' @param models A \pkg{SingleCellExperiment} object containing the fitted GAM smoothers, computed using \code{\link[tradeSeq]{fitGAM}}, with or without \code{conditions} provided.
#' @param predictSmooth.df A \code{data.frame} object containing the average fitted GAM smoothers, computed using \code{\link[tradeSeq]{predictSmooth}} with \code{nPoints}, \code{lineages} and \code{tidy} = \code{TRUE}. If \code{NULL}, average fitted GAM smoothers will be internally computed.
#' @param genes Character. The names of one or several genes to plot the average fitted GAM smoothers from.
#' @param lineages Numeric. The indices of the lineages (for example, c(1, 5), 2:4 etc) to plot the average fitted GAM smoothers and the pseudotime values from. Please note that for gene clustering, if you would like to use a specific lineage as reference, you need to set it as the first element (for example, if you have six lineages and you would like to cluster genes based on lineage 5, you would input \code{lineages} = c(5, 1:4, 6)). This is due to a limitation in how \pkg{ComplexHeatmap} clusters the rows of concatenated \code{\link[ComplexHeatmap]{HeatmapList}} objects (always based on the first element).
#' @param lineages.to.remove Numeric. The lineages to exclude from \code{\link[scales]{rescale}}. Useful if you want to remove the influence of other lineages in a comparison.
#' @param conditions Character. The names of one or several \code{conditions} identities to select. If \code{NULL}, all identities are used, and each unique condition will be plotted, for each of the \code{lineages} provided (for example, if you have four lineages and three conditions, twelve heatmaps will be plotted). Please note that for gene clustering, if you would like to use a specific condition as reference, you need to set it as the first element (for example, if you have three unique conditions (for example, 'control', 'IFNg' and 'TGFb') and you would like to cluster genes based on 'IFNg', you would input \code{conditions} = c('IFNg', 'control', 'TGFb')). This is due to a limitation in how \pkg{ComplexHeatmap} clusters the rows of concatenated \code{\link[ComplexHeatmap]{HeatmapList}} objects (always based on the first element). Ignored if the \code{models} object was computed using \code{\link[tradeSeq]{fitGAM}} without \code{conditions}. Please note that if the \code{models} object was computed using \code{\link[tradeSeq]{fitGAM}} with \code{conditions}, it is not possible to plot global lineages (without \code{conditions}), you would need to compute a new \code{models} object using \code{\link[tradeSeq]{fitGAM}} without \code{conditions}. This is due to a limitation in how \code{\link[tradeSeq]{predictSmooth}} returns average fitted GAM smoothers (if the \code{models} object was computed with \code{conditions}, the function will always mean the fitted GAM smoothers for each lineage and each condition independently).
#' @param clusters Character or Factor. Either the name of a metadata present in the \code{sds} object (for example, 'annotations', 'seurat_clusters', etc) to plot the cell densities from, or the identities, as character or factor, of length equal to the number of cells in the \code{sds} object.
#' @param nPoints Numeric. (from \code{\link[tradeSeq]{predictSmooth}} documentation) The number of points used to create the grid along the smoother for each lineage.
#' @param branch.points Numeric, Character or List. Branching points may be shown on the cell density plot and the heatmap, to partition the lineages and help visualize differences. Either one or more pseudotime values to plot the branching points at, and/or one or several values containing 'knot' followed by a number (for example, 'knot2', 'knot4' etc), which correspond to the knots (\code{k}) input in \code{\link[tradeSeq]{fitGAM}} and which divide each lineage into segments; the function will extract the pseudotime values from each knot number and plot the branching points. Mixing any of the two options is possible (for example, c(4.3, 7, 'knot4')). You may also provide a \code{list} containing any of the two options, and named after one or several identities of a single metadata among 'lineages' or 'conditions' (for example, if \code{conditions} = c('Control', 'Treated'), list('Control' = c(2, 'knot5')) which will only plot branching points on Control cell density plots and heatmaps). For \code{lineages}, the \code{list} names need to be provided as 'Lineage' (with capital L) followed by the index (for example, list('Lineage1' = c(4.3, 'knot4'), 'Lineage2' = c('knot2', 7, 9)), which will only plot branching points on lineages 1 and 2 cell density plots and heatmaps). Finally, you may also provide identities corresponding to several metadata at the same time by pasting them together with '_' (for example, list('Lineage1_Treated' = c(1, 'knot1', 12)), which will only plot branching points on lineage 1 and Treated cell density plots and heatmaps).
#' @param pseudotime.zoom Numeric, Character or List. Oftentimes you may want to zoom into a specific region of a lineage to better visualize the differences. Either two pseudotime values to zoom into, and/or 'min' or 'max' to zoom respectively into the minimum or maximum pseudotime values of the lineage, and/or 'knot' followed by a number (for example, 'knot2', 'knot4' etc), which correspond to the knots (\code{k}) input in \code{\link[tradeSeq]{fitGAM}} and which divide each lineage into segments; the function will extract the pseudotime values from each knot number and zoom into, and/or 'knot' followed by a number, a + or - and another number, which will set the zoom around knots by adding or subtracting from the pseudotime value of a knot (for example, c('knot3-5', 'knot4+2'), which corresponds to the pseudotime value of the third knot minus 5, and to the pseudotime value of the fourth knot plus 2). Mixing any two of the four options is possible (for example, c(5.5, 'knot4'), or c('knot3-1', 'max')). You may also provide a \code{list} containing any two of the four options, and named after one or several identities of a single metadata among 'lineages' or 'conditions' (for example, if \code{conditions} = c('Control', 'Treated'), list('Control' = c(2, 'knot5')) which will zoom on Control cell density plots and heatmaps). For \code{lineages}, the \code{list} names need to be provided as 'Lineage' (with capital L) followed by the index (for example, list('Lineage1' = c('min', 'knot2+3'), 'Lineage2' = c('knot2', 7)), which will only zoom on lineages 1 and 2 cell density plots and heatmaps). Finally, you may also provide identities corresponding to several metadata at the same time by pasting them together with '_' (for example, list('Lineage1_Treated' = c('knot1+1', 'max')), which will only zoom on lineage 1 and Treated cell density plots and heatmaps).
#' @param rescale Logical. If \code{TRUE}, average fitted GAM smoothers will be adjusted using \code{\link[scales]{rescale}} between the first numerical value of \code{rescale.range} (lowest value) and the second numerical value (highest value). This is different than \code{\link[base]{scale}} as this doesn't compare values to any mean or standard deviation and is therefore not a Z-score, it only refits each value (independently for each gene) in order to visualize all \code{genes} in the same dimension regardless of their differences in fitted GAM smoothers.
#' @param rescale.range Numeric. The minimum and maximum values to resize the average fitted GAM smoothers and internally passed to \code{\link[scales]{rescale}}. These values are arbitrary and will not change the visualization, only the values in the legend. Ignored if \code{rescale} = \code{FALSE}.
#' @param kmeans.repeats Numeric. The number of runs to get a consensus K-means clustering. Ignored if \code{genes.kmeans} = 1.
#' @param cluster.genes Logical or Function. If \code{TRUE}, the function will cluster the \code{genes}. You may also pass an \code{hclust} or \code{dendrogram} object which contains clustering.
#' @param genes.kmeans Numeric. The number of slices to use for gene K-means clustering.
#' @param genes.kmeans.numbers.size Numeric. The font size of the gene K-means slice numbers. Set to 0 to remove them.
#' @param outer.border Logical. If \code{TRUE}, the function will display an outer border around each heatmap.
#' @param pseudotime.type Character. Determines the \code{pseudotime.colors} range scale of each lineage: either 'independent', where each lineage's \code{pseudotime.colors} range is independent of the others (i.e., minimum and maximum values are identical for each lineage, regardless of trajectory length differences), or 'relative', where each lineage's \code{pseudotime.colors} range is scaled relative to the highest pseudotime value (i.e., the longest trajectory, reflecting differences in trajectory lengths).
#' @param density.type Character. Determines the cell density scale (height) of each lineage: either 'independent,' where each lineage's cell density is independent of the others (i.e., the maximum density height is identical for each lineage, regardless of differences in cell counts contributing to a trajectory), or 'relative,' where each lineage's density scale is relative to the highest cell density (i.e., the highest cell counts, reflecting differences in the number of cells contributing to a trajectory).
#' @param show.pseudotime Logical. If \code{TRUE}, the pseudotime range will be shown above the heatmap.
#' @param show.density Logical. If \code{TRUE}, the cell density plot will be shown.
#' @param show.density.pseudotime.values Logical. If \code{TRUE}, the pseudotime values will be shown on the cell density plot.
#' @param data.colors Character. Either two color names, corresponding to the lowest and highest values in the average fitted GAM smoothers and internally passed to \code{\link[colorRamp2]{colorRamp2}}, or the name of a palette and internally passed to \code{hcl_palette} in \code{\link[colorRamp2]{colorRamp2}} (such as 'Inferno', 'Berlin', 'Viridis' etc, check \code{\link[grDevices]{hcl.pals}} for all palettes available).
#' @param pseudotime.colors Character. Either the name of a palette and internally passed to \code{hcl_palette} in \code{\link[colorRamp2]{colorRamp2}} (such as 'Inferno', 'Berlin', 'Viridis' etc, check \code{\link[grDevices]{hcl.pals}} for all palettes available), or two or more color names, corresponding to the lowest and highest pseudotime values, and additional color names beyond two will either be spaced at regular intervals if unnamed (for example, c('blue', 'green', 'yellow', 'red') will be spaced at 0%, 33%, 66% and 100%) or if named with a pseudotime value (for example, c('orange', "12" = 'white', 'royalblue')) or if named with 'knot' followed by a number (for example, 'knot2', 'knot4' etc), which correspond to the knots (\code{k}) input in \code{\link[tradeSeq]{fitGAM}} and which divide each lineage into segments; the function will extract the pseudotime values from each knot number (for example, c('firebrick', 'knot3' = 'lightgrey', 'gold')). Note that the first and last color names do not need to be named as they always correspond to the lowest and highest pseudotime values.
#' @param clusters.colors Character. The color names for each identity in \code{clusters}. If \code{NULL}, uses \pkg{Seurat}'s default colors. Ignored if \code{show.density} = \code{FALSE}.
#' @param genes.names.size Numeric. The font size of the gene names. Set to 0 to remove them.
#' @param genes.names.style Character. The font face of the gene names. The \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7494048/}{Gene nomenclature} used by almost all scientific journals require that gene names are italicized, therefore the parameter is by default set to 'italic'. Use 'plain' to revert back to regular font face.
#' @param show.density.legend Logical. If \code{TRUE}, the cell density plot legend will be shown.
#' @param density.legend.name Character. The name of the cell density plot legend.
#' @param heatmap.names Character. You may provide custom names for each heatmap, instead of 'Lineage' followed by the index of the \code{lineages} and \code{conditions} provided.
#' @param heatmap.width Numeric. The width of each heatmap.
#' @param heatmap.height Numeric. The height of each heatmap.
#' @param density.height Numeric. The height of the cell density plot.
#' @param raster Logical. (from \code{\link[ComplexHeatmap]{Heatmap}} documentation) If \code{TRUE}, the function will render the heatmap body as a raster image. It helps to reduce file size when the matrix is huge.
#' @param return.grob Logical. If \code{TRUE}, the function will return a \code{grob} (\code{gTree} object, printable using \code{\link[grid]{grid.draw}}). This allows to capture the cell density plot (which is a \pkg{ggplot2} object) and the heatmap (a \code{\link[ComplexHeatmap]{HeatmapList}} object) as a single object which can be further customized using \pkg{grid} functions or included with other plots in complex layouts using \code{\link[cowplot]{plot_grid}}. If \code{FALSE}, the function will return a \code{\link[ComplexHeatmap]{HeatmapList}} object, which can then be used to extract various parameters, such as the row clustering values using \code{\link[ComplexHeatmap]{row_order}}. Please note that this only affects the returned object when it is assigned to a variable, such as \code{> htmp <- heatmapSmoothers(...)}, and not when the function is called without assignment, such as \code{> heatmapSmoothers(...)}.
#' @param ... Additional arguments to be passed to \code{\link[ComplexHeatmap]{Heatmap}}, such as \code{show_parent_dend_line}, \code{clustering_method_rows}, etc, accepts any parameter that wasn't already internally passed to \code{\link[ComplexHeatmap]{Heatmap}} (for example, \code{outer.border} sets the \code{border} parameter of \code{\link[ComplexHeatmap]{Heatmap}}, so you will get an error if you try to pass the \code{border} parameter in \code{\link[RightOmicsTools]{heatmapSmoothers}}).
#'
#' @return A \code{grob} (\code{gTree} object, printable using \code{\link[grid]{grid.draw}}), or a \code{\link[ComplexHeatmap]{HeatmapList}} object.
#'
#' @import SingleCellExperiment
#' @import slingshot
#' @import tradeSeq
#' @import data.table
#' @import ComplexHeatmap
#' @import colorRamp2
#' @import ggplot2
#' @import scales
#' @importFrom graphics plot.new
#' @importFrom stats na.omit
#' @importFrom S4Vectors metadata
#' @export

heatmapSmoothers = function(sds,
                            models,
                            predictSmooth.df = NULL,
                            genes = if (is.null(predictSmooth.df)) NA else unique(predictSmooth.df$gene),
                            lineages,
                            lineages.to.remove = NULL,
                            conditions = NULL,
                            clusters = NULL,
                            nPoints = 1000,
                            branch.points = NULL,
                            pseudotime.zoom = NULL,
                            rescale = TRUE,
                            rescale.range = c(0, 3),
                            kmeans.repeats = 100,
                            cluster.genes = TRUE,
                            genes.kmeans = 1,
                            genes.kmeans.numbers.size = 11,
                            outer.border = TRUE,
                            pseudotime.type = "relative",
                            density.type = "relative",
                            show.pseudotime = TRUE,
                            show.density = TRUE,
                            show.density.pseudotime.values = TRUE,
                            data.colors = "Inferno",
                            pseudotime.colors = "Viridis",
                            clusters.colors = NULL,
                            genes.names.size = 6,
                            genes.names.style = "italic",
                            show.density.legend = TRUE,
                            density.legend.name = "Clusters",
                            heatmap.names = NULL,
                            heatmap.width = 3,
                            heatmap.height = length(genes),
                            density.height = 1,
                            raster = TRUE,
                            return.grob = TRUE,
                            ...) {

  gene = lineage = xx = yy = count = heatmap.named = NULL

  if (any(is.na(genes))) {
    stop("Please provide the gene names to plot the smoothed expression from")
  }

  if (is.null(clusters) & isTRUE(show.density)) {
    stop("Please provide the clusters, as metadata name or identities, to plot the cell densities from")
  }

  knots = ceiling(unname(metadata(models)$tradeSeq$knots)*100)/100

  data = data.frame(slingPseudotime(sds))
  data = ceiling(data*100)/100
  data = cbind(data.frame(slingReducedDim(sds)), data)
  if (length(clusters) == 1) {
    data$clusters = as.character(colData(sds)[ , clusters])
  }
  else {
    data$clusters = as.character(clusters)
  }
  colnames(data) = gsub("Lineage", "", colnames(data))

  if (!is.data.frame(predictSmooth.df)) {
    predictSmooth.df = predictSmooth(models, gene = genes, nPoints = nPoints, tidy = TRUE)
  }
  if (length(setdiff(lineages.to.remove, lineages)) > 0) {
    predictSmooth.df = predictSmooth.df[predictSmooth.df$lineage != setdiff(lineages.to.remove, lineages), , drop = FALSE]
  }

  if (isFALSE("condition" %in% colnames(predictSmooth.df))) {
    predictSmooth.df$condition = ""
    nocond = TRUE
  }
  else {
    nocond = FALSE
    if (is.null(conditions)) {
      conditions = unique(predictSmooth.df$condition)
      message("Defaulting to '", conditions[1], "' as reference for heatmap clustering")
    }
    lineages = paste(rep(lineages, each = length(conditions)), conditions, sep = "_")
    data$conditions = as.character(colData(models)$tradeSeq$conditions)
  }

  if (!is.null(pseudotime.zoom)) {
    if (!is.list(pseudotime.zoom)) {
      pseudotime.zoom = rep(list(pseudotime.zoom), length(lineages))
      names(pseudotime.zoom) = pseudotime.names = paste0("Lineage", lineages)
    }
    else if (isFALSE(nocond)) {
      pseudotime.names = grep("Lineage", names(pseudotime.zoom), value = T)
      if (length(unlist(strsplit(names(pseudotime.zoom), "_"))) - length(pseudotime.names) == 0) {
        pseudotime.zoom = rep(pseudotime.zoom, each = length(conditions))
        names(pseudotime.zoom) = pseudotime.names = paste(rep(pseudotime.names, each = length(conditions)), conditions, sep = "_")
      }
      else if (length(pseudotime.names) == 0) {
        pseudotime.names = paste(rep(paste0("Lineage", gsub("_.*","",lineages)), each = length(names(pseudotime.zoom))), names(pseudotime.zoom), sep = "_")
        pseudotime.zoom = rep(pseudotime.zoom, each = length(gsub("_.*","",lineages)))
        names(pseudotime.zoom) = pseudotime.names
      }
    }
    else {
      pseudotime.names = names(pseudotime.zoom)
    }
    pseudotime.zoom = lapply(seq_along(lineages), function(i) {
      if (paste0("Lineage", lineages)[i] %in% names(pseudotime.zoom)) {
        pseudotime.zoom[[paste0("Lineage", lineages)[i]]]
      }
      else {
        NULL
      }
    })
    subsetted.pseudo = list()
    for (i in 1:length(lineages)) {
      subsetted.pseudo[[i]] = vector("numeric", 2)
      if (paste0("Lineage", lineages[i]) %in% pseudotime.names) {
        for (j in 1:length(pseudotime.zoom[[i]])) {
          if (isTRUE(grepl("max",pseudotime.zoom[[i]][j]))) {
            if (isFALSE(nocond)) {
              subsetted.pseudo[[i]][j] = max(predictSmooth.df[predictSmooth.df$lineage %in% gsub("_.*","",lineages[i]) & predictSmooth.df$condition %in% sub(".*_","",lineages[i]), "time"], na.rm = T)
            }
            else {
              subsetted.pseudo[[i]][j] = max(predictSmooth.df[predictSmooth.df$lineage %in% lineages[i], "time"], na.rm = T)
            }
          }
          else if (isTRUE(grepl("min",pseudotime.zoom[[i]][j]))) {
            if (isFALSE(nocond)) {
              subsetted.pseudo[[i]][j] = min(predictSmooth.df[predictSmooth.df$lineage %in% gsub("_.*","",lineages[i]) & predictSmooth.df$condition %in% sub(".*_","",lineages[i]), "time"], na.rm = T)
            }
            else {
              subsetted.pseudo[[i]][j] = min(predictSmooth.df[predictSmooth.df$lineage %in% lineages[i], "time"], na.rm = T)
            }
          }
          else if (isTRUE(grepl("knot",pseudotime.zoom[[i]][j]))) {
            pseudotime.zoom[[i]][j] = gsub("knot", "", pseudotime.zoom[[i]][j])
            if (nchar(pseudotime.zoom[[i]][j]) > 1) {
              pseudotime.operation = sub("^[0-9]", "", pseudotime.zoom[[i]][j])
              pseudotime.zoom[[i]][j] = sub("^([0-9]).*","\\1", pseudotime.zoom[[i]][j])
              pseudotime.zoom[[i]][j] = knots[as.numeric(pseudotime.zoom[[i]][j])]
              pseudotime.zoom[[i]][j] = paste(pseudotime.zoom[[i]][j], pseudotime.operation, sep = "")
              subsetted.pseudo[[i]][j] = as.numeric(eval(parse(text = pseudotime.zoom[[i]][j])))
            }
            else {
              subsetted.pseudo[[i]][j] = knots[as.numeric(pseudotime.zoom[[i]][j])]
            }
          }
          else {
            subsetted.pseudo[[i]][j] = as.numeric(pseudotime.zoom[[i]][j])
          }
        }
      }
    }
  }
  else {
    pseudotime.names = paste0("Lineage", lineages)
  }

  if (isTRUE(show.pseudotime)) {
    pseudo.col = list()
    for (i in 1:length(lineages)) {
      if (isFALSE(nocond)) {
        min.lin = min(predictSmooth.df[predictSmooth.df$lineage %in% gsub("_.*","",lineages[i]) & predictSmooth.df$condition %in% sub(".*_","",lineages[i]), "time"], na.rm = T)
        max.lin = max(predictSmooth.df[predictSmooth.df$lineage %in% gsub("_.*","",lineages[i]) & predictSmooth.df$condition %in% sub(".*_","",lineages[i]), "time"], na.rm = T)
      }
      else {
        min.lin = min(predictSmooth.df[predictSmooth.df$lineage %in% lineages[i], "time"], na.rm = T)
        max.lin = max(predictSmooth.df[predictSmooth.df$lineage %in% lineages[i], "time"], na.rm = T)
      }
      if (pseudotime.type == "relative") {
        min.breaks = min(data.frame(slingPseudotime(sds)), na.rm = T)
        max.breaks = max(data.frame(slingPseudotime(sds)), na.rm = T)
      }
      else {
        min.breaks = min.lin
        max.breaks = max.lin
      }
      if (length(pseudotime.colors) == 1) {
        pseudo.col[[i]] = colorRamp2(breaks = c(min.breaks, max.breaks),
                                     hcl_palette = tolower(pseudotime.colors))(seq(min.lin, max.lin, length.out = nPoints))
      }
      else {
        pseudotime.colors = unlist(lapply(seq_along(pseudotime.colors), function(col) {
          if (col == 1) {
            col = setNames(pseudotime.colors[col], min.breaks)
          }
          else if (col == length(pseudotime.colors)) {
            col = setNames(pseudotime.colors[col], max.breaks)
          }
          else if (isTRUE(grepl("knot",names(pseudotime.colors[col])))) {
            col = setNames(pseudotime.colors[col], knots[as.numeric(gsub("knot", "", names(pseudotime.colors[col])))])
          }
          else if (is.null(names(pseudotime.colors[col])) || names(pseudotime.colors[col]) == "") {
            col = setNames(pseudotime.colors[col], as.character((max.breaks-min.breaks)*(col-1)/(length(pseudotime.colors)-1)))
          }
          else {
            col = setNames(pseudotime.colors[col], names(pseudotime.colors[col]))
          }
          return(col)
        }))
        pseudo.col[[i]] = colorRamp2(breaks = as.numeric(names(pseudotime.colors)),
                                     colors = pseudotime.colors)(seq(min.lin, max.lin, length.out = nPoints))
      }
      names(pseudo.col[[i]]) = predictSmooth.df[predictSmooth.df$lineage %in% gsub("_.*","",lineages[i]) & predictSmooth.df$condition %in% sub(".*_","",lineages[i]) & predictSmooth.df$gene == genes[1], "time"]
      if (is.list(pseudotime.zoom) & paste0("Lineage", lineages[i]) %in% pseudotime.names) {
        pseudo.col[[i]] = pseudo.col[[i]][as.numeric(names(pseudo.col[[i]])) > subsetted.pseudo[[i]][1] & as.numeric(names(pseudo.col[[i]])) < subsetted.pseudo[[i]][2]]
      }
    }
  }

  mat = setDF(dcast(setDT(predictSmooth.df),gene~lineage+condition+time, value.var = "yhat"))
  rownames(mat) = mat$gene
  mat$gene = NULL
  mat = log1p(mat)
  mat = mat[match(genes, rownames(mat)), ]

  if (isTRUE(rescale)) {
    mat = t(apply(mat, 1, rescale, to = rescale.range))
  }

  mat.ht = list()
  for (i in 1:length(lineages)) {
    mat.ht[[i]] = mat[ , grepl(paste0("^",lineages[i],"_"), colnames(mat)), drop = FALSE]
    colnames(mat.ht[[i]]) = sapply(strsplit(colnames(mat.ht[[i]]), "_\\s*"), tail, 1)
    if (is.list(pseudotime.zoom) & paste0("Lineage", lineages[i]) %in% pseudotime.names) {
      mat.ht[[i]] = mat.ht[[i]][ , as.numeric(colnames(mat.ht[[i]])) > subsetted.pseudo[[i]][1] & as.numeric(colnames(mat.ht[[i]])) < subsetted.pseudo[[i]][2], drop = FALSE]
    }
  }

  if (is.null(clusters.colors)) {
    clusters.colors = hue_pal()(n = length(levels(as.factor(colData(sds)[ , clusters]))))
    names(clusters.colors) = levels(as.factor(colData(sds)[ , clusters]))
  }
  else if (is.null(names(clusters.colors))) {
    names(clusters.colors) = levels(as.factor(colData(sds)[ , clusters]))
  }

  if (isTRUE(show.density)) {
    p = tmp = list()
    for (i in 1:length(lineages)) {
      if (isFALSE(nocond)) {
        tmp[[i]] = na.omit(data[data$conditions == gsub(".*_","",lineages[i]) , c(gsub("_.*","",lineages[i]), "clusters"), drop = FALSE])
      }
      else {
        tmp[[i]] = na.omit(data[ , c(gsub("_.*","",lineages[i]), "clusters"), drop = FALSE])
      }
      colnames(tmp[[i]]) = c("lineage", "clusters")
      if (is.list(pseudotime.zoom) & paste0("Lineage", lineages[i]) %in% pseudotime.names) {
        tmp[[i]] = tmp[[i]][tmp[[i]]$lineage > subsetted.pseudo[[i]][1] & tmp[[i]]$lineage < subsetted.pseudo[[i]][2], , drop = FALSE]
      }
      p[[i]] = ggplot(data=tmp[[i]], aes(x = lineage)) +
        geom_density(alpha = 0.5, aes(y = after_stat(count), fill = clusters), col = "transparent") +
        geom_density(aes(y = after_stat(count), col = clusters), fill = "transparent", linewidth = 1) +
        theme_bw() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linewidth = 0.2),
              legend.key = element_rect(colour = "black"),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title = element_blank(),
              plot.margin = margin(0, 0, 0, 0)) +
        guides(col = "none", fill ="none") +
        scale_fill_manual(values = clusters.colors[names(clusters.colors) %in% unique(tmp[[i]]$clusters)]) +
        scale_color_manual(values = clusters.colors[names(clusters.colors) %in% unique(tmp[[i]]$clusters)]) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(ceiling(min(tmp[[i]]$lineage, na.rm = TRUE)), max(tmp[[i]]$lineage, na.rm = TRUE), by = ceiling(sqrt(max(tmp[[i]]$lineage, na.rm = TRUE)-min(tmp[[i]]$lineage, na.rm = TRUE)))))
    }

    if (density.type == "relative") {
      max.value = rep(max(unlist(lapply(p, function(x) max(suppressWarnings(layer_scales(x)$y$range$range), na.rm = TRUE)))), length(p))
    }
    if (density.type == "independent") {
      max.value = unlist(lapply(p, function(x) max(suppressWarnings(layer_scales(x)$y$range$range), na.rm = TRUE)))
    }

    ht.split = vector("list", length(lineages))
    if (!is.null(branch.points)) {
      if (!is.list(branch.points)) {
        branch.points = rep(list(branch.points), length(lineages))
        names(branch.points) = branch.names = paste0("Lineage", lineages)
      }
      else if (isFALSE(nocond)) {
        branch.names = grep("Lineage", names(branch.points), value = T)
        if (length(unlist(strsplit(names(branch.points), "_"))) - length(branch.names) == 0) {
          branch.points = rep(branch.points, each = length(conditions))
          names(branch.points) = branch.names = paste(rep(branch.names, each = length(conditions)), conditions, sep = "_")
        }
        else if (length(branch.names) == 0) {
          branch.names = paste(rep(paste0("Lineage", gsub("_.*","",lineages)), each = length(names(branch.points))), names(branch.points), sep = "_")
          branch.points = rep(branch.points, each = length(gsub("_.*","",lineages)))
          names(branch.points) = branch.names
        }
      }
      else {
        branch.names = names(branch.points)
      }
      branch.points = lapply(seq_along(lineages), function(i) {
        if (paste0("Lineage", lineages)[i] %in% names(branch.points)) {
          branch.points[[paste0("Lineage", lineages)[i]]]
        }
        else {
          NULL
        }
      })
      for (i in 1:length(lineages)) {
        if (paste0("Lineage", lineages[i]) %in% branch.names) {
          for (j in 1:length(branch.points[[i]])) {
            if (isTRUE(grepl("knot",branch.points[[i]][j]))) {
              branchpoint = gsub("knot", "", branch.points[[i]][j])
              branchpoint = knots[as.numeric(branchpoint)]
            }
            else {
              branchpoint = as.numeric(branch.points[[i]][j])
            }
            if (is.list(pseudotime.zoom)) {
              branchpoint = ifelse(branchpoint < max(tmp[[i]]$lineage, na.rm = TRUE), branchpoint, ifelse(isTRUE(any(grepl("max", pseudotime.zoom[[i]]))) | pseudotime.zoom[[i]][2] >= max(data[ , gsub("_.*","",lineages[i]), drop = FALSE], na.rm = TRUE) | is.na(pseudotime.zoom[[i]]), max(tmp[[i]]$lineage, na.rm = TRUE), NA))
              branchpoint = ifelse(!is.na(branchpoint), ifelse(branchpoint > min(tmp[[i]]$lineage, na.rm = TRUE), branchpoint, ifelse(isTRUE(any(grepl("min", pseudotime.zoom[[i]]))) | pseudotime.zoom[[i]][2] <= min(data[ , gsub("_.*","",lineages[i]), drop = FALSE], na.rm = TRUE) | is.na(pseudotime.zoom[[i]]), min(tmp[[i]]$lineage, na.rm = TRUE), NA)), NA)
              if (is.na(branchpoint)) {
                warning("A branch point from Lineage ", gsub("_.*","",lineages[i]), " is outside the pseudotime zoom, it will not be displayed")
              }
            }
            else {
              branchpoint = ifelse(branchpoint < max(tmp[[i]]$lineage, na.rm = TRUE), branchpoint, max(tmp[[i]]$lineage, na.rm = TRUE))
              branchpoint = ifelse(branchpoint > min(tmp[[i]]$lineage, na.rm = TRUE), branchpoint, min(tmp[[i]]$lineage, na.rm = TRUE))
            }
            if (!is.na(branchpoint)) {
              if (branchpoint == max(tmp[[i]]$lineage, na.rm = TRUE)) {
                branchpoint.df = data.frame(xx = rep(branchpoint-(0.07*(max(tmp[[i]]$lineage, na.rm = TRUE)-min(tmp[[i]]$lineage, na.rm = TRUE))/17), 5), yy = seq(0, max.value[i]*0.85, length.out = 5))
                branchpoint.df2 = data.frame(xx = branchpoint-(0.65*(max(tmp[[i]]$lineage, na.rm = TRUE)-min(tmp[[i]]$lineage, na.rm = TRUE))/17), yy = max.value[i]*0.85)
                p[[i]] = p[[i]] +
                  geom_path(data=branchpoint.df, aes(x = xx, y = yy), linewidth = 1, col = "black")+
                  geom_point(data=branchpoint.df2, aes(x = xx, y = yy), col = "black", size = 8, shape = 15)+
                  annotate("text", x = branchpoint-(0.65*(max(tmp[[i]]$lineage, na.rm = TRUE)-min(tmp[[i]]$lineage, na.rm = TRUE))/17), y = max.value[i]*0.85,
                           label = ifelse(is.null(names(branch.points[[i]])),"BP",names(branch.points[[i]])[j]),
                           hjust=0.5,
                           size=3,
                           col="white")
              }
              else if (branchpoint == min(tmp[[i]]$lineage, na.rm = TRUE)) {
                branchpoint.df = data.frame(xx = rep(branchpoint+(0.07*(max(tmp[[i]]$lineage, na.rm = TRUE)-min(tmp[[i]]$lineage, na.rm = TRUE))/17), 5), yy = seq(0, max.value[i]*0.85, length.out = 5))
                branchpoint.df2 = data.frame(xx = branchpoint+(0.7*(max(tmp[[i]]$lineage, na.rm = TRUE)-min(tmp[[i]]$lineage, na.rm = TRUE))/17), yy = max.value[i]*0.85)
                p[[i]] = p[[i]] +
                  geom_path(data=branchpoint.df, aes(x = xx, y = yy), linewidth = 1, col = "black")+
                  geom_point(data=branchpoint.df2, aes(x = xx, y = yy), col = "black", size = 8, shape = 15)+
                  annotate("text", x = branchpoint+(0.7*(max(tmp[[i]]$lineage, na.rm = TRUE)-min(tmp[[i]]$lineage, na.rm = TRUE))/17), y = max.value[i]*0.85,
                           label = ifelse(is.null(names(branch.points[[i]])),"BP",names(branch.points[[i]])[j]),
                           hjust=0.5,
                           size=3,
                           col="white")
              }
              else {
                branchpoint.df = data.frame(xx = rep(branchpoint, 5), yy = seq(0, max.value[i]*0.85, length.out = 5))
                branchpoint.df2 = data.frame(xx = branchpoint, yy = max.value[i]*0.85)
                p[[i]] = p[[i]] +
                  geom_path(data=branchpoint.df, aes(x = xx, y = yy), linewidth = 1, col = "black", arrow = arrow(length = unit(0.1, "inches"), ends = "first"))+
                  geom_point(data=branchpoint.df2, aes(x = xx, y = yy), col = "black", size = 8)+
                  annotate("text", x = branchpoint, y = max.value[i]*0.85,
                           label = ifelse(is.null(names(branch.points[[i]])),"BP",names(branch.points[[i]])[j]),
                           hjust=0.5,
                           size=3,
                           col="white")
              }
              if (branchpoint < max(tmp[[i]]$lineage, na.rm = TRUE) & branchpoint > min(tmp[[i]]$lineage, na.rm = TRUE)) {
                ht.split[[i]] = c(ht.split[[i]], round(Position(function(x) x > branchpoint, as.numeric(colnames(mat.ht[[i]])))/ncol(mat.ht[[i]]), 3))
              }
            }
          }
        }
        p[[i]] = p[[i]] + scale_y_continuous(expand = c(0, 0), limits = c(0, max.value[i]*1.05))
      }
    }
    else {
      for (i in 1:length(lineages)) {
        p[[i]] = p[[i]] + scale_y_continuous(expand = c(0, 0), limits = c(0, max.value[i]*1.05))
      }
    }
  }
  if (isFALSE(show.density.pseudotime.values)) {
    for (i in 1:length(lineages)) {
      p[[i]] = p[[i]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
  }

  if (length(data.colors) == 1) {
    data.colors = colorRamp2(breaks = rescale.range, hcl_palette = data.colors)
  }
  else {
    data.colors = colorRamp2(breaks = rescale.range, colors = c(data.colors[1], data.colors[2]))
  }

  density.legend = NULL
  if (isTRUE(show.density.legend) & isTRUE(show.density)) {
    density.legend = Legend(at = names(clusters.colors)[names(clusters.colors) %in% unique(do.call(rbind, tmp)$clusters)],
                            legend_gp = gpar(fill = clusters.colors[names(clusters.colors)[names(clusters.colors) %in% unique(do.call(rbind, tmp)$clusters)]]),
                            title = density.legend.name,
                            gap = unit(0.5, "cm"),
                            border = TRUE)
  }

  ht_opt$ANNOTATION_LEGEND_PADDING = unit(10, "mm")
  ht = NULL
  ht.anno = list()
  for (i in 1:length(lineages)) {
    if (!is.character(heatmap.names)) {
      heatmap.named[i] = ifelse(isFALSE(nocond),paste0("Lineage ", gsub("_.*","",lineages[i]), "\n", gsub(".*_","",lineages[i])),
                                paste0("Lineage ", gsub("_.*","",lineages[i])))
    }
    else {
      heatmap.named[i] = heatmap.names[i]
    }
    if (isTRUE(show.pseudotime) & isTRUE(show.density)) {
      ht.anno[[i]] = HeatmapAnnotation(clust = anno_empty(height = unit(density.height, "inches"), border = F),
                                       pseudo = 1:length(pseudo.col[[i]]),
                                       col = list(pseudo = setNames(pseudo.col[[i]], 1:length(pseudo.col[[i]]))),
                                       annotation_label = c("", ifelse(i == length(lineages), "Pseudotime", "")),
                                       show_legend = F)
      ht.anno[[i]]@anno_list$clust@name = paste0("clust",i)
      ht = ht + Heatmap(mat.ht[[i]],
                        top_annotation = ht.anno[[i]],
                        show_heatmap_legend = ifelse(i == 1, TRUE, FALSE),
                        cluster_columns = F,
                        show_column_names = F,
                        show_row_names = ifelse(genes.names.size > 0, TRUE, FALSE),
                        row_km = ifelse(i == 1, genes.kmeans, 1),
                        row_km_repeats = ifelse(i == 1, kmeans.repeats, 1),
                        row_title_rot = 0,
                        row_title_gp = gpar(fontsize = genes.kmeans.numbers.size, col = ifelse(genes.kmeans.numbers.size > 0, "black", "white")),
                        width = unit(heatmap.width, "inches"),
                        height = unit(heatmap.height/8, "inches"),
                        row_names_gp = gpar(fontsize = genes.names.size, fontface = genes.names.style),
                        column_title = heatmap.named[i],
                        column_title_gp = gpar(fontsize = 20),
                        col = data.colors,
                        name = paste0("ht",i),
                        use_raster = raster,
                        heatmap_legend_param = list(title = "Smoothed Expression",
                                                    legend_direction = "horizontal",
                                                    title_position = "topcenter",
                                                    legend_width = unit(5, "cm")),
                        ...)
    }
    else if (isTRUE(show.pseudotime)) {
      ht.anno[[i]] = HeatmapAnnotation(pseudo = 1:length(pseudo.col[[i]]),
                                       col = list(pseudo = setNames(pseudo.col[[i]], 1:length(pseudo.col[[i]]))),
                                       annotation_label = ifelse(i == length(lineages), "Pseudotime", ""),
                                       show_legend = F)
      ht = ht + Heatmap(mat.ht[[i]],
                        top_annotation = ht.anno[[i]],
                        show_heatmap_legend = ifelse(i == 1, TRUE, FALSE),
                        cluster_rows = cluster.genes,
                        cluster_columns = F,
                        show_column_names = F,
                        show_row_names = ifelse(genes.names.size > 0, TRUE, FALSE),
                        row_km = ifelse(i == 1, genes.kmeans, 1),
                        row_km_repeats = ifelse(i == 1, kmeans.repeats, 1),
                        row_title_rot = 0,
                        row_title_gp = gpar(fontsize = genes.kmeans.numbers.size, col = ifelse(genes.kmeans.numbers.size > 0, "black", "white")),
                        width = unit(heatmap.width, "inches"),
                        height = unit(heatmap.height/8, "inches"),
                        row_names_gp = gpar(fontsize = genes.names.size, fontface = genes.names.style),
                        column_title = heatmap.named[i],
                        column_title_gp = gpar(fontsize = 20),
                        col = data.colors,
                        name = paste0("ht",i),
                        use_raster = raster,
                        heatmap_legend_param = list(title = "Smoothed Expression",
                                                    legend_direction = "horizontal",
                                                    title_position = "topcenter",
                                                    legend_width = unit(5, "cm")),
                        ...)
    }
    else if (isTRUE(show.density)) {
      ht.anno[[i]] = HeatmapAnnotation(clust = anno_empty(height = unit(density.height, "inches"), border = F),
                                       show_legend = F)
      ht.anno[[i]]@anno_list$clust@name = paste0("clust",i)
      ht = ht + Heatmap(mat.ht[[i]],
                        top_annotation = ht.anno[[i]],
                        show_heatmap_legend = ifelse(i == 1, TRUE, FALSE),
                        cluster_rows = cluster.genes,
                        cluster_columns = F,
                        show_column_names = F,
                        show_row_names = ifelse(genes.names.size > 0, TRUE, FALSE),
                        row_km = ifelse(i == 1, genes.kmeans, 1),
                        row_km_repeats = ifelse(i == 1, kmeans.repeats, 1),
                        row_title_rot = 0,
                        row_title_gp = gpar(fontsize = genes.kmeans.numbers.size, col = ifelse(genes.kmeans.numbers.size > 0, "black", "white")),
                        width = unit(heatmap.width, "inches"),
                        height = unit(heatmap.height/8, "inches"),
                        row_names_gp = gpar(fontsize = genes.names.size, fontface = genes.names.style),
                        column_title = heatmap.named[i],
                        column_title_gp = gpar(fontsize = 20),
                        col = data.colors,
                        name = paste0("ht",i),
                        use_raster = raster,
                        heatmap_legend_param = list(title = "Smoothed Expression",
                                                    legend_direction = "horizontal",
                                                    title_position = "topcenter",
                                                    legend_width = unit(5, "cm")),
                        ...)
    }
    else {
      ht = ht + Heatmap(mat.ht[[i]],
                        show_heatmap_legend = ifelse(i == 1, TRUE, FALSE),
                        cluster_rows = cluster.genes,
                        cluster_columns = F,
                        show_column_names = F,
                        show_row_names = ifelse(genes.names.size > 0, TRUE, FALSE),
                        row_km = ifelse(i == 1, genes.kmeans, 1),
                        row_km_repeats = ifelse(i == 1, kmeans.repeats, 1),
                        row_title_rot = 0,
                        row_title_gp = gpar(fontsize = genes.kmeans.numbers.size, col = ifelse(genes.kmeans.numbers.size > 0, "black", "white")),
                        width = unit(heatmap.width, "inches"),
                        height = unit(heatmap.height/8, "inches"),
                        row_names_gp = gpar(fontsize = genes.names.size, fontface = genes.names.style),
                        column_title = heatmap.named[i],
                        column_title_gp = gpar(fontsize = 20),
                        col = data.colors,
                        name = paste0("ht",i),
                        use_raster = raster,
                        heatmap_legend_param = list(title = "Smoothed Expression",
                                                    legend_direction = "horizontal",
                                                    title_position = "topcenter",
                                                    legend_width = unit(5, "cm")),
                        ...)
    }
  }

  if (isTRUE(return.grob)) {
    htt = grid.grabExpr({draw(ht,
                              heatmap_legend_side = "bottom",
                              align_heatmap_legend = "heatmap_center",
                              align_annotation_legend = "heatmap_center",
                              annotation_legend_list = density.legend,
                              legend_grouping = "original")

      for (i in 1:length(lineages)) {
        if (!is.null(ht.split[[i]])) {
          for (j in 1:length(ht.split[[i]])) {
            for (k in 1:genes.kmeans) {
              tryCatch(decorate_heatmap_body(paste0("ht",i), {
                grid.lines(c(ht.split[[i]][j], ht.split[[i]][j]), c(0, 0.99), gp = gpar(lty = 2, lwd = 2, col = "lightblue"))
              }, slice = k), error = function(e) {genes.kmeans <<- k-1})
            }
          }
        }
      }

      if (isTRUE(show.density)) {
        pp = list()
        for (i in 1:length(p)) {
          pp[[i]] = grid.grabExpr(suppressWarnings(print(p[[i]])))
          decorate_annotation(paste0("clust",i), {grid.draw(pp[[i]])})
        }
      }
    })

    plot.new()
    grid.draw(htt)
  }
  else {
    htt = draw(ht,
               heatmap_legend_side = "bottom",
               align_heatmap_legend = "heatmap_center",
               align_annotation_legend = "heatmap_center",
               annotation_legend_list = density.legend,
               legend_grouping = "original")

    for (i in 1:length(lineages)) {
      if (!is.null(ht.split[[i]])) {
        for (j in 1:length(ht.split[[i]])) {
          for (k in 1:genes.kmeans) {
            tryCatch(decorate_heatmap_body(paste0("ht",i), {
              grid.lines(c(ht.split[[i]][j], ht.split[[i]][j]), c(0, 0.99), gp = gpar(lty = 2, lwd = 2, col = "lightblue"))
            }, slice = k), error = function(e) {})
          }
        }
      }
    }

    if (isTRUE(show.density)) {
      pp = list()
      for (i in 1:length(p)) {
        pp[[i]] = grid.grabExpr(suppressWarnings(print(p[[i]])))
        decorate_annotation(paste0("clust",i), {grid.draw(pp[[i]])})
      }
    }
  }

  return(invisible(htt))
}
