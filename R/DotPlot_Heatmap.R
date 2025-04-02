#' @title Dot plot or heatmap of average gene expression in each identity
#'
#' @description This function generates a dot plot or a heatmap to visualize the average expression of features in each identity of the active.ident metadata of a \pkg{Seurat} object. Credits to \pkg{Seurat}'s dev team for the original \code{\link[Seurat]{DotPlot}} from which data processing of this function is derived from and to \href{https://divingintogeneticsandgenomics.com/post/clustered-dotplot-for-single-cell-rnaseq/}{Ming Tang} for the initial idea to use \pkg{ComplexHeatmap} to draw a dot plot and the \code{layer_fun} function that draws the dots. Various new parameters were added to offer more flexibility and customization.
#'
#' @param seurat_object A \pkg{Seurat} object.
#' @param assay Character. The name of an assay containing the \code{layer} with the expression matrix. If the \code{seurat_object} contains multiple 'RNA' assays, you may specify which one to use (for example, 'RNA2' if you have created a second 'RNA' assay you named 'RNA2'. See \href{https://satijalab.org/seurat/articles/seurat5_essential_commands.html#create-seurat-or-assay-objects}{Seurat v5 vignettes} for more information). You may also use another assay, such as 'SCT', to pull feature expression from.
#' @param layer Character. The name of a layer (formerly known as slot) which stores the expression matrix. It is recommended to use 'data'.
#' @param data.are.log Logical. If \code{TRUE}, and if \code{layer} = 'data', the function assumes data are log-transformed and cell expression values will be exponentiated (using \code{\link[base]{expm1}}) so that averaging is done in non-log space (as per \code{\link[Seurat]{DotPlot}} or \code{\link[Seurat]{AverageExpression}}'s default behavior), average expression values are then log-transformed back (using \code{\link[base]{log1p}}). If \code{layer} = 'scale.data' or 'counts', cell expression values are not exponentiated prior to averaging.
#' @param features Character. The names of one or several features to plot the average expression from.
#' @param split.by Character. The name of a metadata (for example, 'orig.ident', 'seurat_clusters', etc) to split the identities of the active.ident metadata by.
#' @param idents Character. The names of one or several identities in the active.ident metadata to select. If \code{NULL}, all identities are used.
#' @param split.idents Character. The names of one or several \code{split.by} identities to select. If \code{NULL}, all identities are used. Ignored if \code{split.by} = \code{NULL}.
#' @param scale Logical. If \code{TRUE}, average expression values will be scaled using \code{\link[base]{scale}} and default parameters. The resulting values will be Z-scores (mean subtracted values divided by standard deviation) and not positive average expression values anymore, which is why there will be positive and negative values displayed, depending on if the average expression in a particular identity is below or above the mean average expression from all identities (which is calculated independently for each feature). Caution should be exercised when interpreting results with low number of identities (typically below five), as small differences in average expression might lead to exacerbated differences when scaled.
#' @param rescale Logical. If \code{TRUE}, average expression values will be adjusted using \code{\link[scales]{rescale}} between the first numerical value of \code{rescale.range} (lowest expression) and the second numerical value (highest expression). This is different than \code{\link[base]{scale}} as this doesn't compare values to any mean and standard deviation and is therefore not a Z-score, it only refits each average expression value (independently for each feature) in order to visualize all \code{features} in the same dimension regardless of their differences in levels of expression. Caution should be exercised when interpreting results with low number of identities (typically below five), as small differences in average expression might lead to exacerbated differences when rescaled. Ignored if \code{scale} = \code{TRUE}.
#' @param rescale.range Numeric. The minimum and maximum values to resize the average expression values and internally passed to \code{\link[scales]{rescale}}. These values are arbitrary and will not change the visualization, only the values in the legend, you need to adjust \code{col.min} and \code{col.max} to influence the color scale. Ignored if \code{rescale} = \code{FALSE} or \code{scale} = \code{TRUE}.
#' @param rotate.axis Logical. If \code{TRUE}, flips the axis, so that \code{features} are displayed as rows and identities as columns.
#' @param dotplot Logical. If \code{TRUE}, the function will display a dot plot, with dots size proportional to the percentage of cells in each identity expressing a feature. If \code{FALSE}, the function will instead display a heatmap.
#' @param dots.type Character. Determines the dot size differences between 0 and 100% expression. Either 'square root' (lower difference) or 'radius' (higher difference). Ignored if \code{dotplot} = \code{FALSE}.
#' @param dots.size Numeric. The size of the dots. Decreasing this parameter helps when displaying a large number of \code{features}. Ignored if \code{dotplot} = \code{FALSE}.
#' @param show.noexpr.dots Logical. If \code{TRUE}, the function will display a small dot for \code{features} with 0% expression, instead of nothing. Ignored if \code{dotplot} = \code{FALSE}.
#' @param col.min Character or Numeric. The minimum value for the \code{breaks} internally passed to \code{\link[colorRamp2]{colorRamp2}}. If character, must be a quantile in the form 'qX' where X is a number between 0 and 100. A value of 'q5' or 'q10' is useful to reduce the effect of outlier values (i.e. a very low value that significantly alters the color scale range of all other values).
#' @param col.max Character or Numeric. The maximum value for the \code{breaks} internally passed to \code{\link[colorRamp2]{colorRamp2}}. If character, must be a quantile in the form 'qX' where X is a number between 0 and 100. A value of 'q95' or 'q90' is useful to reduce the effect of outlier values (i.e. a very high value that significantly alters the color scale range of all other values).
#' @param data.colors Character or Function. Either three color names, corresponding to the lowest, zero (or middle if \code{scale} = \code{FALSE}), and highest values in the expression matrix and internally passed to \code{\link[colorRamp2]{colorRamp2}}, or two color names, corresponding to the lowest and highest values, or the name of a palette and internally passed to \code{hcl_palette} in \code{\link[colorRamp2]{colorRamp2}} (such as 'Inferno', 'Berlin', 'Viridis' etc, check \code{\link[grDevices]{hcl.pals}} for all palettes available), or a custom \code{\link[colorRamp2]{colorRamp2}} function.
#' @param palette.reverse Logical. If \code{TRUE} and if \code{data.colors} is a palette (such as 'Viridis'), the function will reverse its colors.
#' @param na.color Character. The color name for missing values (\code{NA}).
#' @param background.color Character. The color name for the background behind the dots. Ignored if \code{dotplot} = \code{FALSE}.
#' @param idents.colors Character. The color names for each identity of the active.ident metadata or in \code{idents}. If \code{NULL}, uses \pkg{Seurat}'s default colors.
#' @param show.idents.names.colors Logical. If \code{TRUE}, the function will display the colors specified in \code{idents.colors} next to identity names.
#' @param show.idents.oppo.colors Logical. If \code{TRUE}, the function will display the colors specified in \code{idents.colors} on the opposite side of identity names.
#' @param split.colors Character. The color names for each \code{split.by} identity or in \code{split.idents}. If \code{NULL}, uses a custom set of colors from \code{\link[grDevices]{colors}}. Ignored if \code{split.by} = \code{NULL}.
#' @param show.split.names.colors Logical. If \code{TRUE}, the function will display the colors specified in \code{split.colors} next to identity names. Ignored if \code{split.by} = \code{NULL}.
#' @param show.split.oppo.colors Logical. If \code{TRUE}, the function will display the colors specified in \code{split.colors} on the opposite side of identity names. Ignored if \code{split.by} = \code{NULL}.
#' @param order.idents Character or Numeric. Either 'reverse', or the identities (as names or as numerical values corresponding to the indices) of the active.ident metadata or in \code{idents} to order the cells. If \code{cluster.idents} = \code{TRUE} or Function, only the legend names will be ordered.
#' @param order.split Character or Numeric. Either 'reverse', or the \code{split.by} identities (as names or as numerical values corresponding to the indices) or in \code{split.idents} to order the cells. If \code{cluster.idents} = \code{TRUE} or Function, only the legend names will be ordered. Ignored if \code{split.by} = \code{NULL}.
#' @param order.colors Logical. If \code{TRUE}, the \code{data.colors} and \code{split.colors} will automatically be ordered according to \code{order.idents} and \code{order.split}. Ignored if \code{order.idents} and \code{order.split} are \code{NULL}.
#' @param kmeans.repeats Numeric. The number of runs to get a consensus K-means clustering. Ignored if \code{idents.kmeans} and \code{features.kmeans} are equal to 1.
#' @param cluster.idents Logical or Function. If \code{TRUE}, the function will cluster the identities. You may also pass an \code{hclust} or \code{dendrogram} object which contains clustering.
#' @param idents.kmeans Numeric. The number of slices to use for identity K-means clustering.
#' @param idents.kmeans.numbers.size Numeric. The font size of the identity K-means slice numbers. Set to 0 to remove them.
#' @param cluster.features Logical or Function. If \code{TRUE}, the function will cluster the \code{features}. You may also pass an \code{hclust} or \code{dendrogram} object which contains clustering.
#' @param features.kmeans Numeric. The number of slices to use for feature K-means clustering.
#' @param features.kmeans.numbers.size Numeric. The font size of the feature K-means slice numbers. Set to 0 to remove them.
#' @param idents.gap Numeric. The gap between the identity slices. Ignored if \code{idents.kmeans} = 1.
#' @param features.gap Numeric. The gap between the feature slices. Ignored if \code{features.kmeans} = 1.
#' @param idents.names.size Numeric. The font size of the identity names. Set to 0 to remove them.
#' @param features.names.size Numeric. The font size of the feature names. Set to 0 to remove them.
#' @param features.names.style Character. The font face of the feature names. The \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7494048/}{Gene nomenclature} used by almost all scientific journals require that feature names are italicized, therefore the parameter is by default set to 'italic'. Use 'plain' to revert back to regular font face.
#' @param row.names.side Character. The side where the row names will be displayed, either 'left' or 'right'. The dendrogram will be displayed on the opposite side.
#' @param row.names.width Numeric. The width of the row names. Increase this parameter if your row names are truncated.
#' @param column.names.angle Numeric. The angle of rotation of the column names.
#' @param column.names.side Character. The side where the column names will be displayed, either 'top' or 'bottom'. The dendrogram will be displayed on the opposite side.
#' @param column.names.height Numeric. The height of the column names. Increase this parameter if your column names are truncated.
#' @param inner.border Logical. If \code{TRUE}, the function will display a black outline around each dot if \code{dotplot} = \code{TRUE}, or a black border around each cell of the heatmap if \code{dotplot} = \code{FALSE}.
#' @param outer.border Logical. If \code{TRUE}, the function will display an outer border around the plot or around each slice if \code{idents.kmeans} and/or \code{features.kmeans} are higher than 1.
#' @param data.legend.name Character. The name of the data legend.
#' @param data.legend.side Character. The side where the data legend will be displayed, either 'left', 'right', 'top' or 'bottom'.
#' @param data.legend.direction Character. The direction of the data legend, either 'horizontal' or 'vertical'.
#' @param data.legend.position Character. The centering of the data legend name, there are many options, default option from \code{\link[ComplexHeatmap]{Heatmap}} is 'topleft'.
#' @param data.legend.width Numeric. How long the data legend will be, only affects the data legend if \code{data.legend.direction} = 'horizontal'.
#' @param show.data.legend Logical. If \code{TRUE}, the function will display a legend for the average expression data.
#' @param idents.legend.name Character. The name of the active.ident metadata legend. Ignored if \code{show.idents.names.colors} and \code{show.idents.oppo.colors} are \code{FALSE}.
#' @param show.idents.legend Logical. If \code{TRUE}, the function will display a legend for the active.ident metadata identities or \code{idents}. Ignored if \code{show.idents.names.colors} and \code{show.idents.oppo.colors} are \code{FALSE}.
#' @param split.legend.name Character. The name of the \code{split.by} legend. Ignored if \code{split.by} = \code{NULL}. Ignored if \code{show.split.names.colors} and \code{show.split.oppo.colors} are \code{FALSE}.
#' @param show.split.legend Logical. If \code{TRUE}, the function will display a legend for \code{split.by} identities or \code{split.idents}. Ignored if \code{show.split.names.colors} and \code{show.split.oppo.colors} are \code{FALSE}.
#' @param legend.title.size Numeric. The font size of all legend titles.
#' @param legend.text.size Numeric. The font size of all legend texts.
#' @param legend.gap Numeric. The gap between the legends and the dot plot or heatmap. This parameter sets the value in the global options of \code{\link[ComplexHeatmap]{ht_opt}}, so it will affect all \code{\link[ComplexHeatmap]{Heatmap}} objects in the same R session. Use \pkg{ComplexHeatmap}::\code{\link[ComplexHeatmap]{ht_opt}}(RESET = \code{TRUE}) to restore default parameters.
#' @param output.data Logical. If \code{TRUE}, the function will return a \code{list} containing a \code{matrix} object of the average expression data, scaled or not, and another \code{matrix} object containing the percentage of cells expressing each feature, instead of displaying anything.
#' @param ... Additional arguments to be passed to \code{\link[ComplexHeatmap]{Heatmap}}, such as \code{show_parent_dend_line}, \code{clustering_method_rows}, etc, accepts any parameter that wasn't already internally passed to \code{\link[ComplexHeatmap]{Heatmap}} (for example, \code{outer.border} sets the \code{border} parameter of \code{\link[ComplexHeatmap]{Heatmap}}, so you will get an error if you try to pass the \code{border} parameter in \code{\link[RightOmicsTools]{DotPlot_Heatmap}}).
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap}} object, either as a dot plot, or a heatmap, or a \code{list} containing a \code{matrix} object of the average expression data, scaled or not, and another \code{matrix} object containing the percentage of cells expressing each feature.
#'
#' @examples
#' \dontshow{
#' suppressWarnings(suppressPackageStartupMessages(library(Seurat)))
#' }
#' # Prepare data
#' pbmc3k <- Right_Data("pbmc3k")
#' pbmc3k.markers <- c("CCR7", "CD14", "CD40LG",
#'                     "CD79A", "CD8A", "CDKN1C",
#'                     "GNLY", "CLEC10A", "PPBP")
#'
#' # Example 1: default parameters
#' DotPlot_Heatmap(pbmc3k,
#'                 features = pbmc3k.markers)
#' @import Seurat
#' @import SeuratObject
#' @import scales
#' @import grid
#' @import ComplexHeatmap
#' @import grDevices
#' @import colorRamp2
#' @importFrom stats setNames quantile
#' @export

DotPlot_Heatmap = function(seurat_object,
                           assay = "RNA",
                           layer = "data",
                           data.are.log = TRUE,
                           features,
                           split.by = NULL,
                           idents = NULL,
                           split.idents = NULL,
                           scale = TRUE,
                           rescale = FALSE,
                           rescale.range = c(0, 3),
                           rotate.axis = FALSE,
                           dotplot = TRUE,
                           dots.type = "square root",
                           dots.size = 4,
                           show.noexpr.dots = FALSE,
                           col.min = ifelse(isTRUE(scale), -2, 0),
                           col.max = ifelse(isTRUE(scale), 2, "q100"),
                           data.colors = if (isTRUE(scale)) c("#35A5FF","white","red") else "Viridis",
                           palette.reverse = FALSE,
                           na.color = "grey40",
                           background.color = "white",
                           idents.colors = NULL,
                           show.idents.names.colors = FALSE,
                           show.idents.oppo.colors = TRUE,
                           split.colors = NULL,
                           show.split.names.colors = FALSE,
                           show.split.oppo.colors = TRUE,
                           order.idents = NULL,
                           order.split = NULL,
                           order.colors = TRUE,
                           kmeans.repeats = 100,
                           cluster.idents = TRUE,
                           idents.kmeans = 1,
                           idents.kmeans.numbers.size = 11,
                           cluster.features = TRUE,
                           features.kmeans = 1,
                           features.kmeans.numbers.size = 11,
                           idents.gap = 1,
                           features.gap = 1,
                           idents.names.size = 9,
                           features.names.size = 9,
                           features.names.style = "italic",
                           row.names.side = "left",
                           row.names.width = 15,
                           column.names.angle = 45,
                           column.names.side = "bottom",
                           column.names.height = 15,
                           inner.border = TRUE,
                           outer.border = TRUE,
                           data.legend.name = ifelse(isTRUE(scale),"Z-Score","Average Expression"),
                           data.legend.side = "bottom",
                           data.legend.direction = "horizontal",
                           data.legend.position = "topcenter",
                           data.legend.width = 5,
                           show.data.legend = TRUE,
                           idents.legend.name = "Clusters",
                           show.idents.legend = TRUE,
                           split.legend.name = split.by,
                           show.split.legend = TRUE,
                           legend.title.size = 10,
                           legend.text.size = 10,
                           legend.gap = 10,
                           output.data = FALSE,
                           ...) {

  if (isFALSE(any(Assays(seurat_object) %in% assay))) {
    message("Assay '",assay,"' was not found in the Seurat object, using 'RNA' instead")
    assay = "RNA"
  }

  if (isFALSE(any(Layers(seurat_object[[assay]]) %in% layer))) {
    if (isTRUE(any(Layers(seurat_object[[assay]]) %in% "data"))) {
      message("Layer '",layer,"' was not found in the Seurat object's '",assay,"' assay, using 'data' instead")
      layer = "data"
    }
    else {
      message("Layer '",layer,"' was not found in the Seurat object's '",assay,"' assay, using 'counts' instead")
      layer = "counts"
    }
  }

  if (!is.character(split.by)) {
    ident = "ident"
  }
  else {
    ident = c("ident",split.by)
  }

  ident.1 = levels(Idents(seurat_object))

  if (is.character(idents)) {
    ident.1 = ident.1[ident.1 %in% idents]
    if (length(ident.1) == 0) {
      stop("None of the identities supplied to idents were found")
    }
    if (length(ident.1) < length(idents)) {
      message("The following identities supplied to idents were not found:\n", paste0(setdiff(idents, ident.1), collapse = ", "))
    }
  }
  if (is.character(order.idents) | is.numeric(order.idents)) {
    if (length(order.idents) == length(ident.1)) {
      if (is.character(order.idents)) {
        ident.1 = ident.1[order(match(ident.1, order.idents))]
      }
      else {
        ident.1 = ident.1[order.idents]
      }
    }
    else {
      if (length(order.idents) == 1 & any(order.idents == "reverse")) {
        ident.1 = rev(ident.1)
      }
      else {
        stop("order.idents needs to be 'reverse' or a character/numeric vector of same length as the number of identities")
      }
    }
  }
  if (length(ident.1) == 1 & is.null(split.by)) {
    stop("Only one identity in the active.ident identity, cannot compare features expression")
  }

  if (is.character(split.by)) {
    vars = c("ident",split.by,features)
  }
  else {
    vars = c("ident",features)
  }
  DefaultAssay(seurat_object) = assay
  data = suppressWarnings(FetchData(object = seurat_object, vars = vars, layer = layer))
  features.removed = setdiff(features, colnames(data))
  if (is.character(split.by)) {
    if (!split.by %in% colnames(data)) {
      message("The split.by identity was not found in the meta.data slot, data will not be split")
      split.by = NULL
      ident.3 = ident.1
    }
    else {
      Idents(seurat_object) = split.by
      ident.2 = levels(Idents(seurat_object))
      if (is.character(split.idents)) {
        ident.2 = ident.2[ident.2 %in% split.idents]
        if (length(ident.2) == 0) {
          stop("None of the identities supplied to split.idents were found")
        }
        if (length(ident.2) < length(split.idents)) {
          message("The following identities supplied to split.idents were not found:\n", paste0(setdiff(split.idents, ident.2), collapse = ", "))
        }
      }
        if (is.character(order.split) | is.numeric(order.split)) {
        if (length(order.split) == length(ident.2)) {
          if (is.character(order.split)) {
            ident.2 = ident.2[order(match(ident.2, order.split))]
          }
          else {
            ident.2 = ident.2[order.split]
          }
        }
        else {
          if (length(order.split) == 1 & any(order.split == "reverse")) {
            ident.2 = rev(ident.2)
          }
          else {
            stop("order.split needs to be 'reverse' or a character/numeric vector of same length as the number of identities")
          }
        }
      }
      ident.3 = NULL
      for (i in ident.1) {
        for (j in ident.2) {
          tmp = paste(i, "in", j)
          ident.3 = c(ident.3,tmp)
        }
      }
    }
    if (length(ident.3) == 1) {
      stop("Only one identity in the active.ident and split.by identities, cannot compare features expression")
    }
  }
  else {
    ident.3 = ident.1
  }
  if (is.character(split.by)) {
    data$ident = paste(data$ident, "in", data[ , 2])
    data[ , 2] = NULL
  }
  if (ncol(data) == 1) {
    stop("None of the features were found")
  }
  mat = dot = list()
  for (i in 2:ncol(data)) {
    mat[[i-1]] = data[ , c(1,i), drop = FALSE]
    avg = pct = list()
    for (j in ident.3) {
      avg[[j]] = mat[[i-1]][mat[[i-1]]$ident == j, 2, drop = FALSE]
      pct[[j]] = as.data.frame(length(avg[[j]][avg[[j]] > 0])/length(avg[[j]][ , 1]))
      if (isTRUE(data.are.log) & layer == "data") {
        avg[[j]] = as.data.frame(mean(expm1(avg[[j]][ , 1])))
      }
      else {
        avg[[j]] = as.data.frame(mean(avg[[j]][ , 1]))
      }
      colnames(avg[[j]]) = colnames(data)[i]
      rownames(avg[[j]]) = j
      colnames(pct[[j]]) = colnames(data)[i]
      rownames(pct[[j]]) = j
    }
    mat[[i-1]] = do.call(rbind, avg)
    dot[[i-1]] = do.call(rbind, pct)
  }
  mat = as.matrix(as.data.frame(mat))[ident.3, , drop = FALSE]
  colnames(mat) = colnames(data)[2:ncol(data)]
  dot = as.matrix(as.data.frame(dot))[ident.3, , drop = FALSE]
  colnames(dot) = colnames(data)[2:ncol(data)]
  if (isTRUE(data.are.log) & layer == "data") {
    mat = log1p(mat)
  }
  dot = dot[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  mat = mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  if (nrow(mat) < 2) {
    stop("Less than two identities left, cannot compare features expression")
  }
  idents.removed = setdiff(ident.3, rownames(mat))
  dot = dot[ , colSums(mat) > 0, drop = FALSE]
  mat = mat[ , colSums(mat) > 0, drop = FALSE]
  if (ncol(mat) == 0) {
    stop("None of the features were expressed in any cells")
  }
  features.removed = c(features.removed,setdiff(colnames(data), colnames(mat)))
  features.removed = features.removed[!features.removed == "ident"]
  if (length(features.removed) > 0) {
    message("The following features were removed as they were not found or were not expressed in any cells:\n",paste0(features.removed, collapse = ", "))
  }
  if (isTRUE(scale)) {
    if (nrow(mat) < 5 & isFALSE(output.data)) {
      message("Less than 5 identities will be displayed, scaling may produce misleading visualization.\nUsing Cell_Heatmap() might be more appropriate.")
    }
    mat = scale(mat)
  }
  if (isTRUE(rescale) & isFALSE(scale)) {
    if (nrow(mat) < 5 & isFALSE(output.data)) {
      message("Less than 5 identities will be displayed, rescaling may produce misleading visualization.\nUsing Cell_Heatmap() might be more appropriate.")
    }
    mat = apply(mat, 2, function(x) rescale(x, to = rescale.range))
  }

  if (is.character(col.min)) {
    if (isTRUE(grepl("^q[0-9]{1,3}$", col.min, perl = TRUE))) {
      q1 = as.numeric(sub("q","", col.min)) / 100
      q1 = quantile(mat, probs = q1, na.rm = TRUE)
      col.min = q1
    }
  }
  if (is.character(col.max)) {
    if (isTRUE(grepl("^q[0-9]{1,3}$", col.max, perl = TRUE))) {
      q2 = as.numeric(sub("q","", col.max)) / 100
      q2 = quantile(mat, probs = q2, na.rm = TRUE)
      col.max = q2
    }
  }

  if (isTRUE(rotate.axis)) {
    mat = t(mat)
    dot = t(dot)
  }

  if (isTRUE(output.data)) {
    return(list(Expression = round(mat, 4), Percent = round(dot, 4)))
  }

  if (isTRUE(scale)) {
    col.mid = 0
  }
  else {
    col.mid = (col.min+col.max)/2
  }
  if (!is.function(data.colors)) {
    if (length(data.colors) > 3) {
      stop("data.colors must have 2 or 3 colors or must be a palette")
    }
    if (length(data.colors) == 3) {
      data.colors = colorRamp2(breaks = c(col.min,col.mid,col.max),
                               colors = c(data.colors[1],data.colors[2],data.colors[3]))
    }
    else if (length(data.colors) == 2) {
      data.colors = colorRamp2(breaks = c(col.min,col.max),
                               colors = c(data.colors[1],data.colors[2]))
    } else {
      data.colors = colorRamp2(breaks = c(col.min,col.max),
                               hcl_palette = data.colors, reverse = palette.reverse)
    }
  }

  if (!is.character(split.colors) & is.character(split.by)) {
    current.seed = .Random.seed
    set.seed(1234)
    CustomColors = colors()[grep('gr(a|e)y|white|black', colors(), invert = T)]
    split.colors = sample(CustomColors, length(ident.2))
    .Random.seed = current.seed
  }

  if (!is.character(idents.colors)) {
    idents.colors = hue_pal()(n = length(ident.1))
  }
  if (isTRUE(order.colors)) {
    if (!is.null(names(idents.colors))) {
      idents.colors = idents.colors[ident.1]
    }
    if (is.character(order.idents)) {
      if (length(order.idents) > 1) {
        names(idents.colors) = ident.1
        idents.colors = idents.colors[order.idents]
      }
      else {
        if (order.idents == "reverse") {
          idents.colors = rev(idents.colors)
        }
      }
    }
    if (is.character(split.by)) {
      if (!is.null(names(split.colors))) {
        split.colors = split.colors[ident.2]
      }
      if (is.character(order.split)) {
        if (length(order.split) > 1) {
          names(split.colors) = ident.2
          split.colors = split.colors[order.split]
        }
        else {
          if (order.split == "reverse") {
            split.colors = rev(split.colors)
          }
        }
      }
    }
  }
  idents.colors2 = idents.colors
  if (is.character(split.by)) {
    split.colors2 = split.colors
    dup.colors = dup.colors2 = list()
    for (i in 1:length(idents.colors)) {
      dup.colors[[i]] = rep(idents.colors[i], length(ident.3)/length(ident.1))
      dup.colors2[[i]] = split.colors
    }
    idents.colors = unlist(dup.colors)
    split.colors = unlist(dup.colors2)
    names(idents.colors) = names(split.colors) = ident.3
    idents.colors = idents.colors[!names(idents.colors) %in% idents.removed]
    split.colors = split.colors[!names(split.colors) %in% idents.removed]
  }
  else {
    names(idents.colors) = ident.3
    idents.colors = idents.colors[!names(idents.colors) %in% idents.removed]
  }

  idents.legend = split.legend = NULL
  if (isTRUE(show.idents.legend) & (isTRUE(show.idents.names.colors) | isTRUE(show.idents.oppo.colors))) {
    idents.legend = Legend(at = ident.1,
                           legend_gp = gpar(fill = idents.colors2),
                           title = idents.legend.name,
                           gap = unit(0.5, "cm"),
                           border = TRUE,
                           title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
                           labels_gp = gpar(fontsize = legend.text.size))
  }
  if (is.character(split.by) & isTRUE(show.split.legend) & (isTRUE(show.split.names.colors) | isTRUE(show.split.oppo.colors))) {
    split.legend = Legend(at = ident.2,
                          legend_gp = gpar(fill = split.colors2),
                          title = split.legend.name,
                          gap = unit(0.5, "cm"),
                          border = TRUE,
                          title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
                          labels_gp = gpar(fontsize = legend.text.size))
  }
  anno.legend = list()
  if (!is.null(idents.legend)) {
    anno.legend = c(anno.legend, list(idents.legend))
  }
  if (!is.null(split.legend)) {
    anno.legend = c(anno.legend, list(split.legend))
  }
  if (is.character(split.by)) {
    ident.1 = names(idents.colors)
    dup.ident = names(split.colors)
  }
  if (isFALSE(rotate.axis)) {
    annotation.side = "row"
  }
  else {
    annotation.side = "column"
  }
  idents.labels3 = idents.labels4 = NULL
  row.dend.side = row.title.side = "right"
  column.dend.side = column.title.side = "top"
  idents.kmeans.color = ifelse(idents.kmeans.numbers.size > 0, "black", "white")
  features.kmeans.color = ifelse(features.kmeans.numbers.size > 0, "black", "white")
  if (is.character(split.by)) {
    if (isTRUE(show.split.names.colors) & isTRUE(show.idents.names.colors)) {
      if ((isFALSE(rotate.axis) & row.names.side == "left") | (isTRUE(rotate.axis) & column.names.side == "top")) {
        idents.labels = HeatmapAnnotation(
          Identities = ident.1,
          Split = dup.ident,
          col = list(Identities = idents.colors,
                     Split = split.colors),
          na_col = na.color,
          which = annotation.side,
          show_annotation_name = FALSE,
          show_legend = FALSE)
      }
      else {
        idents.labels = HeatmapAnnotation(
          Split = dup.ident,
          Identities = ident.1,
          col = list(Split = split.colors,
                     Identities = idents.colors),
          na_col = na.color,
          which = annotation.side,
          show_annotation_name = FALSE,
          show_legend = FALSE)
      }
    }
    else if (isTRUE(show.idents.names.colors)) {
      idents.labels = HeatmapAnnotation(
        Identities = ident.1,
        col = list(Identities = idents.colors),
        na_col = na.color,
        which = annotation.side,
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else if (isTRUE(show.split.names.colors)) {
      idents.labels = HeatmapAnnotation(
        Split = dup.ident,
        col = list(Split = split.colors),
        na_col = na.color,
        which = annotation.side,
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else {
      idents.labels = NULL
    }
    if(isTRUE(show.split.oppo.colors) & isTRUE(show.idents.oppo.colors)) {
        if ((isFALSE(rotate.axis) & row.names.side == "left") | (isTRUE(rotate.axis) & column.names.side == "top")) {
          idents.labels2 = HeatmapAnnotation(
            Split = dup.ident,
            Identities = ident.1,
            col = list(Split = split.colors,
                       Identities = idents.colors),
            na_col = na.color,
            which = annotation.side,
            show_annotation_name = FALSE,
            show_legend = FALSE)
        }
        else {
          idents.labels2 = HeatmapAnnotation(
            Identities = ident.1,
            Split = dup.ident,
            col = list(Identities = idents.colors,
                       Split = split.colors),
            na_col = na.color,
            which = annotation.side,
            show_annotation_name = FALSE,
            show_legend = FALSE)
        }
      }
    else if (isTRUE(show.idents.oppo.colors)) {
      idents.labels2 = HeatmapAnnotation(
        Identities = ident.1,
        col = list(Identities = idents.colors),
        na_col = na.color,
        which = annotation.side,
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else if (isTRUE(show.split.oppo.colors)) {
      idents.labels2 = HeatmapAnnotation(
        Split = dup.ident,
        col = list(Split = split.colors),
        na_col = na.color,
        which = annotation.side,
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else {
      idents.labels2 = NULL
    }
  }
  else {
    if (isTRUE(show.idents.names.colors)) {
      idents.labels = HeatmapAnnotation(
        Identities = ident.1,
        col = list(Identities = idents.colors),
        na_col = na.color,
        which = annotation.side,
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else {
      idents.labels = NULL
    }
    if (isTRUE(show.idents.oppo.colors)) {
    idents.labels2 = HeatmapAnnotation(
      Identities = ident.1,
      col = list(Identities = idents.colors),
      na_col = na.color,
      which = annotation.side,
      show_annotation_name = FALSE,
      show_legend = FALSE)
    }
    else {
      idents.labels2 = NULL
    }
  }
  if (row.names.side == "right") {
    if (isFALSE(rotate.axis)) {
    idents.labelstemp = idents.labels
    idents.labels = idents.labels2
    idents.labels2 = idents.labelstemp
    }
    row.dend.side = row.title.side = "left"
  }
  idents.names.gp = gpar(fontsize = idents.names.size)
  features.names.gp = gpar(fontsize = features.names.size, fontface = features.names.style)
  idents.title.gp = gpar(fontsize = idents.kmeans.numbers.size, col = idents.kmeans.color)
  features.title.gp = gpar(fontsize = features.kmeans.numbers.size, col = features.kmeans.color)
  if (isTRUE(rotate.axis)) {
    idents.labels3 = idents.labels
    idents.labels = NULL
    idents.labels4 = idents.labels2
    idents.labels2 = NULL
    idents.names.gptmp = idents.names.gp
    idents.names.gp = features.names.gp
    features.names.gp = idents.names.gptmp
    idents.names.sizetmp = idents.names.size
    idents.names.size = features.names.size
    features.names.size = idents.names.sizetmp
    cluster.identstmp = cluster.idents
    cluster.idents = cluster.features
    cluster.features = cluster.identstmp
    idents.kmeanstmp = idents.kmeans
    idents.kmeans = features.kmeans
    features.kmeans = idents.kmeanstmp
    idents.title.gptmp = idents.title.gp
    idents.title.gp = features.title.gp
    features.title.gp = idents.title.gptmp
    idents.gaptmp = idents.gap
    idents.gap = features.gap
    features.gap = idents.gaptmp
  }
  if (column.names.side == "top") {
    if (isTRUE(rotate.axis)) {
    idents.labelstemp = idents.labels3
    idents.labels3 = idents.labels4
    idents.labels4 = idents.labelstemp
    }
    column.dend.side = column.title.side = "bottom"
  }

  if (isFALSE(dotplot)) {
    if (sum(is.na(mat)) > 0) {
      na.label = list(Legend(at = "NA",
                        legend_gp = gpar(fill = na.color),
                        title = NULL,
                        gap = unit(0.5, "cm"),
                        border = TRUE))
    }
    else{
      na.label = list()
    }
    anno.legend = c(na.label,anno.legend)
  }

  ht_opt$ANNOTATION_LEGEND_PADDING = unit(legend.gap, "mm")
  if (data.legend.side == "right" | data.legend.side == "left" | column.names.side == "top") {
    ht_opt$HEATMAP_LEGEND_PADDING = unit(legend.gap, "mm")
  }
  else {
    ht_opt$HEATMAP_LEGEND_PADDING = unit(2, "mm")
  }
  if (isTRUE(dotplot)) {
    zero.dot = 0
    if (isTRUE(show.noexpr.dots)) {
      zero.dot = min(dot[dot > 0], na.rm = TRUE)
      if (zero.dot > 0.001) {
        zero.dot = 0.001
      }
      dot[dot==0] = zero.dot
    }
    if (dots.type == "radius") {
    dots.creation = function(j, i, x, y, w, h, fill){
      grid.rect(x = x, y = y, width = w, height = h,
                gp = gpar(col = NA, fill = background.color))
      ifelse(isTRUE(inner.border), grid.circle(x=x,y=y,r= pindex(dot, i, j) * unit(dots.size, "mm"),
                                               gp = gpar(fill = data.colors(pindex(mat, i, j)), col = "black", lty = 1)),
             grid.circle(x=x,y=y,r= pindex(dot, i, j) * unit(dots.size, "mm"),
                         gp = gpar(fill = data.colors(pindex(mat, i, j)), col = NA)))}
    dots.legend = Legend(labels = c(0,10,25,50,75,100), title = "Percent\nof cells\nexpressing",
                         row_gap = unit(legend.text.size/15, "mm"),
                         title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
                         labels_gp = gpar(fontsize = legend.text.size),
              graphics = list(
                function(x, y, w, h) grid.circle(x = x, y = y, r = zero.dot * unit(legend.text.size/6, "mm"),
                                                 gp = gpar(fill = "black")),
                function(x, y, w, h) grid.circle(x = x, y = y, r = 0.1 * unit(legend.text.size/6, "mm"),
                                                 gp = gpar(fill = "black")),
                function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(legend.text.size/6, "mm"),
                                                 gp = gpar(fill = "black")),
                function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(legend.text.size/6, "mm"),
                                                 gp = gpar(fill = "black")),
                function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(legend.text.size/6, "mm"),
                                                 gp = gpar(fill = "black")),
                function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(legend.text.size/6, "mm"),
                                                 gp = gpar(fill = "black"))))
    }
    if (dots.type == "square root") {
      dots.creation = function(j, i, x, y, w, h, fill){
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = NA, fill = background.color))
        ifelse(isTRUE(inner.border), grid.circle(x=x,y=y,r= sqrt(pindex(dot, i, j)) * unit(dots.size, "mm"),
                    gp = gpar(fill = data.colors(pindex(mat, i, j)), col = "black", lty = 1)),
               grid.circle(x=x,y=y,r= sqrt(pindex(dot, i, j)) * unit(dots.size, "mm"),
                           gp = gpar(fill = data.colors(pindex(mat, i, j)), col = NA)))}
      dots.legend = Legend(labels = c(0,10,25,50,75,100), title = "Percent\nof cells\nexpressing",
                           row_gap = unit(legend.text.size/15, "mm"),
                           title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
                           labels_gp = gpar(fontsize = legend.text.size),
                           graphics = list(
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(zero.dot) * unit(legend.text.size/6, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.1) * unit(legend.text.size/6, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(legend.text.size/6, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(legend.text.size/6, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(legend.text.size/6, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(1) * unit(legend.text.size/6, "mm"),
                                                              gp = gpar(fill = "black"))))
    }
    anno.legend = c(anno.legend,list(dots.legend))
    if (data.legend.side == "right") {
      show.data.legend = FALSE
      if (isTRUE(scale)) {
        data.legend.breaks = c(col.min, col.min/2, 0, col.max/2, col.max)
      }
      else {
        data.legend.breaks = c(col.min, col.min, col.max/2, col.max)
      }
      anno.legend = c(list(Legend(col_fun = data.colors, title = gsub("\\ ", "\n", data.legend.name),
                                  legend_height = unit(data.legend.width, "cm"),
                                  title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
                                  labels_gp = gpar(fontsize = legend.text.size),
                                  grid_width = unit(0.5, "cm"),
                                  title_gap = unit(legend.title.size/4, "mm"))), anno.legend)
    }

    ht = Heatmap(
      mat,
      col = data.colors,
      na_col = na.color,
      name = data.legend.name,
      rect_gp = gpar(type = "none"),
      border = outer.border,
      layer_fun = dots.creation,
      left_annotation = idents.labels,
      right_annotation = idents.labels2,
      bottom_annotation = idents.labels3,
      top_annotation = idents.labels4,
      cluster_rows = cluster.idents,
      cluster_columns = cluster.features,
      row_title_rot = 0,
      row_title_side = row.title.side,
      row_title_gp = idents.title.gp,
      row_km = idents.kmeans,
      row_km_repeats = kmeans.repeats,
      column_km = features.kmeans,
      column_km_repeats = kmeans.repeats,
      show_column_names = ifelse(features.names.size > 0, TRUE, FALSE),
      show_row_names = ifelse(idents.names.size > 0, TRUE, FALSE),
      column_names_gp = features.names.gp,
      row_names_gp = idents.names.gp,
      column_names_rot = column.names.angle,
      column_title_side = column.title.side,
      column_title_gp = features.title.gp,
      row_names_side = row.names.side,
      row_names_max_width = unit(row.names.width, "cm"),
      column_names_side = column.names.side,
      column_names_max_height = unit(column.names.height, "cm"),
      row_dend_side = row.dend.side,
      column_dend_side = column.dend.side,
      show_heatmap_legend = show.data.legend,
      heatmap_legend_param = list(
        legend_direction = data.legend.direction,
        title_position = data.legend.position,
        legend_width = unit(data.legend.width, "cm"),
        title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
        labels_gp = gpar(fontsize = legend.text.size)),
      ...)
  }

  else {
    if (isTRUE(inner.border)) {
      inner.border = gpar(col = "black")
    }
    else {
      inner.border = gpar(col = NA)
    }
    ht = Heatmap(
      mat,
      col = data.colors,
      na_col = na.color,
      name = data.legend.name,
      rect_gp = inner.border,
      border = outer.border,
      left_annotation = idents.labels,
      right_annotation = idents.labels2,
      bottom_annotation = idents.labels3,
      top_annotation = idents.labels4,
      cluster_rows = cluster.idents,
      cluster_columns = cluster.features,
      row_title_rot = 0,
      row_title_side = row.title.side,
      row_title_gp = idents.title.gp,
      row_km = idents.kmeans,
      row_km_repeats = kmeans.repeats,
      column_km = features.kmeans,
      column_km_repeats = kmeans.repeats,
      column_names_gp = features.names.gp,
      row_names_gp = idents.names.gp,
      column_names_rot = column.names.angle,
      column_title_side = column.title.side,
      column_title_gp = features.title.gp,
      row_names_side = row.names.side,
      row_names_max_width = unit(row.names.width, "cm"),
      column_names_side = column.names.side,
      column_names_max_height = unit(column.names.height, "cm"),
      row_dend_side = row.dend.side,
      column_dend_side = column.dend.side,
      row_gap = unit(idents.gap, "mm"),
      column_gap = unit(features.gap, "mm"),
      show_heatmap_legend = show.data.legend,
      heatmap_legend_param = list(
        legend_direction = data.legend.direction,
        title_position = data.legend.position,
        legend_width = unit(data.legend.width, "cm"),
        title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
        labels_gp = gpar(fontsize = legend.text.size)),
      ...)
  }

  plot.new()
  return(draw(ht,
              heatmap_legend_side = data.legend.side,
              align_heatmap_legend = "heatmap_center",
              align_annotation_legend = "heatmap_center",
              annotation_legend_list = anno.legend,
              legend_grouping = "original"))
}
