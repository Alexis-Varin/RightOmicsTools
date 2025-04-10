#' @title Heatmap of gene expression in each cell
#'
#' @description This function generates a heatmap to visualize the expression of features in each cell of a \pkg{Seurat} object. Credits to \href{https://divingintogeneticsandgenomics.com/post/enhancement-of-scrnaseq-heatmap-using-complexheatmap/}{Ming Tang} for the initial idea to replicate \code{\link[Seurat]{DoHeatmap}} using \pkg{ComplexHeatmap}. Various new parameters were added to offer more flexibility and customization.
#'
#' @param seurat_object A \pkg{Seurat} object.
#' @param assay Character. The name of an assay containing the \code{layer} with the expression matrix. If the \code{seurat_object} contains multiple 'RNA' assays, you may specify which one to use (for example, 'RNA2' if you have created a second 'RNA' assay you named 'RNA2'. See \href{https://satijalab.org/seurat/articles/seurat5_essential_commands.html#create-seurat-or-assay-objects}{Seurat v5 vignettes} for more information). You may also use another assay, such as 'SCT', to pull feature expression from.
#' @param layer Character. The name of a layer (formerly known as slot) which stores the expression matrix. It is recommended to use 'data'.
#' @param features Character. The names of one or several features to plot the cell expression from.
#' @param split.by Character. The name of a metadata (for example, 'orig.ident', 'seurat_clusters', etc) to split the identities of the active.ident metadata by.
#' @param idents Character. The names of one or several identities in the active.ident metadata to select. If \code{NULL}, all identities are used.
#' @param split.idents Character. The names of one or several \code{split.by} identities to select. If \code{NULL}, all identities are used. Ignored if \code{split.by} = \code{NULL}.
#' @param scale Logical. If \code{TRUE}, cell expression values will be scaled using \code{\link[base]{scale}} and default parameters. The resulting values will be Z-scores (mean subtracted values divided by standard deviation) and not positive cell expression values anymore, which is why there will be positive and negative values displayed, depending on if the expression in a particular cell is below or above the mean expression from all cells (which is calculated independently for each feature).
#' @param rescale Logical. If \code{TRUE}, cell expression values will be adjusted using \code{\link[scales]{rescale}} between the first numerical value of \code{rescale.range} (lowest expression) and the second numerical value (highest expression). This is different than \code{\link[base]{scale}} as this doesn't compare values to any mean and standard deviation and is therefore not a Z-score, it only refits each cell expression value (independently for each feature) in order to visualize all \code{features} in the same dimension regardless of their differences in levels of expression. Ignored if \code{scale} = \code{TRUE}.
#' @param rescale.range Numeric. The minimum and maximum values to resize the cell expression values and internally passed to \code{\link[scales]{rescale}}. These values are arbitrary and will not change the visualization, only the values in the legend, you need to adjust \code{col.min} and \code{col.max} to influence the color scale. Ignored if \code{rescale} = \code{FALSE} or \code{scale} = \code{TRUE}.
#' @param rotate.axis Logical. If \code{TRUE}, flips the axis, so that \code{features} are displayed as columns and identities as rows.
#' @param col.min Character or Numeric. The minimum value for the \code{breaks} internally passed to \code{\link[colorRamp2]{colorRamp2}}. If character, must be a quantile in the form 'qX' where X is a number between 0 and 100. A value of 'q5' or 'q10' is useful to reduce the effect of outlier values (i.e. a very low value that significantly alters the color scale range of all other values).
#' @param col.max Character or Numeric. The maximum value for the \code{breaks} internally passed to \code{\link[colorRamp2]{colorRamp2}}. If character, must be a quantile in the form 'qX' where X is a number between 0 and 100. A value of 'q95' or 'q90' is useful to reduce the effect of outlier values (i.e. a very high value that significantly alters the color scale range of all other values).
#' @param data.colors Character or Function. Either three color names, corresponding to the lowest, zero (or middle if \code{scale} = \code{FALSE}), and highest values in the expression matrix and internally passed to \code{\link[colorRamp2]{colorRamp2}}, or two color names, corresponding to the lowest and highest values, or the name of a palette and internally passed to \code{hcl_palette} in \code{\link[colorRamp2]{colorRamp2}} (such as 'Inferno', 'Berlin', 'Viridis' etc, check \code{\link[grDevices]{hcl.pals}} for all palettes available), or a custom \code{\link[colorRamp2]{colorRamp2}} function.
#' @param palette.reverse Logical. If \code{TRUE} and if \code{data.colors} is a palette (such as 'Viridis'), the function will reverse its colors.
#' @param na.color Character. The color name for missing values (\code{NA}).
#' @param idents.colors Character. The color names for each identity of the active.ident metadata or in \code{idents}. If \code{NULL}, uses \pkg{Seurat}'s default colors.
#' @param show.idents.names.colors Logical. If \code{TRUE}, the function will display the colors specified in \code{idents.colors} next to identity names.
#' @param show.idents.oppo.colors Logical. If \code{TRUE}, the function will display the colors specified in \code{idents.colors} on the opposite side of identity names.
#' @param split.colors Character. The color names for each \code{split.by} identity or in \code{split.idents}. If \code{NULL}, uses a custom set of colors from \code{\link[grDevices]{colors}}. Ignored if \code{split.by} = \code{NULL}.
#' @param show.split.names.colors Logical. If \code{TRUE}, the function will display the colors specified in \code{split.colors} next to identity names. Ignored if \code{split.by} = \code{NULL}.
#' @param show.split.oppo.colors Logical. If \code{TRUE}, the function will display the colors specified in \code{split.colors} on the opposite side of identity names. Ignored if \code{split.by} = \code{NULL}.
#' @param order.idents Character or Numeric. Either 'reverse', or the identities (as names or as numerical values corresponding to the indices) of the active.ident metadata or in \code{idents} to order the cells.
#' @param order.split Character or Numeric. Either 'reverse', or the \code{split.by} identities (as names or as numerical values corresponding to the indices) or in \code{split.idents} to order the cells. Ignored if \code{split.by} = \code{NULL}.
#' @param order.colors Logical. If \code{TRUE}, the \code{data.colors} and \code{split.colors} will automatically be ordered according to \code{order.idents} and \code{order.split}. Ignored if \code{order.idents} and \code{order.split} are \code{NULL}.
#' @param kmeans.repeats Numeric. The number of runs to get a consensus K-means clustering. Ignored if \code{features.kmeans} = 1.
#' @param shuffle.cells Logical. If \code{TRUE}, the function will randomize the distribution of cells in each identity. Useful to smooth expression, which limits visible batch effect (cells in a merged \pkg{Seurat} object are typically ordered based on the levels of the 'orig.ident' metadata, this might lead to unwanted patterns of expression in another metadata). Note that no values are modified, it only changes the order of cells in each identity. Ignored if \code{cluster.cells} = \code{TRUE}.
#' @param cluster.cells Logical. If \code{TRUE}, the function will cluster the cells within each identity. Will have the opposite effect of \code{shuffle.cells}, as it will order cells based on their expression similarity and will therefore increase batch effect. Useful to visualize if, within an identity, a subset of cells have high expression while the rest of the cells have low expression, or vice versa. Just like \code{shuffle.cells}, no values are modified, it only changes the order of cells in each identity.
#' @param cluster.features Logical or Function. If \code{TRUE}, the function will cluster the \code{features}. You may also pass an \code{hclust} or \code{dendrogram} object which contains clustering.
#' @param features.kmeans Numeric. The number of slices to use for feature K-means clustering.
#' @param features.kmeans.numbers.size Numeric. The font size of the feature K-means slice numbers. Set to 0 to remove them.
#' @param idents.gap Numeric. The gap between the identity slices.
#' @param features.gap Numeric. The gap between the feature slices. Ignored if \code{features.kmeans} = 1.
#' @param idents.names.size Numeric. The font size of the identity names. Set to 0 to remove them.
#' @param features.names.size Numeric. The font size of the feature names. Set to 0 to remove them.
#' @param features.names.style Character. The font face of the feature names. The \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7494048/}{Gene nomenclature} used by almost all scientific journals require that feature names are italicized, therefore the parameter is by default set to 'italic'. Use 'plain' to revert back to regular font face.
#' @param row.names.side Character. The side where the row names will be displayed, either 'left' or 'right'. The dendrogram will be displayed on the opposite side.
#' @param row.names.width Numeric. The width of the row names. Increase this parameter if your row names are truncated.
#' @param column.names.angle Numeric. The angle of rotation of the column names.
#' @param column.names.side Character. The side where the column names will be displayed, either 'top' or 'bottom'. The dendrogram will be displayed on the opposite side.
#' @param column.names.height Numeric. The height of the column names. Increase this parameter if your column names are truncated.
#' @param outer.border Logical. If \code{TRUE}, the function will display an outer border around the heatmap or around each slice if \code{features.kmeans} > 1.
#' @param data.legend.name Character. The name of the data legend.
#' @param data.legend.side Character. The side where the data legend will be displayed, either 'left', 'right', 'top' or 'bottom'.
#' @param data.legend.direction Character. The direction of the data legend, either 'horizontal' or 'vertical'.
#' @param data.legend.position Character. The centering of the data legend name, there are many options, default option from \code{\link[ComplexHeatmap]{Heatmap}} is 'topleft'.
#' @param data.legend.width Numeric. How long the data legend will be, only affects the data legend if \code{data.legend.direction} = 'horizontal'.
#' @param idents.legend.name Character. The name of the active.ident metadata legend. Ignored if \code{show.idents.names.colors} and \code{show.idents.oppo.colors} are \code{FALSE}.
#' @param show.idents.legend Logical. If \code{TRUE}, the function will display a legend for the active.ident metadata identities or \code{idents}. Ignored if \code{show.idents.names.colors} and \code{show.idents.oppo.colors} are \code{FALSE}.
#' @param split.legend.name Character. The name of the \code{split.by} legend. Ignored if \code{split.by} = \code{NULL}. Ignored if \code{show.split.names.colors} and \code{show.split.oppo.colors} are \code{FALSE}.
#' @param show.split.legend Logical. If \code{TRUE}, the function will display a legend for the \code{split.by} identities or \code{split.idents}. Ignored if \code{show.split.names.colors} and \code{show.split.oppo.colors} are \code{FALSE}.
#' @param legend.title.size Numeric. The font size of all legend titles.
#' @param legend.text.size Numeric. The font size of all legend texts.
#' @param legend.gap Numeric. The gap between the legends and the heatmap. This parameter sets the value in the global options of \code{\link[ComplexHeatmap]{ht_opt}}, so it will affect all \code{\link[ComplexHeatmap]{Heatmap}} objects in the same R session. Use \pkg{ComplexHeatmap}::\code{\link[ComplexHeatmap]{ht_opt}}(RESET = \code{TRUE}) to restore default parameters.
#' @param raster Logical. (from \code{\link[ComplexHeatmap]{Heatmap}} documentation) If \code{TRUE}, the function will render the heatmap body as a raster image. It helps to reduce file size when the matrix is huge.
#' @param raster.quality Numeric. The quality of the raster image. A higher value will slow rendering but will lower expression smoothing. Ignored if \code{raster} = \code{FALSE}.
#' @param output.data Logical. If \code{TRUE}, the function will return a \code{matrix} object of the cell expression data, scaled or not, instead of displaying anything.
#' @param ... Additional arguments to be passed to \code{\link[ComplexHeatmap]{Heatmap}}, such as \code{show_parent_dend_line}, \code{clustering_method_rows}, etc, accepts any parameter that wasn't already internally passed to \code{\link[ComplexHeatmap]{Heatmap}} (for example, \code{outer.border} sets the \code{border} parameter of \code{\link[ComplexHeatmap]{Heatmap}}, so you will get an error if you try to pass the \code{border} parameter in \code{\link[RightOmicsTools]{Cell_Heatmap}}).
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap}} object, or a \code{matrix} object of the cell expression data, scaled or not.
#'
#' @examples
#' \dontshow{
#' suppressWarnings(suppressPackageStartupMessages(library(Seurat)))
#' }
#' # Prepare data
#' pbmc3k <- Right_Data("pbmc3k")
#' \dontshow{
#' pbmc3k = suppressWarnings(suppressPackageStartupMessages(subset(pbmc3k, idents = c("DC", "Platelets"), invert = TRUE)))
#' }
#' pbmc3k.markers = c("CCR7", "TCF7", "S100A9", "CD14",
#'                  "CD40LG", "CD2", "CD79A", "TCL1A",
#'                  "CCL5", "CD8A", "CDKN1C", "MS4A4A",
#'                  "GNLY", "GZMB")
#'
#' # Example 1: default parameters
#' Cell_Heatmap(pbmc3k,
#'              features = pbmc3k.markers)
#' @import Seurat
#' @import SeuratObject
#' @import scales
#' @import grid
#' @import ComplexHeatmap
#' @import grDevices
#' @import colorRamp2
#' @importFrom stats setNames quantile hclust dist
#' @export

Cell_Heatmap = function(seurat_object,
                        assay = "RNA",
                        layer = "data",
                        features,
                        split.by = NULL,
                        idents = NULL,
                        split.idents = NULL,
                        scale = TRUE,
                        rescale = FALSE,
                        rescale.range = c(0, 3),
                        rotate.axis = FALSE,
                        col.min = ifelse(isTRUE(scale), -2, 0),
                        col.max = ifelse(isTRUE(scale), 2, "q100"),
                        data.colors = if (isTRUE(scale)) c("#35A5FF","white","red") else "Viridis",
                        palette.reverse = FALSE,
                        na.color = "grey40",
                        idents.colors = NULL,
                        show.idents.names.colors = TRUE,
                        show.idents.oppo.colors = FALSE,
                        split.colors = NULL,
                        show.split.names.colors = TRUE,
                        show.split.oppo.colors = FALSE,
                        order.idents = NULL,
                        order.split = NULL,
                        order.colors = TRUE,
                        kmeans.repeats = 100,
                        shuffle.cells = TRUE,
                        cluster.cells = FALSE,
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
                        column.names.side = "top",
                        column.names.height = 15,
                        outer.border = TRUE,
                        data.legend.name = ifelse(isTRUE(scale),"Z-Score","Expression"),
                        data.legend.side = "bottom",
                        data.legend.direction = "horizontal",
                        data.legend.position = "topcenter",
                        data.legend.width = 5,
                        idents.legend.name = "Clusters",
                        show.idents.legend = TRUE,
                        split.legend.name = split.by,
                        show.split.legend = TRUE,
                        legend.title.size = 10,
                        legend.text.size = 10,
                        legend.gap = 10,
                        raster = ifelse(ncol(seurat_object) > 3000, TRUE, FALSE),
                        raster.quality = 10,
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
  cell.idents = data$ident
  names(cell.idents) = rownames(data)
  data$ident = NULL
  mat = as.matrix(data)
  mat = mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  cell.idents = cell.idents[names(cell.idents) %in% rownames(mat)]
  if (length(unique(cell.idents)) < 2) {
    stop("Less than two identities left, cannot compare features expression")
  }
  idents.removed = setdiff(ident.3, unique(cell.idents))
  mat = mat[ , colSums(mat) > 0, drop = FALSE]
  if (ncol(mat) == 0) {
    stop("None of the features were expressed in any cells")
  }
  features.removed = c(features.removed,setdiff(colnames(data), colnames(mat)))
  if (length(features.removed) > 0) {
    message("The following features were removed as they were not found or were not expressed in any cells:\n",paste0(features.removed, collapse = ", "))
  }
  if (isTRUE(scale)) {
    mat = scale(mat)
  }
  if (isTRUE(rescale) & isFALSE(scale)) {
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

  mat = cbind(cell.idents, as.data.frame(mat))

  mats = list()
  low.cell.call = FALSE
  for (i in ident.3[!ident.3 %in% idents.removed]) {
    mats[[i]] = mat[mat$cell.idents == i, , drop = FALSE]
    rownames(mats[[i]]) = paste(mats[[i]]$cell.idents, rownames(mats[[i]]), sep = "_")
    mats[[i]]$cell.idents = NULL
    if (isTRUE(shuffle.cells) & isFALSE(cluster.cells)) {
      mats[[i]] = t(mats[[i]][sample(nrow(mats[[i]])), , drop = FALSE])
    }
    else {
      if (isTRUE(cluster.cells)) {
      hc = hclust(dist(mats[[i]]))
      mats[[i]] = mats[[i]][hc$order, , drop = FALSE]
      }
      mats[[i]] = t(mats[[i]])
    }
    if (ncol(mats[[i]]) < 0.05*nrow(mat) & isFALSE(output.data) & isFALSE(low.cell.call) & length(ident.3[!ident.3 %in% idents.removed]) <= 20) {
      message("One or more identities account for less than 5% of cells.\nYou might experience trouble visualizing them.\nUsing DotPlot_Heatmap() might be more appropriate.")
      low.cell.call = TRUE
    }
  }

  mat = do.call(cbind, mats)

  cell.idents = colnames(mat)

  if (isTRUE(rotate.axis)) {
    mat = t(mat)
  }

  if (isTRUE(output.data)) {
    return(round(mat, 4))
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
    dup.colors = list()
    dup.colors2 = list()
    for (i in idents.colors) {
      dup.colors[[i]] = rep(i, length(ident.3)/length(ident.1))
      dup.colors2[[i]] = split.colors
    }
    idents.colors = unlist(dup.colors)
    split.colors = unlist(dup.colors2)
    names(idents.colors) = ident.3
    names(split.colors) = ident.3
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
                           border = TRUE)
  }
  if (is.character(split.by) & isTRUE(show.split.legend) & (isTRUE(show.split.names.colors) | isTRUE(show.split.oppo.colors))) {
    split.legend = Legend(at = ident.2,
                          legend_gp = gpar(fill = split.colors2),
                          title = split.legend.name,
                          gap = unit(0.5, "cm"),
                          border = TRUE)
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

  if (sum(is.na(mat)) > 0) {
    na.label = list(Legend(at = "NA",
                           legend_gp = gpar(fill = na.color),
                           title = NULL,
                           gap = unit(0.5, "cm"),
                           border = TRUE))
  }
  else {
    na.label = list()
  }
  anno.legend = c(na.label,anno.legend)

  if (length(mats) > 20) {
    message("More than 20 identities will be displayed.\nYou might experience excessive raster smoothing and/or difficult visualization of small clusters.\nUsing DotPlot_Heatmap() might be more appropriate.")
  }

  cells.color = split.color = idents.slices = NULL
  slices.color = ifelse(idents.names.size > 0, "black", "white")
  kmeans.color = ifelse(features.kmeans.numbers.size > 0, "black", "white")
  for (i in 1:length(mats)) {
    cells.color = c(cells.color,rep(idents.colors[i], ncol(mats[[i]])))
    idents.slices = c(idents.slices, rep(names(mats)[i], ncol(mats[[i]])))
    if (is.character(split.by)) {
      split.color = c(split.color,rep(split.colors[i], ncol(mats[[i]])))
    }
  }
  if (isFALSE(rotate.axis)) {
    annotation.side = "column"
  }
  else {
    annotation.side = "row"
  }
  idents.labels = idents.labels2 = NULL
  if (is.character(split.by)) {
    if (isTRUE(show.idents.names.colors) & isTRUE(show.split.names.colors)) {
      if ((column.names.side == "top" & isFALSE(rotate.axis)) | (row.names.side == "left" & isTRUE(rotate.axis))) {
        idents.labels = HeatmapAnnotation(
          Identities = cell.idents,
          Split = cell.idents,
          col = list(Identities = setNames(cells.color, cell.idents),
                     Split = setNames(split.color, cell.idents)),
          na_col = na.color,
          which = annotation.side,
          show_annotation_name = FALSE,
          show_legend = FALSE)
      }
      else {
        idents.labels = HeatmapAnnotation(
          Split = cell.idents,
          Identities = cell.idents,
          col = list(Split = setNames(split.color, cell.idents),
                     Identities = setNames(cells.color, cell.idents)),
          na_col = na.color,
          which = annotation.side,
          show_annotation_name = FALSE,
          show_legend = FALSE)
      }
    }
    else if (isTRUE(show.idents.names.colors)) {
      idents.labels = HeatmapAnnotation(
        Identities = cell.idents,
        col = list(Identities = setNames(cells.color, cell.idents)),
        na_col = na.color,
        which = annotation.side,
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else if (isTRUE(show.split.names.colors)) {
      idents.labels = HeatmapAnnotation(
        Split = cell.idents,
        col = list(Split = setNames(split.color, cell.idents)),
        na_col = na.color,
        which = annotation.side,
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else {
      idents.labels = NULL
    }
    if (isTRUE(show.idents.oppo.colors) & isTRUE(show.split.oppo.colors)) {
      if ((column.names.side == "top" & isFALSE(rotate.axis)) | (row.names.side == "left" & isTRUE(rotate.axis))) {
        idents.labels2 = HeatmapAnnotation(
          Split = cell.idents,
          Identities = cell.idents,
          col = list(Split = setNames(split.color, cell.idents),
                     Identities = setNames(cells.color, cell.idents)),
          na_col = na.color,
          which = annotation.side,
          show_annotation_name = FALSE,
          show_legend = FALSE)
      }
      else {
        idents.labels2 = HeatmapAnnotation(
          Identities = cell.idents,
          Split = cell.idents,
          col = list(Identities = setNames(cells.color, cell.idents),
                     Split = setNames(split.color, cell.idents)),
          na_col = na.color,
          which = annotation.side,
          show_annotation_name = FALSE,
          show_legend = FALSE)
      }
    }
    else if (isTRUE(show.idents.oppo.colors)) {
      idents.labels2 = HeatmapAnnotation(
        Identities = cell.idents,
        col = list(Identities = setNames(cells.color, cell.idents)),
        na_col = na.color,
        which = annotation.side,
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else if (isTRUE(show.split.oppo.colors)) {
      idents.labels2 = HeatmapAnnotation(
        Split = cell.idents,
        col = list(Split = setNames(split.color, cell.idents)),
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
        Identities = cell.idents,
        col = list(Identities = setNames(cells.color, cell.idents)),
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
        Identities = cell.idents,
        col = list(Identities = setNames(cells.color, cell.idents)),
        na_col = na.color,
        which = annotation.side,
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else {
      idents.labels2 = NULL
    }
  }

  row.dend.side = row.title.side = "right"
  column.dend.side = column.title.side = "bottom"
  idents.names.gp = gpar(fontsize = idents.names.size, col = slices.color)
  features.names.gp = gpar(fontsize = features.names.size, fontface = features.names.style)
  features.title.gp = gpar(fontsize = features.kmeans.numbers.size, col = kmeans.color)
  if ((column.names.side == "bottom" & isFALSE(rotate.axis)) | (row.names.side == "right" & isTRUE(rotate.axis))) {
    idents.labels2tmp = idents.labels2
    idents.labels2 = idents.labels
    idents.labels = idents.labels2tmp
  }
  if (isFALSE(rotate.axis)) {
    column.title.side = column.names.side
    if (row.names.side == "right") {
      row.dend.side = row.title.side = "left"
    }
  }
  else {
    row.title.side = row.names.side
    if (column.names.side == "bottom") {
      column.dend.side = column.title.side = "top"
    }
  }

  ht = NULL
  ht_opt$ANNOTATION_LEGEND_PADDING = unit(legend.gap, "mm")
  if (data.legend.side == "right" | data.legend.side == "left" | (isTRUE(rotate.axis)) & column.names.side == "top") {
    ht_opt$HEATMAP_LEGEND_PADDING = unit(legend.gap, "mm")
  }
  else {
    ht_opt$HEATMAP_LEGEND_PADDING = unit(2, "mm")
  }

  if (isFALSE(rotate.axis)) {
    ht = Heatmap(
      mat,
      col = data.colors,
      na_col = na.color,
      name = data.legend.name,
      border = outer.border,
      top_annotation = idents.labels,
      bottom_annotation = idents.labels2,
      cluster_rows = cluster.features,
      column_split = factor(idents.slices, levels = names(mats)),
      row_km = features.kmeans,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      show_row_names = ifelse(features.names.size > 0, TRUE, FALSE),
      column_title_gp = idents.names.gp,
      row_names_gp = features.names.gp,
      row_title_rot = 0,
      row_title_gp = features.title.gp,
      column_title_rot = column.names.angle,
      row_names_max_width = unit(row.names.width, "cm"),
      column_names_max_height = unit(column.names.height, "cm"),
      row_names_side = row.names.side,
      row_dend_side = row.dend.side,
      row_title_side = row.title.side,
      column_dend_side = column.dend.side,
      column_title_side = column.title.side,
      row_gap = unit(features.gap, "mm"),
      column_gap = unit(idents.gap, "mm"),
      use_raster = raster,
      raster_quality = raster.quality,
      heatmap_legend_param = list(
        legend_direction = data.legend.direction,
        title_position = data.legend.position,
        legend_width = unit(data.legend.width, "cm"),
        title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
        labels_gp = gpar(fontsize = legend.text.size)),
      ...)
  }
  else {
    ht = Heatmap(
      mat,
      col = data.colors,
      na_col = na.color,
      name = data.legend.name,
      border = outer.border,
      left_annotation = idents.labels,
      right_annotation = idents.labels2,
      cluster_rows = FALSE,
      row_split = factor(idents.slices, levels = names(mats)),
      column_km = features.kmeans,
      cluster_columns = cluster.features,
      show_row_names = FALSE,
      show_column_names = ifelse(features.names.size > 0, TRUE, FALSE),
      row_title_gp = idents.names.gp,
      column_names_gp = features.names.gp,
      row_title_rot = 0,
      column_title_gp = features.title.gp,
      column_names_rot = column.names.angle,
      row_names_max_width = unit(row.names.width, "cm"),
      column_names_side = column.names.side,
      column_names_max_height = unit(column.names.height, "cm"),
      row_names_side = row.names.side,
      row_dend_side = row.dend.side,
      row_title_side = row.title.side,
      column_dend_side = column.dend.side,
      column_title_side = column.title.side,
      row_gap = unit(idents.gap, "mm"),
      column_gap = unit(features.gap, "mm"),
      use_raster = raster,
      raster_quality = raster.quality,
      heatmap_legend_param = list(
        legend_direction = data.legend.direction,
        title_position = data.legend.position,
        legend_width = unit(data.legend.width, "cm"),
        title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
        labels_gp = gpar(fontsize = legend.text.size)),
      ...)
  }

  return(draw(ht,
              heatmap_legend_side = data.legend.side,
              align_heatmap_legend = "heatmap_center",
              align_annotation_legend = "heatmap_center",
              annotation_legend_list = anno.legend,
              legend_grouping = "original"))
}
