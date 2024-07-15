#' @title Cell_Heatmap
#'
#' @description This function generates a dotplot or a heatmap to visualize features expression in a Seurat object. Credits to Seurat's dev team for the original DotPlot() function from which data processing of this function is derived from and to https://divingintogeneticsandgenomics.com/post/clustered-dotplot-for-single-cell-rnaseq/ for the initial idea to use ComplexHeatmap to draw a dotplot and the layer_fun function that draws the dots. Slight improvements were implemented here and there for my personal use, and made available to all through this package.
#'
#' @param seurat_object A Seurat object.
#' @param assay Character. If the Seurat object contains multiple RNA assays, you may specify which one to use (for example "RNA2" if you have created a second RNA assay you named "RNA2". See Seurat v5 vignettes for more information). You may also use another assay such as SCT to pull gene expression from.
#' @param layer Character. Formerly known as slot. It is recommended to use 'data'.
#' @param features Character. A vector of features to plot.
#' @param split.by Character. The name of an identity in the metadata slot to split the active.ident identity by.
#' @param idents Character. A vector with one or several identities in the active.ident identity to use if you only want those (instead of subsetting your object). If NULL, all identities will be used.
#' @param split.idents Character. A vector with one or several identities in the split.by identity to use if you only want those. If NULL, all identities will be used.
#' @param scale Logical. If TRUE, average expression values will be scaled using R base scale() function and default parameters. The resulting values will be Z-scores (mean subtracted values divided by standard deviation) and not positive average expression values anymore, which is why there will be positive and negative values displayed, depending on if the average expression in a particular identity is below or above the mean average expression from all identities (which is calculated independently for each feature). Caution should be exercised when interpreting results with low number of identities (typically below 5), as small differences in average expression might lead to exacerbated differences when scaled.
#' @param col.min Character or Numeric. The minimum value for the color scale parameter internally passed to colorRamp2::colorRamp2(). If character, must be a quantile in the form 'qX' where X is a number between 0 and 100. A value of 'q5' or 'q10' is useful to reduce the effect of outlier values (e.g. a very low value that significantly alters the color scale range of all other values).
#' @param col.max Character or Numeric. The maximum value for the color scale parameter internally passed to colorRamp2::colorRamp2(). If character, must be a quantile in the form 'qX' where X is a number between 0 and 100. A value of 'q95' or 'q90' is useful to reduce the effect of outlier values (e.g. a very high value that significantly alters the color scale range of all other values).
#' @param data.colors Character. Either a character vector of exactly 3 colors, corresponding to the lowest, zero (or middle if scale = FALSE), and highest values of the expression matrix and internally passed to colorRamp2::colorRamp2(), or a single character value corresponding to the name of a palette and internally passed to the hcl_palette parameter of colorRamp2::colorRamp2() (such as 'Inferno', 'Berlin', 'Viridis' etc, check grDevices::hcl.pals() for all palettes available).
#' @param palette.reverse Logical. If TRUE and if data.colors is a palette (such as "Viridis"), the function will reverse its colors.
#' @param na.color Character. The color to use for NA values.
#' @param idents.colors Character. A vector of colors to use for the active.ident identity, of same length as the number of identities in the active.ident identity or supplied to the idents parameter. If NULL, uses Seurat's default colors.
#' @param show.idents.names.colors Logical. If TRUE, the function will display the colors specified by the idents.colors parameter next to identity names.
#' @param show.idents.dend.colors Logical. If TRUE, the function will display the colors specified by the idents.colors parameter next to the dendogram. Ignored if cluster.rows and cluster.columns are set to FALSE.
#' @param split.colors Character. A vector of colors to use for the split.by identity, of same length as the number of identities in the split.by identity or supplied to the split.idents parameter. If NULL, uses a custom set of colors from grDevices::colors(). Ignored if split.by = NULL.
#' @param show.split.names.colors Logical. If TRUE, the function will display the colors specified by the split.colors parameter next to identity names. Ignored if split.by = NULL.
#' @param show.split.dend.colors Logical. If TRUE, the function will display the colors specified by the split.colors parameter next to the dendogram. Ignored if split.by = NULL or if cluster.rows and cluster.columns are set to FALSE.
#' @param order.idents Character. A vector specifying either "reverse" or the levels of the active.ident identity to order the cells.
#' @param order.split Character. A vector specifying either "reverse" or the levels of the split.by identity to order the cells. Ignored if split.by is NULL.
#' @param order.colors Logical. If TRUE, the colors for idents and split.idents will automatically be ordered according to order.idents and order.split. Ignored if order.idents and order.split are NULL.
#' @param kmeans.repeats Numeric. The number of k-means runs to get a consensus k-means clustering. Ignored if cluster.rows and cluster.columns are set to FALSE.
#' @param cluster.rows Logical or Function. If TRUE, the function will cluster the rows. You may also pass an hclust or dendogram object which contains clustering.
#' @param row.kmeans Numeric. The number of k-means slices to use for row clustering. Ignored if cluster.rows = FALSE.
#' @param row.names.side Character. The side where the row names will be displayed, either 'left' or 'right'. If cluster.rows = TRUE or Function, the dendogram will be displayed on the opposite side.
#' @param row.names.width Numeric. The width of the row names. Increase this parameter if your row names are truncated.
#' @param cluster.columns Logical or Function. If TRUE, the function will cluster the columns. You may also pass an hclust or dendogram object which contains clustering.
#' @param column.kmeans Numeric. The number of k-means slices to use for column clustering. Ignored if cluster.columns = FALSE.
#' @param column.names.rotation Numeric. The angle of rotation for the column names.
#' @param column.names.side Character. The side where the column names will be displayed, either 'top' or 'bottom'. If cluster.columns = TRUE or Function, the dendogram will be displayed on the opposite side.
#' @param column.names.height Numeric. The height of the column names. Increase this parameter if your column names are truncated.
#' @param outer.border Logical. If TRUE, the function will display an outer border around the plot or around each slice if row and/or column k-means is higher than 1.
#' @param data.legend.name Character. The name of the data legend.
#' @param data.legend.side Character. The side where the data legend will be displayed, either 'left', 'right', 'top' or 'bottom'.
#' @param data.legend.direction Character. The direction of the data legend, either 'horizontal' or 'vertical'.
#' @param data.legend.position Numeric. The centering of the data legend name, there are many options, default option from ComplexHeatmap::Heatmap() is 'topleft'.
#' @param data.legend.width Numeric. How long the data legend will be, only affects the data legend if data.legend.direction = 'horizontal'.
#' @param idents.legend.name Character. The name of the active.ident identity legend.
#' @param show.idents.legend Logical. If TRUE, the function will display a legend for the active.ident identity.
#' @param split.legend.name Character. The name of the split.by identity legend. Ignored if split.by = NULL.
#' @param show.split.legend Logical. If TRUE, the function will display a legend for the split.by identity.
#' @param output.data Logical. If TRUE, the function will return a list composed of the expression matrix and another matrix containing the percent of cells expressing each feature, instead of drawing.
#' @param ... Additional arguments to pass to the ComplexHeatmap::Heatmap() function, such as column_names_gp, clustering_method_columns, etc, accepts any parameter that wasn't already internally passed to ComplexHeatmap::Heatmap() (for example, cluster.columns sets the cluster_columns parameter of the inner function, so you will get an error if you try to pass it again).
#'
#' @return A ComplexHeatmap object, either as a dotplot, or a heatmap, or a list containing a matrix of the expression data and another matrix containing the percent of cells expressing each feature.
#'
#' @import Seurat
#' @import SeuratObject
#' @import grid
#' @import stats
#' @import ComplexHeatmap
#' @import grDevices
#' @import colorRamp2

Cell_Heatmap = function(seurat_object,
                           assay = "RNA",
                           layer = "data",
                           features,
                           split.by = NULL,
                           idents = NULL,
                           split.idents = NULL,
                           scale = TRUE,
                           col.min = ifelse(isTRUE(scale), -2, 0),
                           col.max = ifelse(isTRUE(scale), 2, "q100"),
                           data.colors = if (isTRUE(scale)) c("#35A5FF","white","red") else "Viridis",
                           palette.reverse = FALSE,
                           na.color = "black",
                           idents.colors = NULL,
                           show.idents.names.colors = TRUE,
                           show.idents.dend.colors = TRUE,
                           split.colors = NULL,
                           show.split.names.colors = TRUE,
                           show.split.dend.colors = TRUE,
                           order.idents = NULL,
                           order.split = NULL,
                           order.colors = TRUE,
                           kmeans.repeats = 100,
                           cluster.rows =TRUE,
                           row.kmeans = 1,
                           row.names.side = "left",
                           row.names.width = unit(15, "cm"),
                           cluster.columns = TRUE,
                           column.kmeans = 1,
                           column.names.rotation = 45,
                           column.names.side = "bottom",
                           column.names.height = unit(15, "cm"),
                           outer.border = TRUE,
                           data.legend.name = ifelse(isTRUE(scale),"Z-Score","Average Expression"),
                           data.legend.side = "bottom",
                           data.legend.direction = "horizontal",
                           data.legend.position = "topcenter",
                           data.legend.width = 5,
                           idents.legend.name = "Active identities",
                           show.idents.legend = TRUE,
                           split.legend.name = "Split identities",
                           show.split.legend = TRUE,
                           output.data = FALSE,
                           ...) {

  if (!is.character(split.by)) {
    ident = "ident"
  }
  else {
    ident = c("ident",split.by)
  }

  if (!is.character(idents)) {
    ident.1 = levels(Idents(seurat_object))
  }
  else {
    ident.1 = idents
  }
  if (!is.null(order.idents)) {
    if (is.character(order.idents)) {
      if (length(order.idents) > 1) {
        ident.1 = ident.1[order.idents]
      }
      else {
        if (order.idents == "reverse") {
          ident.1 = rev(ident.1)
        }
        else {
          stop("order.idents needs to be either 'reverse' or a character vector")
        }
      }
    }
    else {
      stop("order.idents needs to be either 'reverse' or a character vector")
    }
  }

  if (is.character(split.by)) {
    vars = c("ident",split.by,features)
  }
  else {
    vars = c("ident",features)
  }
  DefaultAssay(seurat_object) = assay
  data = suppressWarnings(FetchData(object = seurat_object, vars = vars, layer = layer))
  if (is.character(split.by)) {
    if (!is.character(split.idents)) {
      Idents(seurat_object) = split.by
      ident.2 = levels(Idents(seurat_object))
    }
    else {
      ident.2 = split.idents
    }
    if (!is.null(order.split)) {
      if (is.character(order.split)) {
        if (length(order.split) > 1) {
          ident.1 = ident.2[order.split]
        }
        else {
          if (order.split == "reverse") {
            ident.2 = rev(ident.2)
          }
          else {
            stop("order.split needs to be either 'reverse' or a character vector")
          }
        }
      }
      else {
        stop("order.split needs to be either 'reverse' or a character vector")
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
  mat = round(as.matrix(data), 3)
  mat = mat[rowSums(is.na(mat)) < ncol(mat), ]
  cell.idents = cell.idents[names(cell.idents) %in% rownames(mat)]
  idents.removed = setdiff(ident.3, unique(cell.idents))
  mat = mat[ , colSums(mat) > 0]
  removed = setdiff(colnames(data), colnames(mat))
  if (length(removed) > 0) {
    warning("The following features were removed as they were not expressed in any cells:\n ",paste0(" ",removed, collapse = ","))
  }
  if (isTRUE(scale)) {
    mat = scale(mat)
  }

  mat = cbind(cell.idents,as.data.frame(mat))

  mats = list()
  for (i in ident.3[!ident.3 %in% idents.removed]) {
    mats[[i]] = mat[mat$cell.idents == i,]
    mats[[i]]$cell.idents = NULL
    mats[[i]] = t(mats[[i]][sample(nrow(mats[[i]])), ])
  }

  if (isTRUE(output.data)) {
    return(mats)
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

  if (isTRUE(scale)) {
    col.mid = 0
  }
  else {
    col.mid = (col.min+col.max)/2
  }
  if (length(data.colors) == 3) {
    data.colors = colorRamp2(breaks = c(col.min,col.mid,col.max),
                             colors = c(data.colors[1],data.colors[2],data.colors[3]))
  }
  else {
    data.colors = colorRamp2(breaks = c(col.min,col.max),
                             hcl_palette = data.colors, reverse = palette.reverse)
  }

  if (!is.character(split.colors) & is.character(split.by)) {
    current.seed = .Random.seed
    set.seed(1234)
    CustomColors = colors()[grep('gr(a|e)y|white|black', colors(), invert = T)]
    split.colors = sample(CustomColors, length(ident.2))
    .Random.seed = current.seed
  }

  if (!is.character(idents.colors)) {
    SeuratColors = function(n = 6, h = c(0, 360) + 15){
      if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
      hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
    }
    idents.colors = SeuratColors(n = length(ident.1))
  }
  if (isTRUE(order.colors)) {
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
    if (is.character(order.split) & is.character(split.by)) {
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

  idents.legend = NULL
  split.legend = NULL
  if (isTRUE(show.idents.legend)) {
    idents.legend = Legend(at = ident.1,
                           legend_gp = gpar(fill = idents.colors2),
                           title = idents.legend.name,
                           gap = unit(0.5, "cm"),
                           border = TRUE)
  }
  if (is.character(split.by) & isTRUE(show.split.legend)) {
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

  if (row.names.side == "right") {
    idents.labelstemp = idents.labels
    idents.labels = idents.labels2
    idents.labels2 = idents.labelstemp
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
    message("The heatmap contains more than 20 identities, you might experience too much smoothing from raster and/or low visibility of small clusters.\nUsing DotPlot_Heatmap() might be more appropriate.")
  }

  ht = NULL
  for (i in 1:length(mats)) {
    cells.color = rep(idents.colors[i], ncol(mats[[i]]))
    if (is.character(split.by)) {
      split.color = rep(split.colors[i], ncol(mats[[i]]))
      idents.labels = HeatmapAnnotation(
        Identities = colnames(mats[[i]]),
        Split = colnames(mats[[i]]),
        col = list(Identities = setNames(cells.color, colnames(mats[[i]])),
                   Split = setNames(split.color, colnames(mats[[i]]))),
        na_col = "grey40",
        which = "column",
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    else {
      idents.labels = HeatmapAnnotation(
        Identities = colnames(mats[[i]]),
        col = list(Identities = setNames(cells.color, colnames(mats[[i]]))),
        na_col = "grey40",
        which = "column",
        show_annotation_name = FALSE,
        show_legend = FALSE)
    }
    ht = suppressWarnings(ht + suppressWarnings(Heatmap(
      mats[[i]],
      col = data.colors,
      na_col = na.color,
      name = data.legend.name,
      border = outer.border,
      #bottom_annotation = idents.labels2,
      top_annotation = idents.labels,
      cluster_rows = cluster.rows,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      column_title = names(mats)[i],
      column_title_gp = gpar(fontsize = 9),
      row_title_rot = 0,
      row_km = row.kmeans,
      row_km_repeats = kmeans.repeats,
      column_title_rot = column.names.rotation,
      row_names_side = row.names.side,
      row_names_max_width = row.names.width,
      column_names_side = column.names.side,
      column_names_max_height = column.names.height,
      heatmap_legend_param = list(
        legend_direction = data.legend.direction,
        title_position = data.legend.position,
        legend_width = unit(data.legend.width, "cm")),
      ...)))
  }

  return(draw(ht,
              heatmap_legend_side = data.legend.side,
              align_heatmap_legend = "heatmap_center",
              align_annotation_legend = "heatmap_center",
              annotation_legend_list = anno.legend,
              legend_grouping = "original",
              use_raster = TRUE,
              raster_quality = 10,
              gap = unit(0.5, "mm")))
}
