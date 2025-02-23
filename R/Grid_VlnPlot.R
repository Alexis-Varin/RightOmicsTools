#' @title Square stacked violin plot of gene expression in each identity
#'
#' @description This function is a stacked violin plot optimized to display the average expression of features in a \pkg{Seurat} object in a grid fashion (square) instead of a single column like other stacked violin functions available in other packages.
#'
#' @param seurat_object A \pkg{Seurat} object.
#' @param assay Character. The name of an assay containing the \code{layer} with the expression matrix. If the \code{seurat_object} contains multiple 'RNA' assays, you may specify which one to use (for example, 'RNA2' if you have created a second 'RNA' assay you named 'RNA2'. See \href{https://satijalab.org/seurat/articles/seurat5_essential_commands.html#create-seurat-or-assay-objects}{Seurat v5 vignettes} for more information). You may also use another assay, such as 'SCT', to pull feature expression from.
#' @param layer Character. The name of a layer (formerly known as slot) which stores the expression matrix. It is recommended to use 'data'.
#' @param features Character. The names of one or several features to plot the cell expression from.
#' @param idents Character. The names of one or several identities in the active.ident metadata to select. If \code{NULL}, all identities are used.
#' @param scale Logical. If \code{TRUE}, scales the violins to have the same max height between features.
#' @param rotate.axis Logical. If \code{TRUE}, flips the axis, displaying violins vertically instead of horizontally.
#' @param colors Character. The color names for each identity of the active.ident metadata or in \code{idents}. If \code{NULL}, uses \pkg{Seurat}'s default colors.
#' @param order.idents Character or Numeric. Either 'reverse', or the identities (as names or as numeric values corresponding to the indices) of the active.ident metadata or in \code{idents} to order the cells.
#' @param order.colors Logical. If \code{TRUE}, the \code{colors} will automatically be ordered according to \code{order.prop}. Ignored if \code{order.prop} = \code{NULL}.
#' @param idents.text.size Numeric. The font size of the identity names. Ignored if \code{show.idents} = \code{FALSE}.
#' @param show.idents Logical. If \code{TRUE}, shows the identity names on the plot.
#' @param features.text.size Numeric. The font size of the feature names.
#' @param legend.text.size Numeric. The font size of the legend text. Ignored if \code{show.legend} = \code{FALSE}.
#' @param legend.side Character. The side where the legend will be displayed, either 'left', 'right', 'top' or 'bottom'. Ignored if \code{show.legend} = \code{FALSE}.
#' @param show.legend Logical. If \code{TRUE}, shows the legend.
#' @param ncol Numeric. The number of columns to use. If 'square', will display features in a square grid or as close as possible depending on the number of features.
#'
#' @return A \pkg{ggplot2} object.
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
#' Grid_VlnPlot(pbmc3k,
#'              features = pbmc3k.markers)
#' @import Seurat
#' @import SeuratObject
#' @import data.table
#' @import ggplot2
#' @import scales
#' @import grDevices
#' @export

Grid_VlnPlot = function(seurat_object,
                        assay = "RNA",
                        layer = "data",
                        features,
                        idents = NULL,
                        scale = TRUE,
                        rotate.axis = FALSE,
                        colors = NULL,
                        order.idents = NULL,
                        order.colors = TRUE,
                        idents.text.size = 9,
                        show.idents = FALSE,
                        features.text.size = 9,
                        legend.text.size = 7,
                        legend.side = "bottom",
                        show.legend = TRUE,
                        ncol = "square") {

  ident = NULL

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

  DefaultAssay(seurat_object) = assay
  data = suppressWarnings(FetchData(object = seurat_object, vars = c("ident",features), layer = layer))
  if (ncol(data) == 1) {
    stop("None of the features were found or expressed in any cells")
  }
  features.removed = setdiff(features, colnames(data))
  data = data[data$ident %in% ident.1, ]
  data = melt(setDT(data), variable.name = "gene", value.name = "expression", id.vars = 1)
  data$ident = factor(data$ident, levels = ident.1)
  if (!is.character(colors)) {
    colors = hue_pal()(n = length(ident.1))
  }
  if (isTRUE(order.colors)) {
    if (!is.null(names(colors))) {
      colors = colors[ident.1]
    }
    if (is.character(order.idents)) {
      if (length(order.idents) > 1) {
        names(colors) = ident.1
        colors = colors[order.idents]
      }
      else {
        if (order.idents == "reverse") {
          colors = rev(colors)
        }
      }
    }
  }
  if (ncol == "square") {
    ncol = ceiling(sqrt(length(levels(data$gene))))
  }
  if (isFALSE(rotate.axis)) {
    gg = ggplot(data, aes(y=ident, x=expression, fill=ident)) + geom_violin(scale = "width", trim = T, adjust = 1) +
      labs(y="", x="", fill ="") + theme_bw()  +
      theme(legend.position = legend.side, panel.spacing = unit(0, "lines"), legend.text = element_text(size = legend.text.size),
            axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.x = element_text(size = features.text.size),
            panel.grid = element_blank()) +
      scale_fill_manual(values = colors)
    if (isTRUE(scale)) {
      gg = gg + facet_wrap(~gene, ncol=ncol, scales="free_x")
    }
    else {
      gg = gg + facet_wrap(~gene, ncol=ncol)
    }
  }
  else {
    gg = ggplot(data, aes(y=expression, x=ident, fill=ident)) + geom_violin(scale = "width", trim = T, adjust = 1) +
      labs(y="", x="", fill ="") + theme_bw()  +
      theme(legend.position = legend.side, panel.spacing = unit(0, "lines"), legend.text = element_text(size = legend.text.size),
            axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text.x = element_text(size = features.text.size),
            panel.grid = element_blank()) +
      scale_fill_manual(values = colors)
    if (isTRUE(scale)) {
      gg = gg + facet_wrap(~gene, ncol=ncol, scales="free_y")
    }
    else {
      gg = gg + facet_wrap(~gene, ncol=ncol)
    }
  }
  if (isFALSE(show.idents) & isFALSE(rotate.axis)) {
    gg = gg + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }
  if (isFALSE(show.idents) & isTRUE(rotate.axis)) {
    gg = gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  if (isFALSE(show.legend)) {
    gg = gg + theme(legend.position="none")
  }
  return(gg)
}
