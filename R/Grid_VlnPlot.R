#' @title Grid_VlnPlot
#'
#' @description This function is a stacked violin plot optimized to display features expression in a Seurat object in a grid fashion (square) instead of a single column like other stacked violin functions available in other packages, resulting in nicer plots and easier to include in publications.
#'
#' @param seurat_object A Seurat object.
#' @param assay Character. If the Seurat object contains multiple RNA assays, you may specify which one to use (for example "RNA2" if you have created a second RNA assay you named "RNA2". See Seurat v5 vignettes for more information). You may also use another assay such as SCT to pull gene expression from.
#' @param layer Character. Formerly known as slot. It is recommended to use 'data'.
#' @param features Character. A vector of features to plot.
#' @param idents Character. A vector with one or several identities in the active.ident identity to use if you only want those (instead of subsetting your object). If NULL, all identities will be used.
#' @param rotate.axis Logical. If TRUE, flips the axis, effectively displaying violins vertically instead of horizontally.
#' @param colors Character. A vector of colors to use for the active.ident identity, of same length as the number of identities in the active.ident identity or supplied to the idents parameter. If NULL, uses Seurat's default colors.
#' @param order.idents Character. A vector specifying either "reverse" or the levels of the active.ident identity to order the cells.
#' @param order.colors Logical. If TRUE, the colors will automatically be ordered according to order.idents. Ignored if order.idents is NULL.
#' @param idents.text.size Numeric. The size of the axis identities.
#' @param show.idents Logical. If TRUE, shows the identities on the plot.
#' @param features.text.size Numeric. The size of the axis features.
#' @param legend.text.size Numeric. The size of the legend text. Ignored if legend is FALSE.
#' @param legend.position Character. Which side to display the legend, "left", "right", "top" or "bottom". Ignored if legend is FALSE.
#' @param legend Logical. If TRUE, shows the legend.
#' @param ncol Numeric. Number of columns to use. If "square", will display features in a square grid or as close as possible depending on number of features.
#'
#' @return A ggplot object.
#'
#' @import Seurat
#' @import SeuratObject
#' @import data.table
#' @import ggplot2
#' @import grDevices
#' @export

Grid_VlnPlot = function(seurat_object,
                        assay = "RNA",
                        layer = "data",
                        features = NULL,
                        idents = NULL,
                        rotate.axis = FALSE,
                        colors = NULL,
                        order.idents = NULL,
                        order.colors = TRUE,
                        idents.text.size = 9,
                        show.idents = FALSE,
                        features.text.size = 11,
                        legend.text.size = 12,
                        legend.position = "bottom",
                        legend = TRUE,
                        ncol = "square") {

  if (!is.character(features)) {
    stop("Please provide features to plot.")
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

  DefaultAssay(seurat_object) = assay
  data = FetchData(object = seurat_object, vars = c("ident",features), layer = layer)
  data = data[data$ident %in% ident.1, ]
  data = melt(setDT(data), variable.name = "gene", value.name = "expression", id.vars = 1)
  if (!is.character(colors)) {
    SeuratColors = function(n = 6, h = c(0, 360) + 15){
      if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
      hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
    }
    colors = SeuratColors(n = length(ident.1))
  }
  if (isTRUE(order.colors)) {
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
      facet_wrap(~gene, ncol=ncol, scales="free_x") + labs(y="", x="", fill ="") + theme_bw()  +
      theme(legend.position=legend.position, panel.spacing = unit(0, "lines"), legend.text = element_text(size = legend.text.size),
            axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text.x = element_text(size = features.text.size),
            panel.grid = element_blank()) +
      scale_fill_manual(values = colors)
  }
  else {
    gg = ggplot(data, aes(y=expression, x=ident, fill=ident)) + geom_violin(scale = "width", trim = T, adjust = 1) +
      facet_wrap(~gene, ncol=ncol, scales="free_y") + labs(y="", x="", fill ="") + theme_bw()  +
      theme(legend.position=legend.position, panel.spacing = unit(0, "lines"), legend.text = element_text(size = legend.text.size),
            axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text.x = element_text(size = features.text.size),
            panel.grid = element_blank()) +
      scale_fill_manual(values = colors)
  }
  if (isFALSE(show.idents) & isFALSE(rotate.axis)) {
    gg = gg + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }
  if (isFALSE(show.idents) & isTRUE(rotate.axis)) {
    gg = gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  if (isFALSE(legend)) {
    gg = gg + theme(legend.position="none")
  }
  return(gg)
}
