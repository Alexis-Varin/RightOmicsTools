#' @title ggplot annotation for ComplexHeatmap
#'
#' @description This function creates several plotting regions (boxes) which corresponds to each row/column of a \code{\link[ComplexHeatmap]{Heatmap}}. While \code{\link[ComplexHeatmap]{anno_zoom}} and \code{\link[ComplexHeatmap]{anno_link}} are primarily designed to add ggplot2 objects to slices (row_split, column_split, km etc), adding them to each row/column is possible but tedious, and the annotation cannot be decorated if sliced. This function is designed to function similarly to them, with the same parameters, but with a focus on adding the objects, especially \pkg{ggplot2} objects, to each row/column instead of slices, and to be decorated and sliced just like any other annotation.
#'
#' @param panel_fun Function. (from \code{\link[ComplexHeatmap]{anno_zoom}} documentation) A self-defined function that defines how to draw graphics in the box. The function must have a index argument which is the indices for the rows/columns that the box corresponds to. It can have second argument nm which is the "name" of the selected part in the heatmap.
#' @param which Character. (from \code{\link[ComplexHeatmap]{anno_zoom}} documentation) Whether it is a column annotation or a row annotation?
#' @param outer.border Logical. Whether to draw border around the annotation.
#' @param inner.border Logical. Whether to draw border around each box.
#' @param side Character. (from \code{\link[ComplexHeatmap]{anno_zoom}} documentation) Side of the boxes If it is a column annotation, valid values are "top" and "bottom"; If it is a row annotation, valid values are "left" and "right".
#' @param size Numeric or Unit. (from \code{\link[ComplexHeatmap]{anno_zoom}} documentation) The size of boxes. It can be pure numeric that they are treated as relative fractions of the total height/width of the heatmap. The value of size can also be absolute units.
#' @param gap Unit. (from \code{\link[ComplexHeatmap]{anno_zoom}} documentation) Gaps between boxes.
#' @param width Unit. (from \code{\link[ComplexHeatmap]{anno_zoom}} documentation) Width of the annotation. The value should be an absolute unit. Width is not allowed to be set for column annotation.
#' @param height Unit. (from \code{\link[ComplexHeatmap]{anno_zoom}} documentation) Height of the annotation. The value should be an absolute unit. Height is not allowed to be set for row annotation.
#'
#' @return An annotation function which can be used in HeatmapAnnotation.
#'
#' @import grid
#' @import ComplexHeatmap
#' @export

anno_ggplot = function(panel_fun = function(index, nm = NULL) { grid.rect() },
                       which = c("column", "row"),
                       outer.border = TRUE,
                       inner.border = FALSE,
                       side = ifelse(which == "column", "top", "right"),
                       size = NULL,
                       gap = unit(1, "mm"),
                       width = NULL,
                       height = NULL) {

  which = match.arg(which)

  AnnotationFunction(fun = function(index, k, n) {
    n = length(index)
    if (which == "row") index = rev(index)
    if (isTRUE(outer.border)) grid.rect()
    gg = lapply(index, panel_fun)
    lapply(seq_along(gg), function(i) {
      g = grid.grabExpr(print(gg[[i]]))
      pushViewport(viewport(x = ifelse(which == "row", 0.5, (1+2*(i-1))/(2*n)),
                            y = ifelse(which == "row", (1+2*(i-1))/(2*n), 0.5),
                            width = ifelse(which == "row", 1, 1/n),
                            height = ifelse(which == "row", 1/n, 1)))
      if (isTRUE(inner.border)) grid.rect()
      grid.draw(g)
      popViewport()
    })
  },
  var_import = list(panel_fun = panel_fun, which = which, outer.border = outer.border, inner.border = inner.border,
                    side = side, size = size, gap = gap, width = width, height = height),
  which = which,
  subsettable = TRUE,
  width = width,
  height = height)
}
