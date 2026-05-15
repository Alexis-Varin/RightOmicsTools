# ggplot annotation for ComplexHeatmap

This function creates several plotting regions (boxes) which corresponds
to each row/column of a
[`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html). While
[`anno_zoom`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_zoom.html) and
[`anno_link`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_link.html) are
primarily designed to add ggplot2 objects to slices (row_split,
column_split, km etc), adding them to each row/column is possible but
tedious, and the annotation cannot be decorated if sliced. This function
is designed to function similarly to them, with the same parameters, but
with a focus on adding the objects, especially ggplot2 objects, to each
row/column instead of slices, and to be decorated and sliced just like
any other annotation.

## Usage

``` r
anno_ggplot(
  panel_fun = function(index, nm = NULL) {
     grid.rect()
 },
  which = c("column", "row"),
  outer.border = TRUE,
  inner.border = FALSE,
  side = ifelse(which == "column", "top", "right"),
  size = NULL,
  gap = unit(1, "mm"),
  width = NULL,
  height = NULL
)
```

## Arguments

- panel_fun:

  Function. (from
  [`anno_zoom`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_zoom.html)
  documentation) A self-defined function that defines how to draw
  graphics in the box. The function must have a index argument which is
  the indices for the rows/columns that the box corresponds to. It can
  have second argument nm which is the "name" of the selected part in
  the heatmap.

- which:

  Character. (from
  [`anno_zoom`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_zoom.html)
  documentation) Whether it is a column annotation or a row annotation?

- outer.border:

  Logical. Whether to draw border around the annotation.

- inner.border:

  Logical. Whether to draw border around each box.

- side:

  Character. (from
  [`anno_zoom`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_zoom.html)
  documentation) Side of the boxes If it is a column annotation, valid
  values are "top" and "bottom"; If it is a row annotation, valid values
  are "left" and "right".

- size:

  Numeric or Unit. (from
  [`anno_zoom`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_zoom.html)
  documentation) The size of boxes. It can be pure numeric that they are
  treated as relative fractions of the total height/width of the
  heatmap. The value of size can also be absolute units.

- gap:

  Unit. (from
  [`anno_zoom`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_zoom.html)
  documentation) Gaps between boxes.

- width:

  Unit. (from
  [`anno_zoom`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_zoom.html)
  documentation) Width of the annotation. The value should be an
  absolute unit. Width is not allowed to be set for column annotation.

- height:

  Unit. (from
  [`anno_zoom`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_zoom.html)
  documentation) Height of the annotation. The value should be an
  absolute unit. Height is not allowed to be set for row annotation.

## Value

An annotation function which can be used in HeatmapAnnotation.
