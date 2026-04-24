# Scatterplot of the dimensionality reduction and fitted GAM smoothers for each cell

This function generates a scatterplot of the log-transformed counts, and
displays on top the average fitted GAM smoother curves, computed using
[`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html), along
pseudotime, and is a reworked version of
[`plotSmoothers`](https://rdrr.io/pkg/tradeSeq/man/plotSmoothers.html),
it allows for the plotting of multiple genes simultaneously, as well as
selected lineages and conditions. Additionally, the scatterplot is
rasterized to reduce the size of the output file and branching points
may be displayed on each curve. Finally, the log-transformed counts and
average fitted GAM smoother curves are rescaled to the same range
between genes to allow for a better comparison.

## Usage

``` r
dimSmoothers(
  sds,
  models = NULL,
  sds.conditions = NULL,
  genes = NULL,
  lineages = NULL,
  clusters = NULL,
  conditions = NULL,
  knots = NULL,
  show.knots.labels = FALSE,
  rescale = FALSE,
  rescale.range = c(0, 3),
  pseudotime.type = "relative",
  points.size = 1,
  lineages.width = 1,
  lineages.arrow = arrow(),
  knots.size = lineages.width * 3,
  cells.colors = if (is.character(genes)) "turbo" else if (any(clusters == "pseudotime"))
    "viridis" else NULL,
  na.color = "lightgrey",
  show.na = TRUE,
  lineages.colors = "grey50",
  knots.colors = lineages.colors,
  knots.labels.text.colors = "white",
  knots.labels.fill.colors = lineages.colors,
  axis.text.size = 9,
  axis.title.size = 11,
  plot.title.size = 16,
  lineages.legend.names = NULL,
  legend.text.size = 9,
  raster = FALSE
)
```

## Arguments

- sds:

  A SingleCellExperiment object containing the pseudotime values of one
  or several lineages, computed using
  [`slingshot`](https://rdrr.io/pkg/slingshot/man/slingshot.html)
  (usually, the input object to
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html)).

- models:

  A SingleCellExperiment object containing the fitted GAM smoothers,
  computed using
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html), with or
  without `conditions` provided.

- sds.conditions:

  A list of SingleCellExperiment objects obtained using
  `slingshot_conditions`. The conditions used to compute the pseudotime
  values must be the same as those used to fit the GAM smoothers in the
  `models` object, if provided.

- genes:

  Character. The names of one or several genes to plot the fitted GAM
  smoothers for each cell from. If the `models` object is not provided,
  the function will plot the log-transformed counts.

- lineages:

  Numeric. The indices of the lineages (for example, c(1, 5), 2:4 etc)
  to plot the pseudotime values from. Set to `NA` to remove all
  lineages.

- clusters:

  Character or Factor. Either the name of a metadata present in the
  `sds` object (for example, 'annotations', 'seurat_clusters', etc) to
  color the cells, or the identities, as character or factor, of length
  equal to the number of cells in the `sds` object. Set to 'pseudotime'
  to color the cells by their pseudotime values.

- conditions:

  Character. The names of one or several `conditions` identities to
  select. If `NULL`, all identities are used, and each unique condition
  will be plotted, for each of the `lineages` provided (for example, if
  you have four lineages and three conditions, twelve curves will be
  plotted). Ignored if the `models` object was computed using
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) without
  `conditions`. Please note that if the `models` object was computed
  using [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) with
  `conditions`, it is not possible to plot global lineages (without
  `conditions`), you would need to compute a new `models` object using
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) without
  `conditions`. This is due to a limitation in how
  [`predictSmooth`](https://rdrr.io/pkg/tradeSeq/man/predictSmooth.html)
  returns average fitted GAM smoothers (if the `models` object was
  computed with `conditions`, the function will always mean the fitted
  GAM smoothers for each lineage and each condition independently).

- knots:

  Numeric. The indices of the knots to display on the `lineages`. Set to
  `NA` to remove all knots. Ignored if the `models` object is not
  provided.

- show.knots.labels:

  Logical. If `TRUE`, the knot numbers will be displayed using
  [`geom_label_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Ignored if the `models` object is not provided.

- rescale:

  Logical. If `TRUE`, fitted GAM smoothers will be adjusted using
  [`rescale`](https://scales.r-lib.org/reference/rescale.html) between
  the first numerical value of `rescale.range` (lowest value) and the
  second numerical value (highest value). This is different than
  [`scale`](https://rdrr.io/r/base/scale.html) as this doesn't compare
  values to any mean or standard deviation and is therefore not a
  Z-score, it only refits each value (independently for each gene) in
  order to visualize all `genes` in the same dimension regardless of
  their differences in fitted GAM smoothers.

- rescale.range:

  Numeric. The minimum and maximum values to resize the fitted GAM
  smoothers and internally passed to
  [`rescale`](https://scales.r-lib.org/reference/rescale.html). These
  values are arbitrary and will not change the visualization, only the
  values in the legend. Ignored if `rescale` = `FALSE`.

- pseudotime.type:

  Character. Determines the `cells.colors` range scale when `clusters` =
  'pseudotime' and not all `lineages` are displayed: either 'absolute',
  where the `cells.colors` range is based on the highest pseudotime
  value among the displayed `lineages`, or 'relative', where the
  `cells.colors` range is based on the highest pseudotime value of all
  the `lineages`, including the non-displayed ones.

- points.size:

  Numeric. The size of the cells.

- lineages.width:

  Numeric. The width of the lineage curves.

- lineages.arrow:

  An [`arrow`](https://rdrr.io/r/grid/arrow.html) object added at the
  end of the lineage curves. Use the function's parameters to modify its
  appearance. Set to `NULL` to remove the arrow.

- knots.size:

  Numeric. The size of the knots on the lineage curves. Ignored if the
  `models` object is not provided.

- cells.colors:

  Character. If `genes` is provided, either two color names to use for
  the cells, corresponding to the lowest and highest values in the
  fitted GAM smoothers and internally passed to
  [`scale_color_gradient`](https://ggplot2.tidyverse.org/reference/scale_gradient.html),
  or the name of a palette and internally passed to
  [`scale_color_viridis_c`](https://ggplot2.tidyverse.org/reference/scale_viridis.html)
  (such as 'turbo', 'inferno', 'magma', 'viridis' etc, check the
  function for all palettes available). If `clusters` is provided, the
  color names for each identity in `clusters`. If `NULL`, uses Seurat's
  default colors. If `clusters` = 'pseudotime', the name of a palette
  and internally passed to
  [`scale_color_viridis_c`](https://ggplot2.tidyverse.org/reference/scale_viridis.html).
  If neither `genes` nor `clusters` are provided, the name of a color to
  use for all cells. If `NULL`, defaults to 'lightgrey'.

- na.color:

  Character. The name of a color to use for cells not participating in
  the `lineages` (for example, if you have four lineages and you only
  plot lineages 1 and 3, the cells participating in lineages 2 and 4
  will be colored with `na.color`). Set to 'transparent' to not display
  these cells.

- show.na:

  Logical. If `FALSE`, ignores `na.color` and colors all cells with
  `cells.colors`, even those not participating in the `lineages`.

- lineages.colors:

  Character. either a single color name to use for all lineages, or the
  color names for each of the `lineages` displayed, or for the number of
  `lineages` times the number of `conditions` if the `models` object was
  computed with `conditions`.

- knots.colors:

  Character. either a single color name to use for all knots, or a
  vector of color names of length equal to the number of `lineages`.
  Ignored if the `models` object is not provided.

- knots.labels.text.colors:

  Character. either a single color name to use for all knot labels, or a
  vector of color names of length equal to the number of `lineages`.
  Ignored if the `models` object is not provided or if
  `show.knots.labels` = `FALSE`.

- knots.labels.fill.colors:

  Character. either a single color name to use for all knot label
  backgrounds, or a vector of color names of length equal to the number
  of `lineages`. Ignored if the `models` object is not provided or if
  `show.knots.labels` = `FALSE`.

- axis.text.size:

  Numeric. The font size of the dimensionality reduction values.

- axis.title.size:

  Numeric. The font size of the axis title.

- plot.title.size:

  Numeric. The font size of the plot title.

- lineages.legend.names:

  Character. You may provide custom names for the lineages displayed in
  the legend. Ignored if `lineages.colors` is a single color.

- legend.text.size:

  Numeric. The font size of the legend text.

- raster:

  Logical. If `TRUE`, the cells will be rasterized using
  [`rasterize`](https://rdrr.io/pkg/ggrastr/man/rasterize.html) to
  reduce the size of the output file. The points may appear slightly
  blurry.

## Value

A ggplot2 object, or a `patchwork` object containing ggplot2 objects.
