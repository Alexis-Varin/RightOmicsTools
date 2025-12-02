# Scatterplot of log-transformed counts, and average fitted GAM smoother curves along pseudotime

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
curveSmoothers(
  models,
  predictSmooth.df = NULL,
  genes = if (is.null(predictSmooth.df)) NA else unique(predictSmooth.df$gene),
  lineages,
  lineages.to.remove = NULL,
  conditions = NULL,
  facets = c("genes", "conditions"),
  nPoints = 100,
  branch.points = NULL,
  rescale = TRUE,
  rescale.range = c(0, 3),
  points.size = 3,
  points.alpha = 0.8,
  curves.width = 1,
  curves.alpha = 1,
  facets.scales = "fixed",
  facets.axes = "all",
  colors = NULL,
  axis.text.size = 9,
  axis.title.size = 11,
  facets.title.size = 9,
  legend.names = NULL,
  legend.text.size = 9,
  nrow = floor(sqrt(length(genes)))
)
```

## Arguments

- models:

  A SingleCellExperiment object containing the fitted GAM smoothers,
  computed using
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html), with or
  without `conditions` provided.

- predictSmooth.df:

  A `data.frame` object containing the average fitted GAM smoothers,
  computed using
  [`predictSmooth`](https://rdrr.io/pkg/tradeSeq/man/predictSmooth.html)
  with `nPoints`, `lineages` and `tidy` = `TRUE`. If `NULL`, average
  fitted GAM smoothers will be internally computed.

- genes:

  Character. The names of one or several genes to plot the
  log-transformed counts and the average fitted GAM smoothers from.

- lineages:

  Numeric. The indices of the lineages (for example, c(1, 5), 2:4 etc)
  to plot the log-transformed counts, the average fitted GAM smoothers
  and the pseudotime values from.

- lineages.to.remove:

  Numeric. The lineages to exclude from
  [`rescale`](https://scales.r-lib.org/reference/rescale.html). Useful
  if you want to remove the influence of other lineages in a comparison.

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

- facets:

  Character. The names of one or two metadata among 'genes', 'lineages'
  and/or 'conditions' to facet the plot. If one metadata is provided,
  [`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  will be used, and if two are provided,
  [`facet_grid`](https://ggplot2.tidyverse.org/reference/facet_grid.html)
  will be used, with the first metadata as rows and the second as
  columns. If the `models` object was computed using
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) without
  `conditions`, you cannot facet by 'conditions'. Please note that if
  the `models` object was computed using
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) with
  `conditions`, it is not possible to plot global lineages (without
  `conditions`), you would need to compute a new `models` object using
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) without
  `conditions`. This is due to a limitation in how
  [`predictSmooth`](https://rdrr.io/pkg/tradeSeq/man/predictSmooth.html)
  returns average fitted GAM smoothers (if the `models` object was
  computed with `conditions`, the function will always mean the fitted
  GAM smoothers for each lineage and each condition independently).

- nPoints:

  Numeric. (from
  [`predictSmooth`](https://rdrr.io/pkg/tradeSeq/man/predictSmooth.html)
  documentation) The number of points used to create the grid along the
  smoother for each lineage.

- branch.points:

  Numeric, Character or List. Branching points may be shown on the
  fitted GAM smoother curves, to partition the lineages and help
  visualize differences. Either one or more pseudotime values to plot
  the branching points at, and/or one or several values containing
  'knot' followed by a number (for example, 'knot2', 'knot4' etc), which
  correspond to the knots (`k`) input in
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) and which
  divide each lineage into segments; the function will extract the
  pseudotime values from each knot number and plot the branching points.
  Mixing any of the two options is possible (for example, c(4.3, 7,
  'knot4')). You may also provide a `list` containing any of the two
  options, and named after one or several identities of a single
  metadata among 'genes', 'lineages' or 'conditions' (for example,
  list('CD4' = c(2, 'knot5')) which will only plot branching points on
  CD4 curves). For `lineages`, the `list` names need to be provided as
  'Lineage' (with capital L) followed by the index (for example,
  list('Lineage1' = c(4.3, 'knot4'), 'Lineage2' = c('knot2', 7, 9)),
  which will only plot branching points on lineage 1 and 2 curves).
  Finally, you may also provide identities corresponding to several
  metadata at the same time by pasting them together with '\_' (for
  example, list('Lineage1_CD4' = c(1, 'knot1', 12)), which will only
  plot branching points on lineage 1 and CD4 curves).

- rescale:

  Logical. If `TRUE`, average fitted GAM smoothers will be adjusted
  using [`rescale`](https://scales.r-lib.org/reference/rescale.html)
  between the first numerical value of `rescale.range` (lowest value)
  and the second numerical value (highest value). This is different than
  [`scale`](https://rdrr.io/r/base/scale.html) as this doesn't compare
  values to any mean or standard deviation and is therefore not a
  Z-score, it only refits each value (independently for each gene) in
  order to visualize all `genes` in the same dimension regardless of
  their differences in fitted GAM smoothers.

- rescale.range:

  Numeric. The minimum and maximum values to resize the average fitted
  GAM smoothers and internally passed to
  [`rescale`](https://scales.r-lib.org/reference/rescale.html). These
  values are arbitrary and will not change the visualization, only the
  values in the legend. Ignored if `rescale` = `FALSE`.

- points.size:

  Numeric. The size of the log-transformed count points.

- points.alpha:

  Numeric. The transparency of the log-transformed count points.

- curves.width:

  Numeric. The width of the fitted GAM smoother curves.

- curves.alpha:

  Numeric. The transparency of the fitted GAM smoother curves.

- facets.scales:

  Character. (from
  [`facet_grid`](https://ggplot2.tidyverse.org/reference/facet_grid.html)
  documentation) Are scales shared across all facets (the default,
  "`fixed`"), or do they vary across rows ("`free_x`"), columns
  ("`free_y`"), or both rows and columns ("`free`")?

- facets.axes:

  Character. (from
  [`facet_grid`](https://ggplot2.tidyverse.org/reference/facet_grid.html)
  documentation) Determines which axes will be drawn. When "`margins`"
  (default), axes will be drawn at the exterior margins. "`all_x`" and
  "`all_y`" will draw the respective axes at the interior panels too,
  whereas "`all`" will draw all axes at all panels.

- colors:

  Character. The color names for each identity of the metadata excluded
  from `facets` among 'genes', 'lineages' and/or 'conditions' (for
  example, if `facets` = c('genes', 'conditions'), the `colors` will be
  for each of the `lineages` provided, or if `facets` = 'lineages', the
  `colors` will be for each of the `genes` times `conditions` provided).
  If `NULL`, uses Seurat's default colors.

- axis.text.size:

  Numeric. The font size of the pseudotime and log-transformed counts.

- axis.title.size:

  Numeric. The font size of the axis title.

- facets.title.size:

  Numeric. The font size of the facet titles.

- legend.names:

  Character. You may provide custom names for the colors displayed in
  the legend, instead of the default names (for example, if `facets` =
  c('genes', 'conditions'), the color names will be 'Lineage' followed
  by the index of the `lineages` provided, you may replace these with
  any other names).

- legend.text.size:

  Numeric. The font size of the legend text.

- nrow:

  Numeric. The number of rows to use for the facets. Ignored if two
  metadata are provided to `facets`.

## Value

A ggplot2 object.
