# Heatmap of average fitted GAM smoothers along pseudotime

This function generates a heatmap to visualize gene fitted GAM
smoothers, computed using
[`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html), along
pseudotime.

## Usage

``` r
heatmapSmoothers(
  sds,
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
  ...
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

- predictSmooth.df:

  A `data.frame` object containing the average fitted GAM smoothers,
  computed using
  [`predictSmooth`](https://rdrr.io/pkg/tradeSeq/man/predictSmooth.html)
  with `nPoints`, `lineages` and `tidy` = `TRUE`. If `NULL`, average
  fitted GAM smoothers will be internally computed.

- genes:

  Character. The names of one or several genes to plot the average
  fitted GAM smoothers from.

- lineages:

  Numeric. The indices of the lineages (for example, c(1, 5), 2:4 etc)
  to plot the average fitted GAM smoothers and the pseudotime values
  from. Please note that for gene clustering, if you would like to use a
  specific lineage as reference, you need to set it as the first element
  (for example, if you have six lineages and you would like to cluster
  genes based on lineage 5, you would input `lineages` = c(5, 1:4, 6)).
  This is due to a limitation in how ComplexHeatmap clusters the rows of
  concatenated
  [`HeatmapList`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapList.html)
  objects (always based on the first element).

- lineages.to.remove:

  Numeric. The lineages to exclude from
  [`rescale`](https://scales.r-lib.org/reference/rescale.html). Useful
  if you want to remove the influence of other lineages in a comparison.

- conditions:

  Character. The names of one or several `conditions` identities to
  select. If `NULL`, all identities are used, and each unique condition
  will be plotted, for each of the `lineages` provided (for example, if
  you have four lineages and three conditions, twelve heatmaps will be
  plotted). Please note that for gene clustering, if you would like to
  use a specific condition as reference, you need to set it as the first
  element (for example, if you have three unique conditions (for
  example, 'control', 'IFNg' and 'TGFb') and you would like to cluster
  genes based on 'IFNg', you would input `conditions` = c('IFNg',
  'control', 'TGFb')). This is due to a limitation in how ComplexHeatmap
  clusters the rows of concatenated
  [`HeatmapList`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapList.html)
  objects (always based on the first element). Ignored if the `models`
  object was computed using
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

- clusters:

  Character or Factor. Either the name of a metadata present in the
  `sds` object (for example, 'annotations', 'seurat_clusters', etc) to
  plot the cell densities from, or the identities, as character or
  factor, of length equal to the number of cells in the `sds` object.

- nPoints:

  Numeric. (from
  [`predictSmooth`](https://rdrr.io/pkg/tradeSeq/man/predictSmooth.html)
  documentation) The number of points used to create the grid along the
  smoother for each lineage.

- branch.points:

  Numeric, Character or List. Branching points may be shown on the cell
  density plot and the heatmap, to partition the lineages and help
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
  metadata among 'lineages' or 'conditions' (for example, if
  `conditions` = c('Control', 'Treated'), list('Control' = c(2,
  'knot5')) which will only plot branching points on Control cell
  density plots and heatmaps). For `lineages`, the `list` names need to
  be provided as 'Lineage' (with capital L) followed by the index (for
  example, list('Lineage1' = c(4.3, 'knot4'), 'Lineage2' = c('knot2', 7,
  9)), which will only plot branching points on lineages 1 and 2 cell
  density plots and heatmaps). Finally, you may also provide identities
  corresponding to several metadata at the same time by pasting them
  together with '\_' (for example, list('Lineage1_Treated' = c(1,
  'knot1', 12)), which will only plot branching points on lineage 1 and
  Treated cell density plots and heatmaps).

- pseudotime.zoom:

  Numeric, Character or List. Oftentimes you may want to zoom into a
  specific region of a lineage to better visualize the differences.
  Either two pseudotime values to zoom into, and/or 'min' or 'max' to
  zoom respectively into the minimum or maximum pseudotime values of the
  lineage, and/or 'knot' followed by a number (for example, 'knot2',
  'knot4' etc), which correspond to the knots (`k`) input in
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) and which
  divide each lineage into segments; the function will extract the
  pseudotime values from each knot number and zoom into, and/or 'knot'
  followed by a number, a + or - and another number, which will set the
  zoom around knots by adding or subtracting from the pseudotime value
  of a knot (for example, c('knot3-5', 'knot4+2'), which corresponds to
  the pseudotime value of the third knot minus 5, and to the pseudotime
  value of the fourth knot plus 2). Mixing any two of the four options
  is possible (for example, c(5.5, 'knot4'), or c('knot3-1', 'max')).
  You may also provide a `list` containing any two of the four options,
  and named after one or several identities of a single metadata among
  'lineages' or 'conditions' (for example, if `conditions` =
  c('Control', 'Treated'), list('Control' = c(2, 'knot5')) which will
  zoom on Control cell density plots and heatmaps). For `lineages`, the
  `list` names need to be provided as 'Lineage' (with capital L)
  followed by the index (for example, list('Lineage1' = c('min',
  'knot2+3'), 'Lineage2' = c('knot2', 7)), which will only zoom on
  lineages 1 and 2 cell density plots and heatmaps). Finally, you may
  also provide identities corresponding to several metadata at the same
  time by pasting them together with '\_' (for example,
  list('Lineage1_Treated' = c('knot1+1', 'max')), which will only zoom
  on lineage 1 and Treated cell density plots and heatmaps).

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

- kmeans.repeats:

  Numeric. The number of runs to get a consensus K-means clustering.
  Ignored if `genes.kmeans` = 1.

- cluster.genes:

  Logical or Function. If `TRUE`, the function will cluster the `genes`.
  You may also pass an `hclust` or `dendrogram` object which contains
  clustering.

- genes.kmeans:

  Numeric. The number of slices to use for gene K-means clustering.

- genes.kmeans.numbers.size:

  Numeric. The font size of the gene K-means slice numbers. Set to 0 to
  remove them.

- outer.border:

  Logical. If `TRUE`, the function will display an outer border around
  each heatmap.

- pseudotime.type:

  Character. Determines the `pseudotime.colors` range scale of each
  lineage: either 'independent', where each lineage's
  `pseudotime.colors` range is independent of the others (i.e., minimum
  and maximum values are identical for each lineage, regardless of
  trajectory length differences), or 'relative', where each lineage's
  `pseudotime.colors` range is scaled relative to the highest pseudotime
  value (i.e., the longest trajectory, reflecting differences in
  trajectory lengths).

- density.type:

  Character. Determines the cell density scale (height) of each lineage:
  either 'independent,' where each lineage's cell density is independent
  of the others (i.e., the maximum density height is identical for each
  lineage, regardless of differences in cell counts contributing to a
  trajectory), or 'relative,' where each lineage's density scale is
  relative to the highest cell density (i.e., the highest cell counts,
  reflecting differences in the number of cells contributing to a
  trajectory).

- show.pseudotime:

  Logical. If `TRUE`, the pseudotime range will be shown above the
  heatmap.

- show.density:

  Logical. If `TRUE`, the cell density plot will be shown.

- show.density.pseudotime.values:

  Logical. If `TRUE`, the pseudotime values will be shown on the cell
  density plot.

- data.colors:

  Character. Either two color names, corresponding to the lowest and
  highest values in the average fitted GAM smoothers and internally
  passed to
  [`colorRamp2`](https://rdrr.io/pkg/colorRamp2/man/colorRamp2.html), or
  the name of a palette and internally passed to `hcl_palette` in
  [`colorRamp2`](https://rdrr.io/pkg/colorRamp2/man/colorRamp2.html)
  (such as 'Inferno', 'Berlin', 'Viridis' etc, check
  [`hcl.pals`](https://rdrr.io/r/grDevices/palettes.html) for all
  palettes available).

- pseudotime.colors:

  Character. Either the name of a palette and internally passed to
  `hcl_palette` in
  [`colorRamp2`](https://rdrr.io/pkg/colorRamp2/man/colorRamp2.html)
  (such as 'Inferno', 'Berlin', 'Viridis' etc, check
  [`hcl.pals`](https://rdrr.io/r/grDevices/palettes.html) for all
  palettes available), or two or more color names, corresponding to the
  lowest and highest pseudotime values, and additional color names
  beyond two will either be spaced at regular intervals if unnamed (for
  example, c('blue', 'green', 'yellow', 'red') will be spaced at 0%,
  33%, 66% and 100%) or if named with a pseudotime value (for example,
  c('orange', "12" = 'white', 'royalblue')) or if named with 'knot'
  followed by a number (for example, 'knot2', 'knot4' etc), which
  correspond to the knots (`k`) input in
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) and which
  divide each lineage into segments; the function will extract the
  pseudotime values from each knot number (for example, c('firebrick',
  'knot3' = 'lightgrey', 'gold')). Note that the first and last color
  names do not need to be named as they always correspond to the lowest
  and highest pseudotime values.

- clusters.colors:

  Character. The color names for each identity in `clusters`. If `NULL`,
  uses Seurat's default colors. Ignored if `show.density` = `FALSE`.

- genes.names.size:

  Numeric. The font size of the gene names. Set to 0 to remove them.

- genes.names.style:

  Character. The font face of the gene names. The [Gene
  nomenclature](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7494048/)
  used by almost all scientific journals require that gene names are
  italicized, therefore the parameter is by default set to 'italic'. Use
  'plain' to revert back to regular font face.

- show.density.legend:

  Logical. If `TRUE`, the cell density plot legend will be shown.

- density.legend.name:

  Character. The name of the cell density plot legend.

- heatmap.names:

  Character. You may provide custom names for each heatmap, instead of
  'Lineage' followed by the index of the `lineages` and `conditions`
  provided.

- heatmap.width:

  Numeric. The width of each heatmap.

- heatmap.height:

  Numeric. The height of each heatmap.

- density.height:

  Numeric. The height of the cell density plot.

- raster:

  Logical. (from
  [`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  documentation) If `TRUE`, the function will render the heatmap body as
  a raster image. It helps to reduce file size when the matrix is huge.

- return.grob:

  Logical. If `TRUE`, the function will return a `grob` (`gTree` object,
  printable using [`grid.draw`](https://rdrr.io/r/grid/grid.draw.html)).
  This allows to capture the cell density plot (which is a ggplot2
  object) and the heatmap (a
  [`HeatmapList`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapList.html)
  object) as a single object which can be further customized using grid
  functions or included with other plots in complex layouts using
  [`plot_grid`](https://wilkelab.org/cowplot/reference/plot_grid.html).
  If `FALSE`, the function will return a
  [`HeatmapList`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapList.html)
  object, which can then be used to extract various parameters, such as
  the row clustering values using
  [`row_order`](https://rdrr.io/pkg/ComplexHeatmap/man/row_order-dispatch.html).
  Please note that this only affects the returned object when it is
  assigned to a variable, such as `> htmp <- heatmapSmoothers(...)`, and
  not when the function is called without assignment, such as
  `> heatmapSmoothers(...)`.

- ...:

  Additional arguments to be passed to
  [`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html), such
  as `show_parent_dend_line`, `clustering_method_rows`, etc, accepts any
  parameter that wasn't already internally passed to
  [`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html) (for
  example, `outer.border` sets the `border` parameter of
  [`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html), so
  you will get an error if you try to pass the `border` parameter in
  `heatmapSmoothers`).

## Value

A `grob` (`gTree` object, printable using
[`grid.draw`](https://rdrr.io/r/grid/grid.draw.html)), or a
[`HeatmapList`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapList.html)
object.
