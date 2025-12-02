# Gene signatures from GSEA-MSigDB pathways

Following gene set enrichment analysis (GSEA), one often would like to
explore the expression of genes comprised in enriched pathways. For this
purpose, this function builds a pathway database from
[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) and creates signatures
(module scores calculated from UCell or Seurat's respective functions)
from features found in a Seurat object and extracted from supplied
pathways. It also returns the feature names which can be used to
visualize their individual expression, using for example
[`DotPlot_Heatmap`](https://alexis-varin.github.io/RightOmicsTools/reference/DotPlot_Heatmap.md).

## Usage

``` r
GSEA_Signatures(
  seurat_object,
  assay = "RNA",
  layer = "data",
  species = "Homo sapiens",
  category = NULL,
  subcategory = NULL,
  pathways,
  min.features = 2,
  signatures.names = "name",
  method = "UCell",
  only.features = FALSE,
  fail.safe = 10,
  verbose = TRUE,
  ...
)
```

## Arguments

- seurat_object:

  A Seurat object.

- assay:

  Character. The name of an assay containing the `layer` with the
  expression matrix. If the `seurat_object` contains multiple 'RNA'
  assays, you may specify which one to use (for example, 'RNA2' if you
  have created a second 'RNA' assay you named 'RNA2'. See [Seurat v5
  vignettes](https://satijalab.org/seurat/articles/seurat5_essential_commands.html#create-seurat-or-assay-objects)
  for more information). You may also use another assay, such as 'SCT',
  to pull feature expression from.

- layer:

  Character. The name of a layer (formerly known as slot) which stores
  the expression matrix. If the `seurat_object` contains split layers,
  the function will always join them before searching features and
  adding the signatures.

- species:

  Character. The species name to be internally passed to
  [`msigdbr`](https://igordot.github.io/msigdbr/reference/msigdbr.html)
  to build the pathway database. Use
  msigdbr::[`msigdbr_species`](https://igordot.github.io/msigdbr/reference/msigdbr_species.html)
  for the names of available species.

- category:

  Character. The names of one or several categories to be internally
  passed to
  [`msigdbr`](https://igordot.github.io/msigdbr/reference/msigdbr.html)
  to build the pathway database. Use
  msigdbr::[`msigdbr_collections`](https://igordot.github.io/msigdbr/reference/msigdbr_collections.html)
  for the names of available categories (gs_cat column). If `NULL`, all
  categories will be used.

- subcategory:

  Character. The names of one or several subcategories to be internally
  passed to
  [`msigdbr`](https://igordot.github.io/msigdbr/reference/msigdbr.html)
  to build the pathway database. Use
  msigdbr::[`msigdbr_collections`](https://igordot.github.io/msigdbr/reference/msigdbr_collections.html)
  for the names of available subcategories (gs_subcat column). If
  `NULL`, all subcategories will be used.

- pathways:

  Character. The names of one or several pathways to be searched in the
  pathway database and added as signatures. You may provide either a
  pathway id (for example, 'GO:0006574') or a name matching the pattern
  found in the gs_name column (uppercase letters and underscores between
  words). Please note that you may also provide a partial match (for
  example, 'TYPE_I_INTERFERON') and the function will find all pathways
  containing this partial pattern. Beware that this may result in a
  large number of pathways to be added as signatures (using
  `only.features` = `TRUE` is highly recommended) but is very handy to
  explore all pathways of interest in a particular biological process.

- min.features:

  Numeric. The minimum number of features present in the `seurat_object`
  for a pathway to be added as a signature.

- signatures.names:

  Character. Either 'id', which will add the ids of the `pathways` as
  signature names (for example, 'GO:0004657', 'hsa05200' etc), or
  'name', which will add the names of the `pathways` as signature names
  (for example, 'GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY'). You may also
  provide custom names to be used as signature names, whose length must
  match the length of `pathways` supplied. If multiple results are found
  for a pathway, the function will append a number to the corresponding
  custom signature name for each result.

- method:

  Character. The method used to calculate the module scores, either
  'UCell' or 'Seurat'.

- only.features:

  Logical. If `TRUE`, the function will not add any signature to the
  `seurat_object` and will only return the feature names from the
  `pathways` found in the `seurat_object` and the feature names present
  in the `seurat_object`.

- fail.safe:

  Numeric. The maximum number of signatures the function will attempt to
  add to the `seurat_object`. If the number of signatures found is
  higher than this number, the function will not add any signature, and
  will instead return the `seurat_object` as well as the feature names
  from the `pathways` found in the `seurat_object` and the feature names
  present in the `seurat_object`. This prevents the function from adding
  a large number of signatures to the `seurat_object` by mistake.

- verbose:

  Logical. If `FALSE`, does not print progress messages and output, but
  warnings and errors will still be printed.

- ...:

  Additional arguments to be passed to
  [`AddModuleScore_UCell`](https://rdrr.io/pkg/UCell/man/AddModuleScore_UCell.html)
  or
  [`AddModuleScore`](https://satijalab.org/seurat/reference/AddModuleScore.html),
  such as `nbin` or `maxRank`.

## Value

A `list` containing the `seurat_object` with added signatures, all
feature names from the `pathways` found in the `seurat_object`, the
feature names present in the `seurat_object` and the signature names. If
`only.features` = `TRUE`, the function will instead return a `list`
containing the feature names from the `pathways` found in the
`seurat_object` and the feature names present in the `seurat_object`.
