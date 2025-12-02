# Reduce the size of a Seurat object

As of Seurat v5 release,
[`DietSeurat`](https://satijalab.org/seurat/reference/DietSeurat.html)
does not remove data and scale.data layers, resulting in objects with
little to no reduction in size. This function was created with the
purpose to restore
[`DietSeurat`](https://satijalab.org/seurat/reference/DietSeurat.html)'s
proper behavior until the function is fixed by Seurat's dev team, as
well as to offer a few new options such as the ability to subset cells
or keep variable features.

## Usage

``` r
Right_DietSeurat(
  seurat_object,
  idents = NULL,
  cells = NULL,
  features = NULL,
  dimreducs = NULL,
  graphs = NULL,
  variable.features = FALSE,
  misc = TRUE,
  split.counts.layer = NULL,
  data.layer = FALSE,
  scale.layer = FALSE,
  SCTAssay = FALSE
)
```

## Arguments

- seurat_object:

  A Seurat object.

- idents:

  Character. The names of one or several metadata (for example,
  'orig.ident', 'seurat_clusters', etc) to keep in the diet object. If
  `NULL`, all metadata will be kept.

- cells:

  Character. The names of one or several cells to keep in the diet
  object. If `NULL`, all cells will be kept.

- features:

  Character. The names of one or several features to keep in the diet
  object. If `NULL`, all features will be kept.

- dimreducs:

  Character. The names of one or several reductions to keep in the diet
  object. If `NULL`, all reductions will be removed. Please note that if
  you subset features and remove all from which principal components
  (PC) were calculated, the 'PCA' reduction will not be kept, even if it
  is provided, as the feature.loadings slot will be empty. This is
  intended behavior.

- graphs:

  Character. The names of one or several graphs to keep in the diet
  object. If `NULL`, all graphs will be removed.

- variable.features:

  Logical. If `TRUE`, the variable features slot will be kept in the
  diet object. Please note that if you keep only a subset of features,
  the variable features will also be subset, and if you remove all
  features from which variable features were calculated, the variable
  features will not be kept even if `TRUE`. This is intended behavior.

- misc:

  Logical. If `TRUE`, the misc slot will be kept in the diet object.

- split.counts.layer:

  Character. The name of a metadata (for example, 'orig.ident',
  'seurat_clusters', etc) to split the 'counts' layer in the 'RNA' assay
  by if you need to in your downstream analysis. If `NULL`, the diet
  object will have a single 'counts' layer, even if the original object
  had split 'counts' layers.

- data.layer:

  Logical. If `TRUE`, the 'data' layer in the 'RNA' assay will be kept
  in the diet object if it is present. As with the 'counts' layer, if
  there are split 'data' layers, they will be joined into a single
  'data' layer unless `split.counts.layer` is provided.

- scale.layer:

  Logical. If `TRUE`, the 'scale.data' layer in the 'RNA' assay will be
  kept in the diet object if it is present.

- SCTAssay:

  Logical. If `TRUE`, the 'SCT' assay will be kept in the diet object if
  it is present.

## Value

A Seurat object, hopefully smaller, with class Assay5 'RNA' assay and
specified layers and slots.
