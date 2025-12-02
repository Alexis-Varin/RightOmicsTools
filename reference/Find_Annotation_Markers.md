# Get the top markers for fast annotation

This function is a wrapper around
[`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
that allows for parallelization and filtering of mitochondrial,
ribosomal and non-coding RNA features in human, as well as filtering of
pseudogenes in mouse. It may also directly return the top X markers for
each identity.

## Usage

``` r
Find_Annotation_Markers(
  seurat_object,
  ident.1 = NULL,
  ident.2 = NULL,
  min.pct = 0.25,
  top.markers = 5,
  unique.markers = TRUE,
  filter.mito = TRUE,
  filter.ribo = TRUE,
  filter.ncRNA = TRUE,
  species = "human",
  parallelized = FALSE,
  BPPARAM = NULL,
  name.markers = FALSE,
  output.df = FALSE,
  output.list = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- seurat_object:

  A Seurat object.

- ident.1:

  Character. (from
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  documentation) Identity class to define markers for; pass an object of
  class `phylo` or 'clustertree' to find markers for a node in a cluster
  tree; passing 'clustertree' requires
  [`BuildClusterTree`](https://satijalab.org/seurat/reference/BuildClusterTree.html)
  to have been run. Leave `NULL` to find markers for all clusters.

- ident.2:

  Character. (from
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  documentation) A second identity class for comparison; if `NULL`, use
  all other cells for comparison; if an object of class `phylo` or
  'clustertree' is passed to `ident.1`, must pass a node to find markers
  for.

- min.pct:

  Numeric. (from
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  documentation) Only test features that are detected in a minimum
  fraction of min.pct cells in either of the two populations. Meant to
  speed up the function by not testing features that are very
  infrequently expressed.

- top.markers:

  Numeric. The number of top markers to return. If `Inf`, all markers
  will be returned.

- unique.markers:

  Logical. If `TRUE`, unique markers will be returned for each identity
  in order to prevent features repeated multiple times.

- filter.mito:

  Logical. If `TRUE`, mitochondrial features will be filtered out.

- filter.ribo:

  Logical. If `TRUE`, ribosomal features will be filtered out.

- filter.ncRNA:

  Logical. If `TRUE`, non-coding RNA features will be filtered out.

- species:

  Character. The species name to filter out non-coding RNA features. If
  'human', a dataset named
  [ncRNA_human](https://alexis-varin.github.io/RightOmicsTools/reference/ncRNA_human.html)
  built from [genenames
  database](https://www.genenames.org/data/genegroup/#!/group/475) will
  be used as reference. If 'mouse', only pseudogenes will be filtered
  out based on a dataset named
  [pseudogenes_mouse](https://alexis-varin.github.io/RightOmicsTools/reference/pseudogenes_mouse.html)
  and built from [dreamBase2
  database](https://rna.sysu.edu.cn/dreamBase2/scrna.php?SClade=mammal&SOrganism=mm10&SDataId=0&SProteinID=0).
  These datasets are loaded with RightOmicsTools and may be checked for
  more information.

- parallelized:

  Logical. If `TRUE`,
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  will be parallelized using BiocParallel. Please note that
  parallelization is complex and depends on your operating system
  (Windows users might not see a gain or might even experience a
  slowdown).

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object to be used for parallelization. If `NULL` and `parallelized` =
  `TRUE`, the function will use a
  [`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  object configured to use a single worker (core) and is therefore not
  parallelized, in order to prevent accidental use of large computation
  resources. Ignored if `parallelized` = `FALSE`.

- name.markers:

  Logical. If `TRUE`, each marker will be named with its associated
  identity. Ignored if `output.df` = `TRUE`.

- output.df:

  Logical. If `TRUE`, a `data.frame` object containing the
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  results, ordered by decreasing log fold change will be returned. If
  `FALSE`, feature names will be returned.

- output.list:

  Logical. If `TRUE`, and if `output.df` = `TRUE`, returns a `list` of
  `data.frame` objects containing
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  results for each identity, or if `output.df` = `FALSE`, returns a
  `list` of feature names for each identity.

- verbose:

  Logical. If `FALSE`, does not print progress messages and output, but
  warnings and errors will still be printed.

- ...:

  Additional arguments to be passed to
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html),
  such as `test.use`, or passed to other methods and to specific DE
  methods.

## Value

A `data.frame` object or a `list` of `data.frame` objects containing the
[`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
results, or a `list` of features names, or feature names.

## Examples

``` r
# Prepare data
pbmc3k <- Right_Data("pbmc3k")

# Example 1: default parameters and origin of markers
pbmc3k.markers <- Find_Annotation_Markers(pbmc3k,
                                          name.markers = TRUE)
#> Finding markers for cluster Naive CD4 T against all other clusters
#> Finding markers for cluster CD14+ Mono against all other clusters
#> Finding markers for cluster Memory CD4 T against all other clusters
#> Finding markers for cluster B against all other clusters
#> Finding markers for cluster CD8 T against all other clusters
#> Finding markers for cluster FCGR3A+ Mono against all other clusters
#> Finding markers for cluster NK against all other clusters
#> Finding markers for cluster DC against all other clusters
#> Finding markers for cluster Platelets against all other clusters
pbmc3k.markers
#>     Naive CD4 T     Naive CD4 T     Naive CD4 T     Naive CD4 T     Naive CD4 T 
#>          "CCR7"          "LEF1"           "MAL"       "PIK3IP1"          "TCF7" 
#>      CD14+ Mono      CD14+ Mono      CD14+ Mono      CD14+ Mono      CD14+ Mono 
#>         "FOLR3"       "S100A12"        "S100A8"        "S100A9"          "CD14" 
#>    Memory CD4 T    Memory CD4 T    Memory CD4 T    Memory CD4 T    Memory CD4 T 
#>        "CD40LG"          "AQP3"         "SUSD3"           "CD2"         "TRAT1" 
#>               B               B               B               B               B 
#>        "VPREB3"         "CD79A"         "FCRLA"         "TCL1A"         "FCER2" 
#>           CD8 T           CD8 T           CD8 T           CD8 T           CD8 T 
#>          "GZMK"          "GZMH"          "CD8A"          "CCL5"         "KLRG1" 
#>    FCGR3A+ Mono    FCGR3A+ Mono    FCGR3A+ Mono    FCGR3A+ Mono    FCGR3A+ Mono 
#>           "CKB"        "CDKN1C"        "MS4A4A"          "HES4"         "BATF3" 
#>              NK              NK              NK              NK              NK 
#>        "AKR1C3"          "GZMB"        "SH2D1B"         "SPON2"        "FGFBP2" 
#>              DC              DC              DC              DC              DC 
#>      "SERPINF1"        "FCER1A"         "CLIC2"       "CLEC10A"          "ENHO" 
#>       Platelets       Platelets       Platelets       Platelets       Platelets 
#>        "LY6G6F" "RP11-879F14.2"         "CLDN5"           "GP9"        "ITGA2B" 

# Example 2: parallelized FindAllMarkers
BPPARAM <- BiocParallel::registered()[[1]]
if (BPPARAM$workers > 4) BPPARAM$workers <- 4
pbmc3k.markers <- Find_Annotation_Markers(pbmc3k,
                                          min.pct = 0.01,
                                          top.markers = Inf,
                                          unique.markers = FALSE,
                                          filter.mito = FALSE,
                                          filter.ribo = FALSE,
                                          filter.ncRNA = FALSE,
                                          parallelized = TRUE,
                                          BPPARAM = BPPARAM,
                                          output.df = TRUE)
#> Finding markers for cluster FCGR3A+ Mono against all other clustersFinding markers for cluster NK against all other clustersFinding markers for cluster DC against all other clustersFinding markers for cluster Platelets against all other clusters
#> Finding markers for cluster Naive CD4 T against all other clustersFinding markers for cluster CD14+ Mono against all other clustersFinding markers for cluster Memory CD4 T against all other clustersFinding markers for cluster B against all other clustersFinding markers for cluster CD8 T against all other clusters
head(pbmc3k.markers)
#>                p_val avg_log2FC pct.1 pct.2    p_val_adj     cluster feature
#> GTSCR1  1.720280e-08   7.163733 0.016 0.000 2.359193e-04 Naive CD4 T  GTSCR1
#> REG4    4.484775e-10   5.902153 0.022 0.000 6.150421e-06 Naive CD4 T    REG4
#> C2orf40 1.024936e-10   5.644114 0.023 0.000 1.405597e-06 Naive CD4 T C2orf40
#> MMP28   1.729283e-07   4.380462 0.016 0.000 2.371539e-03 Naive CD4 T   MMP28
#> NOG     5.025979e-10   4.148065 0.032 0.003 6.892627e-06 Naive CD4 T     NOG
#> FAM153A 8.848417e-08   3.880304 0.020 0.001 1.213472e-03 Naive CD4 T FAM153A
```
