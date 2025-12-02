# Data loader for vignettes

This function downloads and prepares data used for vignettes. The main
advantage of using a data loader is that nothing is stored in the
package, therefore it remains lightweight even if you use large objects
in your tutorials. Anyone interested in adding a dataset to this
function for use in their own package's vignettes (provided the raw data
can be downloaded from a reputable source, [NCBI's GEO FTP
site](https://ftp.ncbi.nlm.nih.gov/geo/datasets/) for example) can open
an [issue on
GitHub](https://github.com/Alexis-Varin/RightOmicsTools/issues).

## Usage

``` r
Right_Data(dataset = NULL, list.datasets = FALSE)
```

## Arguments

- dataset:

  Character. The name of the dataset to load.

- list.datasets:

  Logical. If `TRUE`, prints the names and descriptions of available
  datasets.

## Value

An object of various classes, depending on the dataset, or a message
containing the names and descriptions of available datasets.

## Examples

``` r
Right_Data(list.datasets = TRUE)
#> Available datasets:
#> "pbmc3k"   a Seurat v4 object of 2,700 cells by 13,714 genes
#>            PBMCs of a healthy donor from 10XGenomics
#> 
#> "monolps"   a Seurat v5 object of 3,693 cells by 23,798 genes
#>             Monocytes stimulated or not with LPS from GSE226488

pbmc3k = Right_Data("pbmc3k")
pbmc3k
#> An object of class Seurat 
#> 13714 features across 2700 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  2 layers present: counts, data
#>  1 dimensional reduction calculated: umap
```
