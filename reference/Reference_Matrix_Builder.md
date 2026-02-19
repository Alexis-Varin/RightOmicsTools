# Builder of CIBERSORTx Reference Matrix

This function builds a Reference Matrix for
[CIBERSORTx](https://cibersortx.stanford.edu/index.php) from a Seurat
object.

## Usage

``` r
Reference_Matrix_Builder(
  seurat_object,
  assay = "RNA",
  layer = "counts",
  ident.1 = NULL,
  ident.2 = NULL,
  double.ident = TRUE,
  double.ident.sep = "_",
  reverse.double.ident = FALSE,
  subset.ident.1 = NULL,
  subset.ident.1.invert = FALSE,
  subset.ident.2 = NULL,
  subset.ident.2.invert = FALSE,
  downsample.object.first = NULL,
  downsample.object.last = NULL,
  downsample.threshold = 2,
  downsample.cluster = NULL,
  automatic.downsample = FALSE,
  check.size = FALSE,
  max.matrix.size = 500,
  cell.barcodes = FALSE,
  output.cell.barcodes = FALSE,
  file.name = "Reference_Matrix",
  file.format = "txt",
  file.sep = "\t",
  write.path = NULL,
  write = TRUE,
  verbose = TRUE
)
```

## Arguments

- seurat_object:

  A Seurat object.

- assay:

  Character. The name of an assay containing the `layer` with the
  expression matrix. If the `seurat_object` contains multiple 'RNA'
  assays, you may specify which one to use (for example 'RNA2' if you
  have created a second 'RNA' assay you named 'RNA2'. See [Seurat v5
  vignettes](https://satijalab.org/seurat/articles/seurat5_essential_commands.html#create-seurat-or-assay-objects)
  for more information). You may also use another assay, such as 'SCT',
  if you want to extract the expression matrix for projects other than
  CIBERSORTx.

- layer:

  Character. The name of a layer (formerly known as slot) which stores
  the expression matrix. Must be 'counts' layer for CIBERSORTx. You may
  also specify another layer such as 'data' (log-normalized counts) or
  any other present in the specified `assay` if you want to extract the
  expression matrix for projects other than CIBERSORTx. If you have
  split layers the function will always join them before extracting the
  expression matrix. The raw counts matrix will be normalized to Counts
  Per Million (CPM) for CIBERSORTx.

- ident.1:

  Character. The name of a metadata (for example, 'orig.ident',
  'seurat_clusters', etc) from which to extract identities as column
  names in the Reference Matrix. If `NULL`, uses the active.ident
  metadata.

- ident.2:

  Character. The name of a metadata from which to split `ident.1`
  identities by (for example if `ident.1` is 'seurat_clusters' and
  `ident.2` is 'orig.ident', identity 0 will be separated into several
  0_identities, each corresponding to an orig.ident: 0_sample1,
  0_sample2 etc).

- double.ident:

  Logical. If `TRUE`, column names of the Reference Matrix will be set
  with dual identities (ident.1_ident.2). If `FALSE`, column names will
  be set as `ident.1` identities. Ignored if `ident.2` = `NULL`.

- double.ident.sep:

  Character. The separator to use between the identity names of
  `ident.1` and `ident.2`. Ignored if `ident.2` = `NULL` or
  `double.ident` = `FALSE`.

- reverse.double.ident:

  Logical. If `TRUE`, the identity names of `ident.1` and `ident.2` will
  be reversed (ident.2_ident.1). Ignored if `ident.2` = `NULL` or
  `double.ident` = `FALSE`.

- subset.ident.1:

  Character. The names of one or several `ident.1` identities to select.
  If `NULL`, all identities are used.

- subset.ident.1.invert:

  Logical. If `TRUE`, inverts the selection from `subset.ident.1` (for
  example, if `subset.ident.1` = c('0','1','2') and
  `subset.ident.1.invert` = `TRUE`, all identities except 0, 1 and 2
  will be kept). Ignored if `subset.ident.1` = `NULL`.

- subset.ident.2:

  Character. The names of one or several `ident.2` identities to select.
  If `NULL`, all identities are used. Ignored if `ident.2` = `NULL`.

- subset.ident.2.invert:

  Logical. If `TRUE`, inverts the selection from `subset.ident.2` (for
  example, if `subset.ident.2` = c('patient1','patient4') and
  `subset.ident.2.invert` = `TRUE`, all identities except patient 1 and
  patient 4 will be kept). Ignored if `subset.ident.2` = `NULL` or
  `ident.2` = `NULL`.

- downsample.object.first:

  Numeric. The number of cells to downsample from the entire
  `seurat_object` (not from each identity) before subsetting with
  `subset.ident.1` and/or `subset.ident.2`. Please note that the
  relative proportion of each identity will be kept (for example if an
  identity is 10% of total cells, it will remain 10% after
  downsampling). If `NULL`, all cells are used.

- downsample.object.last:

  Numeric. The number of cells to downsample from the entire
  `seurat_object` after subsetting with `subset.ident.1` and/or
  `subset.ident.2`. Please note that the relative proportion of each
  identity will be kept (for example if an identity is 10% of total
  cells, it will remain 10% after downsampling). If `NULL`, all cells
  are used.

- downsample.threshold:

  Numeric. The minimum number of cells (1 or more) or the minimum
  proportion of cells (between 0 and 1) to take in each cluster when
  downsampling. Ignored for `downsample.cluster`.

- downsample.cluster:

  Numeric. The number of cells to downsample from each `ident.1`
  identity. Will be performed after `downsample.object.last`. If `NULL`,
  all cells are used.

- automatic.downsample:

  Logical. If `TRUE`, automatically downsamples the `seurat_object`
  (respecting the relative proportion of each `ident.1` identity, it is
  therefore similar to `downsample.object.first` and
  `downsample.object.last`) so that the Reference Matrix file written to
  disk would be just under the `max.matrix.size` limit (empirically
  estimated). Performed last. Ignored if `check.size` = `TRUE` or
  `write` = `FALSE`. Please [report an
  issue](https://github.com/Alexis-Varin/RightOmicsTools/issues) if you
  see a significant difference between the file size written to disk and
  `max.matrix.size` (for example, `max.matrix.size` is set to 200 MB but
  the file is 400 MB or 100 MB).

- check.size:

  Logical. If `TRUE`, prints the estimated size of the Reference Matrix
  file that would be written to disk and the number of cells to
  downsample if need be.

- max.matrix.size:

  Numeric. The maximum size of the Reference Matrix file written to disk
  in MB. Will stop the function if the Reference Matrix file written to
  disk is estimated to be over this limit, or if `automatic.downsample`
  = `TRUE`, will downsample the `seurat_object` instead so that the
  Reference Matrix written to disk is under the size limit. Ignored if
  `write` = `FALSE`.

- cell.barcodes:

  Logical. Must be `FALSE` for CIBERSORTx. If `TRUE`, keeps the cell
  barcodes and does not rename with cell identities, if you want to
  extract the expression matrix for projects other than CIBERSORTx.

- output.cell.barcodes:

  Logical. If `TRUE`, outputs the cell barcodes used to build the
  Reference Matrix file.

- file.name:

  Character. The name of the Reference Matrix file written to disk. Must
  not contain any space for CIBERSORTx, the function will automatically
  replace any space with an underscore. Ignored if `check.size` = `TRUE`
  or `write` = `FALSE`.

- file.format:

  Character. The format of the Reference Matrix file written to disk.
  Must be 'txt' or 'tsv' for CIBERSORTx, but you may also specify 'csv'
  for example, if you want to extract the expression matrix for projects
  other than CIBERSORTx. Accepts any format the
  [`fwrite`](https://rdrr.io/pkg/data.table/man/fwrite.html) function
  would accept. Ignored if `check.size` = `TRUE` or `write` = `FALSE`.

- file.sep:

  Character. The separator to use in the Reference Matrix file written
  to disk. Must be tabulation for CIBERSORTx, but you may also specify a
  comma for example, if you want to extract the expression matrix for
  projects other than CIBERSORTx. Accepts any separator the
  [`fwrite`](https://rdrr.io/pkg/data.table/man/fwrite.html) function
  would accept. Ignored if `check.size` = `TRUE` or `write` = `FALSE`.

- write.path:

  Character. The path to write the Reference Matrix into. If `NULL`,
  uses current working directory. Ignored if `check.size` = `TRUE` or
  `write` = `FALSE`.

- write:

  Logical. If `TRUE`, writes to disk the Reference Matrix file. Ignored
  if `check.size` = `TRUE`.

- verbose:

  Logical. If `FALSE`, does not print progress messages and output, but
  warnings and errors will still be printed.

## Value

A `data.table` object containing the normalized counts to CPM from the
`seurat_object` or any other specified `assay` and `layer`, with cell
identities or barcodes as column names and feature names as first
column, or a `list` containing the `data.table` object as well as the
cell barcodes used to build the Reference Matrix file. If `write` =
`TRUE`, the `data.table` object is also written to disk. If `check.size`
= `TRUE`, will instead print the estimated size of the Reference Matrix
file that would be written to disk and the number of cells to downsample
if need be.

## Examples

``` r
# Prepare data
pbmc3k <- Right_Data("pbmc3k")

# Example 1:
refmat <- Reference_Matrix_Builder(pbmc3k,
                                   ident.1 = "seurat_annotations",
                                   ident.2 = "orig.ident",
                                   downsample.object.first = 1500,
                                   subset.ident.1 = c("Naive CD4 T",
                                                      "Memory CD4 T",
                                                      "CD8 T",
                                                      "NK"),
                                   subset.ident.2 = "Donor_2",
                                   subset.ident.2.invert = TRUE,
                                   write = FALSE)
#> Starting... 
#> Downsampling... 
#> Subsetting... 
#> Subsetting... 
#> Normalizing and extracting the expression matrix... 
#> Building the data.table... 
#> Cleaning... 
#> Done. 
refmat[1:5, 1:5]
#>             Gene Memory CD4 T_Donor_3 NK_Donor_3 Memory CD4 T_Donor_1
#>           <char>                <num>      <num>                <num>
#> 1:    AL627309.1                    0          0                    0
#> 2:    AP006222.2                    0          0                    0
#> 3: RP11-206L10.2                    0          0                    0
#> 4: RP11-206L10.9                    0          0                    0
#> 5:     LINC00115                    0          0                    0
#>    CD8 T_Donor_1
#>            <num>
#> 1:             0
#> 2:             0
#> 3:             0
#> 4:             0
#> 5:             0

# Example 2:
Reference_Matrix_Builder(pbmc3k,
                         check.size = TRUE,
                         max.matrix.size = 40)
#> Starting... 
#> Normalizing and extracting the expression matrix... 
#> Current estimated Reference Matrix file size on CIBERSORTx web portal is 109 MB :
#> Matrix of 2700 cells by 13714 features
#> You may want to downsample your Seurat object to 988 cells for a Reference Matrix file under 40 MB
#> [1] 988
```
