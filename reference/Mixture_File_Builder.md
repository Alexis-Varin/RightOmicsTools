# Builder of CIBERSORTx Mixture File

This function builds a Mixture File for
[CIBERSORTx](https://cibersortx.stanford.edu/index.php) from bulk
RNA-seq and/or microarray files and/or R objects.

## Usage

``` r
Mixture_File_Builder(
  objects = NULL,
  files.path = NULL,
  file.name = "Mixture_File",
  file.format = "txt",
  file.sep = "\t",
  write.path = NULL,
  write = TRUE,
  verbose = TRUE
)
```

## Arguments

- objects:

  A `data.frame`, `data.table` or `matrix` object or a mixed `list` of
  `data.frame`, `data.table` and/or `matrix` objects corresponding to
  bulk RNA-seq and/or microarray experiments, with genes as rows and
  gene names as row names or first column. A single object may contain
  multiple samples and/or experiments as columns as long as it contains
  a single gene names column as row names or first column (i.e. \| gene
  \| sample1 \| sample2 \| sample3 \| etc). You may mix objects with
  different gene names order and number (for example, an object
  containing 13567 genes and another with 18458 genes).

- files.path:

  Character. The path to read txt, csv and/or tsv files corresponding to
  bulk RNA-seq or microarray experiments, with genes as rows and gene
  names as first column. A single file may contain multiple samples
  and/or experiments as columns as long as it contains a single gene
  names column as first column. You may mix files with different genes
  names order and number (for example, a file containing 21848 genes and
  another with 19457 genes).

- file.name:

  Character. The name of the Mixture File written to disk. Must not
  contain any space for CIBERSORTx, the function will automatically
  replace any space with an underscore. Ignored if `write` = `FALSE`.

- file.format:

  Character. The format of the Mixture File written to disk. Must be
  'txt' or 'tsv' for CIBERSORTx, but you may also specify 'csv' for
  example, for projects other than CIBERSORTx. Accepts any format the
  [`fwrite`](https://rdatatable.gitlab.io/data.table/reference/fwrite.html)
  function would accept. Ignored if `write` = `FALSE`.

- file.sep:

  Character. The separator to use in the Mixture File written to disk.
  Must be tabulation for CIBERSORTx, but you may also specify a comma
  for example, for projects other than CIBERSORTx. Accepts any separator
  the
  [`fwrite`](https://rdatatable.gitlab.io/data.table/reference/fwrite.html)
  function would accept. Ignored if `write` = `FALSE`.

- write.path:

  Character. The path to write the Mixture File into. If `NULL`, uses
  current working directory. Ignored if `write` = `FALSE`.

- write:

  Logical. If `TRUE`, writes to disk the Mixture File.

- verbose:

  Logical. If `FALSE`, does not print progress messages and output, but
  warnings and errors will still be printed.

## Value

A `data.table` object containing the counts of the objects and/or files
provided, with each column being a sample. If `write` = `TRUE`, the
`data.table` object is also written to disk.
