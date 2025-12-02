# Remove genes with very low expression for tradeSeq analysis

The functions [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html)
and [`evaluateK`](https://rdrr.io/pkg/tradeSeq/man/evaluateK.html), used
to analyze gene expression along trajectories, usually require both
considerable calculation power (access to computer clusters) and time
(up to days) to complete. While one might assume that computing time is
uniform across genes, it appears however that model fitting takes a few
seconds or less for the vast majority of genes, and that only a small
fraction disproportionately affects calculation time, with model fitting
taking up to hours or even days for a single gene. Further examination
shows that all these genes are among the lowest expressed, and that
there is an inverse correlation between the number of cells with a
non-zero count for a particular gene and model fitting time, which
becomes exponentially longer as the number of expressing cells tends
towards zero. Based on these observations, this function was built to
considerably speed up computing time / lower computing power
requirements while retaining as much information as possible. It filters
out genes below a certain expression threshold (10 expressing cells by
default), effectively removing genes whose expression is likely just
noise (and thus the cause of lengthy and difficult model fitting),
rather than a signal of biological significance. This is especially
common when cells are subset from a larger object for trajectory
inference, and the expression matrix contains many zero or near-zero
expression genes, as these would be primarily expressed in excluded
cells. Additionally, it offers a second method of filtering, if
computing resources are still a limiting factor and further reduction in
the number of genes is needed, by retaining only the top genes according
to deviance, which provides a better information estimate than selecting
variable genes, as the latter tends to still include a few genes with
near-zero expressing cells. Finally, if
[`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) and
[`evaluateK`](https://rdrr.io/pkg/tradeSeq/man/evaluateK.html) will be
run with a condition, it also allows filtering based on the specific
metadata that will be used.

## Usage

``` r
tradeSeqPreprocess(
  sds,
  min.cells = 10,
  nb.genes = NULL,
  condition = NULL,
  filter.mito = FALSE,
  filter.ribo = FALSE,
  filter.ncRNA = FALSE,
  species = "human",
  plot.genes = TRUE,
  output.data = FALSE
)
```

## Arguments

- sds:

  A SingleCellExperiment object containing the pseudotime values of one
  or several lineages, computed using
  [`slingshot`](https://rdrr.io/pkg/slingshot/man/slingshot.html)
  (usually, the input object to
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html)).

- min.cells:

  Numeric. The minimum number of cells with 1 count or more to keep a
  gene.

- nb.genes:

  Numeric. The maximum number of genes to keep according to deviance, as
  computed by
  [`devianceFeatureSelection`](https://rdrr.io/pkg/scry/man/devianceFeatureSelection.html)
  internally. If `NULL`, all genes are kept.

- condition:

  Character or Factor. Either the name of a metadata present in `sds`
  (for example, 'treatment', 'disease', etc) to subset the object and
  separately remove lowly expressed genes, or the identities, as
  character or factor, of length equal to the number of cells in `sds`.

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

- plot.genes:

  Logical. If `TRUE`, a bar plot showing the number of genes removed
  according to the number of cells with at least 1 count or more is
  displayed.

- output.data:

  Logical. If `TRUE`, the function will return a `list` containing the
  filtered SingleCellExperiment object, a `data.frame` object with the
  number of cells with at least 1 count or more for each removed gene,
  and a `data.frame` object with the data used to build the bar plot.

## Value

A filtered SingleCellExperiment object, or a `list` containing the
filtered SingleCellExperiment object, a `data.frame` object with the
number of cells with at least 1 count or more for each removed gene, and
a `data.frame` object with the data used to build the bar plot.
