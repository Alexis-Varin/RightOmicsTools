# Perform tradeSeq tests in parallel

This function is a wrapper around tradeSeq test functions, such as
[`associationTest`](https://rdrr.io/pkg/tradeSeq/man/associationTest.html)
or
[`conditionTest`](https://rdrr.io/pkg/tradeSeq/man/conditionTest.html),
and allows multiple tests to be performed simultaneously, with all
necessary parameters specified internally. It also organizes the output
of each test into separate `list` elements, making it easy to access
results. Finally, it calculates an adjusted p-value and orders based on
decreasing Waldstat.

## Usage

``` r
tradeSeqTests(
  models,
  tests = c("assoc", "pattern", "diffend", "early", "startvsend", "condition"),
  global = TRUE,
  pairwise = TRUE,
  lineages = TRUE,
  l2fc = 0,
  eigen.thresh = 0.01,
  n.points = 2 * nknots(models),
  knots = NULL,
  pseudotime.values = NULL,
  parallelized = FALSE,
  BPPARAM = NULL,
  tidy = TRUE,
  raw = FALSE,
  verbose = TRUE
)
```

## Arguments

- models:

  A SingleCellExperiment object containing the fitted GAM smoothers,
  computed using
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html), with or
  without conditions provided.

- tests:

  Character. The names of one or several tests to perform. You may
  provide partial names, such as 'assoc' for
  [`associationTest`](https://rdrr.io/pkg/tradeSeq/man/associationTest.html),
  as long as a single match is found. The available tests are
  [`associationTest`](https://rdrr.io/pkg/tradeSeq/man/associationTest.html),
  [`patternTest`](https://rdrr.io/pkg/tradeSeq/man/patternTest.html),
  [`diffEndTest`](https://rdrr.io/pkg/tradeSeq/man/diffEndTest.html),
  [`earlyDETest`](https://rdrr.io/pkg/tradeSeq/man/earlyDETest.html),
  [`startVsEndTest`](https://rdrr.io/pkg/tradeSeq/man/startVsEndTest.html)
  and
  [`conditionTest`](https://rdrr.io/pkg/tradeSeq/man/conditionTest.html).
  [`conditionTest`](https://rdrr.io/pkg/tradeSeq/man/conditionTest.html)
  is ignored if `models` was computed using
  [`fitGAM`](https://rdrr.io/pkg/tradeSeq/man/fitGAM.html) without
  conditions.

- global:

  Logical. (from
  [`earlyDETest`](https://rdrr.io/pkg/tradeSeq/man/earlyDETest.html)/[`startVsEndTest`](https://rdrr.io/pkg/tradeSeq/man/startVsEndTest.html)
  documentation) If `TRUE`, test for all pairwise comparisons/lineages
  simultaneously.

- pairwise:

  Logical. (from
  [`earlyDETest`](https://rdrr.io/pkg/tradeSeq/man/earlyDETest.html)
  documentation) If `TRUE`, test for all pairwise comparisons
  independently. Ignored in tests that do not have `pairwise`.

- lineages:

  Logical. (from
  [`startVsEndTest`](https://rdrr.io/pkg/tradeSeq/man/startVsEndTest.html)
  documentation) If `TRUE`, test for all lineages independently. Ignored
  in tests that do not have `lineages`.

- l2fc:

  Numeric. (from
  [`earlyDETest`](https://rdrr.io/pkg/tradeSeq/man/earlyDETest.html)
  documentation) The log2 fold change threshold to test against. Note,
  that this will affect both the global test and the pairwise
  comparisons.

- eigen.thresh:

  Numeric. (from
  [`earlyDETest`](https://rdrr.io/pkg/tradeSeq/man/earlyDETest.html)
  documentation) Eigenvalue threshold for inverting the
  variance-covariance matrix of the coefficients to use for calculating
  the Wald test statistics. Lower values are more lenient to adding more
  information but also decrease computational stability. This argument
  should in general not be changed by the user but is provided for
  back-compatability. Set to 1e-8 to reproduce results of older version
  of 'tradeSeq'. Ignored in tests that do not have `eigen.thresh`.

- n.points:

  Numeric. (from
  [`earlyDETest`](https://rdrr.io/pkg/tradeSeq/man/earlyDETest.html)
  documentation) The number of points to be compared between lineages.
  Defaults to twice the number of knots. Ignored in tests that do not
  have `n.points`.

- knots:

  Numeric or List. (from
  [`earlyDETest`](https://rdrr.io/pkg/tradeSeq/man/earlyDETest.html)
  documentation) A vector of length 2 specifying the knots at the start
  and end of the region of interest. You may also provide a `list` of
  multiple elements (for example, list(c(2,4), c(3,4))) to repeat the
  test for multiple sets of knots. Ignored in tests that do not have
  `knots`.

- pseudotime.values:

  Numeric or List. (from
  [`startVsEndTest`](https://rdrr.io/pkg/tradeSeq/man/startVsEndTest.html)
  documentation) A vector of length 2, specifying two pseudotime values
  to be compared against each other, for every lineage of the
  trajectory. @details Note that this test assumes that all lineages
  start at a pseudotime value of zero, which is the starting point
  against which the end point is compared. You may also provide a `list`
  of multiple elements (for example, list(c(8,12), c(6,14))) to repeat
  the test for multiple sets of pseudotime values. Ignored in tests that
  do not have `pseudotime.values`.

- parallelized:

  Logical. If `TRUE`, the `tests` will be parallelized using
  BiocParallel. Please note that parallelization is complex and depends
  on your operating system (Windows users might not see a gain or might
  even experience a slowdown).

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object to be used for parallelization. If `NULL` and `parallelized` =
  `TRUE`, the function will use a
  [`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  object configured to use a single worker (core) and is therefore not
  parallelized, in order to prevent accidental use of large computation
  resources. Ignored if `parallelized` = `FALSE`.

- tidy:

  Logical. If `TRUE`, a `list` is returned, with `data.frame` objects
  corresponding to each test divided into each `global`, `pairwise`
  and/or `lineages` comparison results. An adjusted p-value (False
  Discovery Rate) is also calculated and each `data.frame` object is
  ordered based on decreasing Waldstat.

- raw:

  Logical. If `TRUE`, a `list` is returned, with `data.frame` objects
  corresponding to each test results.

- verbose:

  Logical. If `FALSE`, does not print progress messages and output, but
  warnings and errors will still be printed.

## Value

A `list`, with `data.frame` objects corresponding to each test results
and/or `data.frame` objects corresponding to each test and each
`global`, `pairwise` and/or `lineages` comparison results.
