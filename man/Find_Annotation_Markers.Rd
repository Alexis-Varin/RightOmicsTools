% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Find_Annotation_Markers.R
\name{Find_Annotation_Markers}
\alias{Find_Annotation_Markers}
\title{Get the top markers for fast annotation}
\usage{
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
}
\arguments{
\item{seurat_object}{A \pkg{Seurat} object.}

\item{ident.1}{Character. (from \code{\link[Seurat]{FindMarkers}} documentation) Identity class to define markers for; pass an object of class \code{phylo} or 'clustertree' to find markers for a node in a cluster tree; passing 'clustertree' requires \code{\link[Seurat]{BuildClusterTree}} to have been run. Leave \code{NULL} to find markers for all clusters.}

\item{ident.2}{Character. (from \code{\link[Seurat]{FindMarkers}} documentation) A second identity class for comparison; if \code{NULL}, use all other cells for comparison; if an object of class \code{phylo} or 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for.}

\item{min.pct}{Numeric. (from \code{\link[Seurat]{FindMarkers}} documentation) Only test features that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing features that are very infrequently expressed.}

\item{top.markers}{Numeric. The number of top markers to return. If \code{Inf}, all markers will be returned.}

\item{unique.markers}{Logical. If \code{TRUE}, unique markers will be returned for each identity in order to prevent features repeated multiple times.}

\item{filter.mito}{Logical. If \code{TRUE}, mitochondrial features will be filtered out.}

\item{filter.ribo}{Logical. If \code{TRUE}, ribosomal features will be filtered out.}

\item{filter.ncRNA}{Logical. If \code{TRUE}, non-coding RNA features will be filtered out.}

\item{species}{Character. The species name to filter out non-coding RNA features. If 'human', a dataset named \href{https://alexis-varin.github.io/RightOmicsTools/reference/ncRNA_human.html}{ncRNA_human} built from \href{https://www.genenames.org/data/genegroup/#!/group/475}{genenames database} will be used as reference. If 'mouse', only pseudogenes will be filtered out based on a dataset named \href{https://alexis-varin.github.io/RightOmicsTools/reference/pseudogenes_mouse.html}{pseudogenes_mouse} and built from \href{https://rna.sysu.edu.cn/dreamBase2/scrna.php?SClade=mammal&SOrganism=mm10&SDataId=0&SProteinID=0}{dreamBase2 database}. These datasets are loaded with \pkg{RightOmicsTools} and may be checked for more information.}

\item{parallelized}{Logical. If \code{TRUE}, \code{\link[Seurat]{FindMarkers}} will be parallelized using \pkg{BiocParallel}. Please note that parallelization is complex and depends on your operating system (Windows users might not see a gain or might even experience a slowdown).}

\item{BPPARAM}{A \code{\link[BiocParallel]{BiocParallelParam}} object to be used for parallelization. If \code{NULL} and \code{parallelized} = \code{TRUE}, the function will use a \code{\link[BiocParallel]{SerialParam}} object configured to use a single worker (core) and is therefore not parallelized, in order to prevent accidental use of large computation resources. Ignored if \code{parallelized} = \code{FALSE}.}

\item{name.markers}{Logical. If \code{TRUE}, each marker will be named with its associated identity. Ignored if \code{output.df} = \code{TRUE}.}

\item{output.df}{Logical. If \code{TRUE}, a \code{data.frame} object containing the \code{\link[Seurat]{FindMarkers}} results, ordered by decreasing log fold change will be returned. If \code{FALSE}, feature names will be returned.}

\item{output.list}{Logical. If \code{TRUE}, and if \code{output.df} = \code{TRUE}, returns a \code{list} of \code{data.frame} objects containing \code{\link[Seurat]{FindMarkers}} results for each identity, or if \code{output.df} = \code{FALSE}, returns a \code{list} of feature names for each identity.}

\item{verbose}{Logical. If \code{FALSE}, does not print progress messages and output, but warnings and errors will still be printed.}

\item{...}{Additional arguments to be passed to \code{\link[Seurat]{FindMarkers}}, such as \code{test.use}, or passed to other methods and to specific DE methods.}
}
\value{
A \code{data.frame} object or a \code{list} of \code{data.frame} objects containing the \code{\link[Seurat]{FindMarkers}} results, or a \code{list} of features names, or feature names.
}
\description{
This function is a wrapper around \code{\link[Seurat]{FindMarkers}} that allows for parallelization and filtering of mitochondrial, ribosomal and non-coding RNA features in human, as well as filtering of pseudogenes in mouse. It may also directly return the top X markers for each identity.
}
\examples{
\dontshow{
suppressWarnings(suppressPackageStartupMessages(library(Seurat)))
options(Seurat.presto.wilcox.msg = FALSE)
}
# Prepare data
pbmc3k <- Right_Data("pbmc3k")

# Example 1: default parameters and origin of markers
pbmc3k.markers <- Find_Annotation_Markers(pbmc3k,
                                          name.markers = TRUE)
pbmc3k.markers

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
head(pbmc3k.markers)
}
