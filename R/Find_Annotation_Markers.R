#' @title Get the top markers for each identity in a Seurat object
#'
#' @description This function is a wrapper around \code{\link[Seurat]{FindMarkers}} that allows for parallelization and filtering of mitochondrial, ribosomal and non-coding RNA features in human, as well as filtering of pseudogenes in mouse. It will also directly return the top X markers for each identity.
#'
#' @param seurat_object A \pkg{Seurat} object.
#' @param ident.1 Character. (from \code{\link[Seurat]{FindMarkers}} documentation) Identity class to define markers for; pass an object of class \code{phylo} or 'clustertree' to find markers for a node in a cluster tree; passing 'clustertree' requires \code{\link[Seurat]{BuildClusterTree}} to have been run. Leave \code{NULL} to find markers for all clusters.
#' @param ident.2 Character. (from \code{\link[Seurat]{FindMarkers}} documentation) A second identity class for comparison; if \code{NULL}, use all other cells for comparison; if an object of class \code{phylo} or 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for.
#' @param min.pct Numeric. (from \code{\link[Seurat]{FindMarkers}} documentation) Only test features that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing features that are very infrequently expressed.
#' @param top.markers Numeric. The number of top markers to return. If set to \code{Inf}, all markers will be returned.
#' @param unique.markers Logical. If \code{TRUE}, unique markers will be returned for each identity in order to prevent features repeated multiple times.
#' @param filter.mito Logical. If \code{TRUE}, mitochondrial features will be filtered out.
#' @param filter.ribo Logical. If \code{TRUE}, ribosomal features will be filtered out.
#' @param filter.ncRNA Logical. If \code{TRUE}, non-coding RNA features will be filtered out.
#' @param species Character. The species from which to pull data from to filter out features. If 'human', non-coding RNA features will be filtered out from a dataset named ncRNA_human built from \href{https://www.genenames.org/data/genegroup/#!/group/475}{genenames database}. If 'mouse', only pseudogenes will be filtered out based on a dataset named pseudogenes_mouse and built from \href{https://rna.sysu.edu.cn/dreamBase2/scrna.php?SClade=mammal&SOrganism=mm10&SDataId=0&SProteinID=0}{dreamBase2 database}. These datasets are loaded with \pkg{RightSeuratTools} and may be checked for more information.
#' @param parallelized Logical. If \code{TRUE}, \code{\link[Seurat]{FindMarkers}} will be parallelized using \pkg{BiocParallel}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object to be used for parallelization. If \code{NULL}, the function will set this parameter to \code{\link[BiocParallel]{SerialParam}}, which uses a single worker (core) and is therefore not parallelized, in order to prevent accidental use of large computation resources. Ignored if \code{parallelized} = \code{FALSE}.
#' @param output.df Logical. If \code{TRUE}, a data frame of features names and associated statistics will be returned. If \code{FALSE}, a character vector of features names will be returned.
#' @param output.list Logical. If \code{TRUE}, a list of data frames for each identity with features names and statistics or a list of character vectors containing features names if \code{output.df} = \code{FALSE} will be returned.
#' @param verbose Logical. If \code{FALSE}, does not print progress messages and output, but warnings and errors will still be printed.
#' @param ... Additional arguments to be passed to \code{\link[Seurat]{FindMarkers}}, such as \code{test.use}, or passed to other methods and to specific DE methods.
#'
#' @return A data frame or a list of data frames with features names and associated statistics, or a character vector or a list of character vectors with features names.
#'
#' @examples
#' library(Seurat)
#'
#' # Original object from SeuratObject:
#'
#' pbmc_small
#'
#' # Example of parallelized FindAllMarkers
#' # (provided you set your future::plan() to multisession or multicore):
#'
#' Find_Annotation_Markers(pbmc_small,
#'                        min.pct = 0.01,
#'                        top.markers = Inf,
#'                        unique.markers = FALSE,
#'                        filter.mito = FALSE,
#'                        filter.ribo = FALSE,
#'                        filter.ncRNA = FALSE,
#'                        parallelized = TRUE,
#'                        output.list = FALSE)
#' @import Seurat
#' @import SeuratObject
#' @import BiocParallel
#' @export

Find_Annotation_Markers = function(seurat_object,
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
                                   output.df = FALSE,
                                   output.list = FALSE,
                                   verbose = TRUE,
                                   ...) {

  if (species != "human" & species != "mouse") {
    stop("Currently the function only accepts either 'human' or 'mouse' as species")
  }

  if (species == "human") {
    to.remove = RightSeuratTools::ncRNA_human
  }
  if (species == "mouse") {
    to.remove = RightSeuratTools::pseudogenes_mouse
  }

  if (isTRUE(parallelized)) {
    if (is.null(BPPARAM)) {
      warning("No BPPARAM parameter provided, using BiocParallel::SerialParam(), which is not parallelized")
      BPPARAM = SerialParam()
    }
    if (is.null(ident.1) & is.null(ident.2)) {
      all.markers2 = list()
      idents = levels(Idents(seurat_object))
      all.markers2 = bplapply(idents, function(x) {
        if (isTRUE(verbose)) {
          cat("Finding markers for cluster ",x," against all other clusters","\n",sep="")
        }
        FindMarkers(object = seurat_object, ident.1 = x, min.pct = min.pct, ...)
      }, BPPARAM = BPPARAM)
    }

    if (is.null(ident.1) & !is.null(ident.2)) {
      all.markers2 = list()
      idents = levels(Idents(seurat_object))
      all.markers2 = bplapply(idents, function(x) {
        if (isTRUE(verbose)) {
          cat("Finding markers for cluster ",x," against cluster ",ident.2,"\n",sep="")
        }
        FindMarkers(object = seurat_object, ident.1 = x, ident.2 = ident.2, min.pct = min.pct, ...)
      }, BPPARAM = BPPARAM)
    }

    if (is.null(ident.2) & !is.null(ident.1)) {
      all.markers2 = list()
      idents = levels(Idents(seurat_object))
      all.markers2 = bplapply(idents, function(x) {
        if (isTRUE(verbose)) {
          cat("Finding markers for cluster ",ident.1," against cluster ",x,"\n",sep="")
        }
        FindMarkers(object = seurat_object, ident.1 = ident.1, ident.2 = x, min.pct = min.pct, ...)
      }, BPPARAM = BPPARAM)
    }

    if (!is.null(ident.2) & !is.null(ident.1)) {
      all.markers2 = list()
      idents = levels(Idents(seurat_object))
      if (isTRUE(verbose)) {
        cat("Finding markers between cluster ",ident.1," and cluster ",ident.2,"\n",sep="")
      }
      all.markers2 = FindMarkers(object = seurat_object, ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, ...)
    }
  }

  else {
    if (is.null(ident.1) & is.null(ident.2)) {
      all.markers2 = list()
      idents = levels(Idents(seurat_object))
      all.markers2 = lapply(idents, function(x) {
        if (isTRUE(verbose)) {
          cat("Finding markers for cluster ",x," against all other clusters","\n",sep="")
        }
        FindMarkers(object = seurat_object, ident.1 = x, min.pct = min.pct, ...)
      })
    }

    if (is.null(ident.1) & !is.null(ident.2)) {
      all.markers2 = list()
      idents = levels(Idents(seurat_object))
      all.markers2 = lapply(idents, function(x) {
        if (isTRUE(verbose)) {
          cat("Finding markers for cluster ",x," against cluster ",ident.2,"\n",sep="")
        }
        FindMarkers(object = seurat_object, ident.1 = x, ident.2 = ident.2, min.pct = min.pct, ...)
      })
    }

    if (is.null(ident.2) & !is.null(ident.1)) {
      all.markers2 = list()
      idents = levels(Idents(seurat_object))
      all.markers2 = lapply(idents, function(x) {
        if (isTRUE(verbose)) {
          cat("Finding markers for cluster ",ident.1," against cluster ",x,"\n",sep="")
        }
        FindMarkers(object = seurat_object, ident.1 = ident.1, ident.2 = x, min.pct = min.pct, ...)
      })
    }

    if (!is.null(ident.2) & !is.null(ident.1)) {
      all.markers2 = list()
      idents = levels(Idents(seurat_object))
      if (isTRUE(verbose)) {
        cat("Finding markers between cluster ",ident.1," and cluster ",ident.2,"\n",sep="")
      }
      all.markers2 = FindMarkers(object = seurat_object, ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, ...)
    }
  }

  all.markers = data.frame()
  top.markers.df = data.frame()
  for (i in 1:length(all.markers2)) {
    all.markers2[[i]]$cluster = levels(Idents(seurat_object))[i]
    all.markers2[[i]]$feature = gsub("\\..[0-9]*","",rownames(all.markers2[[i]]))
    all.markers2[[i]] = all.markers2[[i]][order(all.markers2[[i]]$avg_log2FC, decreasing = T),]

    if (isTRUE(filter.ncRNA)) {
      all.markers2[[i]] = all.markers2[[i]][!all.markers2[[i]]$feature %in% to.remove,]
      if (species == "human") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^A[C,L,P][0-9]|^LINC[0-9]|^LNC", x = all.markers2[[i]]$feature),]
      }
    }

    if (isTRUE(filter.mito)) {
      if (species == "human") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^MT-", x = all.markers2[[i]]$feature),]
      }
      if (species == "mouse") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^mt-", x = all.markers2[[i]]$feature),]
      }
    }

    if (isTRUE(filter.ribo)) {
      if (species == "human") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^RP[SL]", x = all.markers2[[i]]$feature),]
      }
      if (species == "mouse") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^Rp[sl]", x = all.markers2[[i]]$feature),]
      }
    }

    all.markers = rbind(all.markers, all.markers2[[i]])

    if (!is.infinite(top.markers)) {
      if (isFALSE(unique.markers)) {
        select.top.markers = all.markers2[[i]]
      }
      else {
        select.top.markers = all.markers2[[i]][setdiff(all.markers2[[i]]$feature,top.markers.df$feature),]
      }
      top.markers.df = rbind(top.markers.df,select.top.markers[1:top.markers,])
      all.markers2[[i]] = select.top.markers[1:top.markers,]
    }
  }

  if (isTRUE(output.df) & isTRUE(output.list)) {
    if (!is.infinite(top.markers)) {
      final.list = list("data frame" = top.markers.df, "list" = all.markers2)
    }
    else {
      final.list = list("data frame" = all.markers, "list" = all.markers2)
    }
    return(final.list)
  }
  if (isTRUE(output.df) & isFALSE(output.list)) {
    if (!is.infinite(top.markers)) {
      return(top.markers.df)
    }
    else {
      return(all.markers)
    }
  }
  if (isFALSE(output.df) & isTRUE(output.list)) {
    if (!is.infinite(top.markers)) {
      final.list = list("features" = top.markers.df$feature, "list" = all.markers2)
    }
    else {
      final.list = list("features" = all.markers$feature, "list" = all.markers2)
    }
    return(final.list)
  }
  else {
    if (!is.infinite(top.markers)) {
      return(top.markers.df$feature)
    }
    else {
      return(all.markers$feature)
    }
  }
}
