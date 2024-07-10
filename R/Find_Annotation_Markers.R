#' @title Find_Annotation_Markers
#'
#' @description This function is a wrapper around Seurat's FindMarkers function that allows for parallelization and filtering of mitochondrial, ribosomal and non-coding RNA genes in human, as well as filtering of pseudogenes in mouse. It will also directly give the top X markers for each identity to use for plotting with DotPlot_Heatmap() for example.
#'
#' @param seurat_object A Seurat object.
#' @param ident.1 Character. Identity class to define markers for; pass an object of class phylo or 'clustertree' to find markers for a node in a cluster tree; passing 'clustertree' requires BuildClusterTree to have been run. Leave NULL to find markers for all clusters.
#' @param ident.2 Character. A second identity class for comparison; if NULL, use all other cells for comparison; if an object of class phylo or 'clustertree' is passed to ident.1, must pass a node to find markers for.
#' @param min.pct Numeric. Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
#' @param top.markers Numeric. The number of top markers to return. If set to Inf, all markers will be returned.
#' @param unique.markers Logical. If TRUE, unique markers will be returned for each cluster in order to prevent markers repeated in multiple clusters.
#' @param filter.mito Logical. If TRUE, mitochondrial genes will be filtered out.
#' @param filter.ribo Logical. If TRUE, ribosomal genes will be filtered out.
#' @param filter.ncRNA Logical. If TRUE, non-coding RNA genes will be filtered out.
#' @param species Character. The species from which to pull data to filter out genes. If "human", non-coding RNA genes will be filtered out from a dataset named ncRNA_human built from genenames database. If "mouse", only pseudogenes will be filtered out based on a dataset named pseudogenes_mouse and built from dreamBase2 database. These datasets are loaded with RightSeuratTools package and may be checked for more information.
#' @param parallelized Logical. If TRUE, FindMarkers will be parallelized using BiocParallel.
#' @param BPPARAM A BiocParallelParam object to be used for parallelization. If NULL, will use SerialParam() which is not parallelized. Ignored if parallelized = FALSE.
#' @param output.df Logical. If TRUE, a data frame of gene names and associated statistics will be returned. If FALSE, a character vector of gene names will be returned.
#' @param output.list Logical. If TRUE, a list of gene names with or without statistics for each cluster will be returned.
#' @param verbose Logical. If FALSE, does not print progress messages and output, but warnings and errors will still be printed.
#' @param ... Additional arguments to be passed to FindMarkers, such as test.use, or other methods and to specific DE methods downstream of FindMarkers.
#'
#' @return A data frame or a list of data frames with gene names and associated statistics, or a character vector or a list of character vectors with gene names.
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
                                   output.df = TRUE,
                                   output.list = FALSE,
                                   verbose = TRUE,
                                   ...) {

  if (species != "human" & species != "mouse") {
    stop("Currently the function only accepts either 'human' or 'mouse' as species")
  }

  if (species == "human") {
    to.remove = ncRNA_human
  }
  if (species == "mouse") {
    to.remove = pseudogenes_mouse
  }

  if (isTRUE(parallelized)) {
    if (is.null(BPPARAM)) {
      warning("No BPPARAM parameter provided, using SerialParam(), which is not parallelized")
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
    all.markers2[[i]]$gene = gsub("\\..[0-9]*","",rownames(all.markers2[[i]]))
    all.markers2[[i]] = all.markers2[[i]][order(all.markers2[[i]]$avg_log2FC, decreasing = T),]

    if (isTRUE(filter.ncRNA)) {
      all.markers2[[i]] = setdiff(all.markers2[[i]], all.markers2[[i]][which(all.markers2[[i]] %in% to.remove)])
      if (species == "human") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^A[C,L,P][0-9]|^LINC[0-9]|^LNC", x = all.markers2[[i]]$gene),]
      }
    }

    if (isTRUE(filter.mito)) {
      if (species == "human") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^MT-", x = all.markers2[[i]]$gene),]
      }
      if (species == "mouse") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^mt-", x = all.markers2[[i]]$gene),]
      }
    }

    if (isTRUE(filter.ribo)) {
      if (species == "human") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^RP[SL]", x = all.markers2[[i]]$gene),]
      }
      if (species == "mouse") {
        all.markers2[[i]] = all.markers2[[i]][!grepl(pattern = "^Rp[sl]", x = all.markers2[[i]]$gene),]
      }
    }

    all.markers = rbind(all.markers, all.markers2[[i]])

    if (!is.infinite(top.markers)) {
      if (isFALSE(unique.markers)) {
        top.markers.df = rbind(top.markers.df,all.markers2[[i]][1:length(top.markers[,1]),])
      }
      if (isTRUE(unique.markers)) {
        j = 1
        k = 0
        while (j < top.markers+1) {
          k = k + 1
          if (isFALSE(any(grepl(all.markers2[[i]]$gene[k],top.markers.df)))) {
            top.markers.df = rbind(top.markers.df,all.markers2[[i]][k,])
            j = j + 1
          }
        }
      }
      all.markers2[[i]] = top.markers.df[length(top.markers.df[,1])-(top.markers-1):length(top.markers.df[,1]),]
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
      final.list = list("features" = top.markers.df$gene, "list" = all.markers2)
    }
    else {
      final.list = list("features" = all.markers$gene, "list" = all.markers2)
    }
    return(final.list)
  }

  if (isFALSE(output.df) & isFALSE(output.list)) {
    if (!is.infinite(top.markers)) {
      return(top.markers.df$gene)
    }
    else {
      return(all.markers$gene)
    }
  }
}
