#' @title Right_DietSeurat
#'
#' @description As of \pkg{Seurat} v5 release, \code{\link[Seurat]{DietSeurat}} does not remove data and scale.data layers, resulting in objects with little to no reduction in size. This function was created with the purpose to restore \code{\link[Seurat]{DietSeurat}}'s proper behavior until the function is fixed by \pkg{Seurat}'s dev team, as well as to offer a few new options such as the ability to subset cells or keep variable features.
#'
#' @param seurat_object A \pkg{Seurat} object.
#' @param idents Character. A vector with one or several identity names in the meta.data slot to keep in the diet object. If \code{NULL}, all identities will be kept.
#' @param cells Character. A vector with one or several cell barcodes to keep in the diet object. If \code{NULL}, all cells will be kept.
#' @param features Character. A vector with one or several feature names to keep in the diet object. If \code{NULL}, all features will be kept.
#' @param dimreducs Character. A vector with one or several reduction names to keep in the diet object. If \code{NULL}, all DimReducs will be removed. Note that if you subset features and remove all from which PCs were calculated, the PCA will not be kept, even if it is set in this parameter, as the feature.loadings slot will be empty, this is intended behavior.
#' @param graphs Character. A vector with one or several graph names to keep in the diet object. If \code{NULL}, all Graphs will be removed.
#' @param variable.features Logical. If \code{TRUE}, the variable features slot will be kept in the diet object. Note that if you keep only a subset of features, the variable features will also be subset, and if you remove all features from which variable features were calculated, the variable features will not be kept even if this parameter is set to \code{TRUE}, this is intended behavior.
#' @param misc Logical. If \code{TRUE}, the misc slot will be kept in the diet object.
#' @param split.counts.layer Character. The name of an identity in the meta.data slot to split the counts layer by if you need to in your downstream analysis. If \code{NULL}, the diet object will have a single counts layer, even if the original object had split counts layers.
#' @param data.layer Logical. If \code{TRUE}, the data layer of the RNA assay will be kept in the diet object if it is present. As with the counts layer, if there are split data layers, they will be joined into a single data layer unless the \code{split.counts.layer} parameter is set.
#' @param scale.layer Logical. If \code{TRUE}, the scale.data layer of the RNA assay will be kept in the diet object if it is present.
#' @param SCTAssay Logical. If \code{TRUE}, the SCT assay will be kept in the diet object if it is present.
#'
#' @return A \pkg{Seurat} object, hopefully smaller, with class Assay5 RNA assay and specified layers and slots.
#'
#' @examples
#' library(Seurat)
#'
#' # Original object from SeuratObject:
#'
#' pbmc_small
#'
#' # Example 1:
#'
#' Right_DietSeurat(pbmc_small,
#'                 idents = c("letter.idents","groups"),
#'                 features = c("CD3D", "NOSIP", "SAFB2",
#'                  "CD2", "IL7R", "PIK3IP1", "MAL"),
#'                 graphs = "RNA_snn",
#'                 split.counts.layer = "letter.idents")
#'
#' # Example 2:
#'
#' pbmc_cells = WhichCells(pbmc_small, idents = c("0","1"))
#'
#' Right_DietSeurat(pbmc_small,
#'                 cells = pbmc_cells,
#'                 variable.features = TRUE,
#'                 dimreducs = "pca",
#'                 data.layer = TRUE)
#'
#' @import Seurat
#' @import SeuratObject
#' @export

Right_DietSeurat = function(seurat_object,
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
                            SCTAssay = FALSE) {

  for (i in 1:length(seurat_object@meta.data)) {
    if (i == 1) {
      idents.df = data.frame("1" = seurat_object@meta.data[i])
    }
    else {
      idents.df = cbind(idents.df, seurat_object@meta.data[i])
    }
  }
  colnames(idents.df) = names(seurat_object@meta.data)

  if (is.character(idents)) {
    idents.df = idents.df[ ,colnames(idents.df) %in% idents]
  }

  if (length(Layers(seurat_object[["RNA"]], search = "counts")) > 1) {
    seurat_object[["RNA"]] = JoinLayers(seurat_object[["RNA"]])
  }

  new_object = CreateSeuratObject(counts = LayerData(seurat_object, assay = "RNA", layer = "counts"),
                                  project = paste0(seurat_object@project.name,"_Diet"),
                                  meta.data = idents.df, min.cells = 0, min.features = 0)

  Idents(new_object) = Idents(seurat_object)

  if (isTRUE(misc)) {
    new_object@misc = seurat_object@misc
  }

  if (isTRUE(data.layer)) {
    new_object[["RNA"]]$data = seurat_object[["RNA"]]$data
  }

  if (isTRUE(scale.layer)) {
    new_object[["RNA"]]$scale.data = seurat_object[["RNA"]]$scale.data
  }

  if (isTRUE(SCTAssay)) {
    new_object[["SCT"]] = seurat_object[["SCT"]]
  }

  if (is.character(c(dimreducs, graphs))) {
    dimreducs = c(dimreducs, graphs)
    for (i in 1:length(dimreducs)) {
      new_object[[dimreducs[i]]] = seurat_object[[dimreducs[i]]]
    }
  }

  if (is.character(features)) {
    new_object = subset(new_object, features = features)
  }

  if (isTRUE(variable.features)) {
    VariableFeatures(new_object) = intersect(VariableFeatures(seurat_object),rownames(new_object[["RNA"]]$counts))
  }

  if (is.character(cells)) {
    new_object = subset(new_object, cells = cells)
  }

  if (is.character(split.counts.layer)) {
    new_object[["RNA"]] = split(new_object[["RNA"]], f = new_object@meta.data[,split.counts.layer])
  }

  return(new_object)

}
