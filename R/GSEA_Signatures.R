#' @title GSEA_Signatures
#'
#' @description This function creates signatures (module scores calculated from \pkg{UCell} or \pkg{Seurat}'s respective functions) from features found in the \pkg{Seurat} object and extracted from supplied MSigDB's pathways.
#'
#' @param seurat_object A \pkg{Seurat} object.
#' @param assay Character. If the \pkg{Seurat} object contains multiple RNA assays, you may specify which one to use (for example 'RNA2' if you have created a second RNA assay you named 'RNA2'. See \href{https://satijalab.org/seurat/articles/seurat5_essential_commands.html}{Seurat v5 vignettes} for more information). You may also use another assay such as 'SCT' to pull features expression from.
#' @param layer Character. Formerly known as slot. If you have split layers the function will always join them before searching features and adding the signatures.
#' @param species Character. The species name to be passed to \code{\link[msigdbr]{msigdbr}} to build the pathways database. Use \code{\link[msigdbr]{msigdbr_species}} for the list of available species.
#' @param category Character. The category or categories to be passed to \code{\link[msigdbr]{msigdbr}} to build the pathways database. Use \code{\link[msigdbr]{msigdbr_collections}} for the list of available categories (gs_cat column). If \code{NULL}, all categories will be used.
#' @param subcategory Character. The subcategory or subcategories to be passed to \code{\link[msigdbr]{msigdbr}} to build the pathways database. Use \code{\link[msigdbr]{msigdbr_collections}} for the list of available subcategories (gs_subcat column). If \code{NULL}, all subcategories will be used.
#' @param pathways Character. The names of the pathways to be used to create the signatures. You may provide either a pathway id (for example, 'GO:0006574') or a name matching the pattern found in gs_name column (all caps and underscores between words). Please note that you may also provide a partial match (for example, 'TYPE_I_INTERFERON') and the function will find all pathways containing this partial string. Beware that this may result in a large number of pathways to be added as signatures (using \code{only.features} = \code{TRUE} is highly recommended) but is very handy to explore all pathways of interest in a particular biological process.
#' @param min.features Numeric. The minimum number of features present in the \pkg{Seurat} object for a pathway to be considered.
#' @param signatures.names Character. Either 'id' which will add pathways ids as signatures (e.g. GO:0004657, hsa05200 etc), or 'name' which will add pathways names as signatures, which might be very long. You may also provide a vector of names to be used as signatures, whose length must match the number of pathways found and kept. It is recommended to use \code{only.features} = \code{TRUE} and set \code{signatures.names} = 'name' or 'id' to get the number of signatures you need to provide names for.
#' @param method Character. The method you want to use to calculate the module scores, either 'UCell' or 'Seurat'.
#' @param only.features Logical. If \code{TRUE}, the function will not add any signature to the \pkg{Seurat} object and will only return the \pkg{Seurat} object as well as the features from the pathways found in the \pkg{Seurat} object and the features present in the \pkg{Seurat} object.
#' @param fail.safe Numeric. The maximum number of signatures the function will attempt to add to the \pkg{Seurat} object. If the number of signatures found is higher than this number, the function will not add any signature to the \pkg{Seurat} object and will only return the \pkg{Seurat} object as well as the features from the pathways found in the \pkg{Seurat} object and the features present in the \pkg{Seurat} object.
#' @param verbose Logical. If \code{FALSE}, does not print progress messages and output, but warnings and errors will still be printed.
#' @param ... Additional arguments to be passed to \code{\link[UCell]{AddModuleScore_UCell}} or \code{\link[Seurat]{AddModuleScore}}.
#'
#' @return A list containing the \pkg{Seurat} object with the added signatures if \code{only.features} = \code{FALSE}, the features from the pathways, the features present in the \pkg{Seurat} object and the names of the signatures in the \pkg{Seurat} object.
#'
#' @import Seurat
#' @import SeuratObject
#' @import msigdbr
#' @import UCell
#' @export

GSEA_Signatures = function(seurat_object,
                           assay = "RNA",
                           layer = "data",
                           species = "Homo sapiens",
                           category = NULL,
                           subcategory = NULL,
                           pathways,
                           min.features = 2,
                           signatures.names = "name",
                           method = "UCell",
                           only.features = FALSE,
                           fail.safe = 10,
                           verbose = TRUE,
                           ...) {

  if (isTRUE(verbose)) {
    cat("Starting...","\n",sep="")
  }

  if (length(Layers(seurat_object, search = "counts")) > 1) {
    if (isTRUE(verbose)) {
      cat("Joining layers...","\n")
    }
    seurat_object[[assay]] = JoinLayers(seurat_object[[assay]])
  }

  if (isTRUE(verbose)) {
    cat("Building MSigDB's pathways database with provided parameters...","\n",sep="")
  }
  msigdbgmt <- msigdbr(species = species, category = category, subcategory = subcategory)
  msigdbgmt$gs_name_source=paste0(msigdbgmt$gs_name, "_", msigdbgmt$gs_exact_source)

  if (is.character(pathways)) {
    features.from.pathways.list = list()
    found.features.list = list()
    for (i in 1:length(pathways)) {
      if (isTRUE(verbose)) {
        cat("Searching for '",pathways[i],"'...","\n",sep="")
      }
      features.from.pathways.list[[i]] = msigdbgmt[grep(pathways[i], msigdbgmt$gs_name_source),]
      if (length(features.from.pathways.list[[i]]$gs_name) == 0) {
        message("'",pathways[i],"' not found in MSigDB's pathways database, please check the spelling or the species/category/subcategory parameters",call. = FALSE, immediate. = TRUE)
        found.features.list[[i]] = features.from.pathways.list[[i]]
      }
      else {
        if (length(signatures.names) == 1) {
          if (signatures.names == "id") {
            pathways.found = unique(as.vector(features.from.pathways.list[[i]]$gs_exact_source))
          }
          else {
            pathways.found = unique(as.vector(features.from.pathways.list[[i]]$gs_name))
          }
        }
        else {
          pathways.found = unique(as.vector(features.from.pathways.list[[i]]$gs_name))
        }
        if (length(pathways.found) > 1) {
          message("'",pathways[i],"' has returned multiple results")
          features.multiple.pathways = list()
          found.multiple.pathways = list()
          for (j in 1:length(pathways.found)) {
            features.multiple.pathways[[j]] = features.from.pathways.list[[i]][grep(pathways.found[j], features.from.pathways.list[[i]]$gs_name_source),]
            features.multiple.pathways[[j]] = unique(as.vector(t(features.multiple.pathways[[j]]$human_feature_symbol)))
            found.multiple.pathways[[j]] = intersect(features.multiple.pathways[[j]],rownames(LayerData(seurat_object, assay = assay, layer = layer)))
            names(features.multiple.pathways)[j] = pathways.found[j]
            names(found.multiple.pathways)[j] = pathways.found[j]
            if (isTRUE(verbose)) {
              cat("Found ",names(features.multiple.pathways)[j],"\n",length(found.multiple.pathways[[j]])," feature(s) present / ",length(features.multiple.pathways[[j]])," total features in pathway","\n",sep="")
            }
          }
          features.from.pathways.list[[i]] = features.multiple.pathways
          found.features.list[[i]] = found.multiple.pathways
        }
        else {
          features.from.pathways.list[[i]] = unique(as.vector(t(features.from.pathways.list[[i]]$human_feature_symbol)))
          found.features.list[[i]] = intersect(features.from.pathways.list[[i]],rownames(LayerData(seurat_object, assay = assay, layer = layer)))
          names(features.from.pathways.list)[i] = pathways.found
          names(found.features.list)[i] = pathways.found
          if (isTRUE(verbose)) {
            cat("Found ",names(features.from.pathways.list)[i],"\n",length(found.features.list[[i]])," feature(s) present / ",length(features.from.pathways.list[[i]])," total features in pathway","\n",sep="")
          }
        }
      }
    }
    if (isTRUE(verbose)) {
      cat("Removing pathways with less than ",min.features," features present in Seurat object...","\n",sep="")
    }
    k = length(features.from.pathways.list)
    i = 1
    while (i <= k) {
      if (length(features.from.pathways.list[[i]]) == 0) {
        features.from.pathways.list[[i]] = NULL
        found.features.list[[i]] = NULL
        i = i - 1
      }
      if (is.list(features.from.pathways.list[[i]])) {
        features.from.pathways.list = c(features.from.pathways.list[1:(i-1)],features.from.pathways.list[[i]],features.from.pathways.list[(i+1):length(features.from.pathways.list)])
        found.features.list = c(found.features.list[1:(i-1)],found.features.list[[i]],found.features.list[(i+1):length(found.features.list)])
        features.from.pathways.list[[i]] = NULL
        found.features.list[[i]] = NULL
      }
      else {
        i = i + 1
      }
    }
    k = length(found.features.list)
    i = 1
    while (i <= k) {
      if (length(found.features.list[[i]]) < min.features) {
        found.features.list[[i]] = NULL
        k = k - 1
      }
      else {
        i = i + 1
      }
    }
    if (isTRUE(verbose)) {
      cat("Removed ",length(features.from.pathways.list) - length(found.features.list)," pathways, ",length(found.features.list)," remaining","\n",sep="")
    }
      if (length(found.features.list) == 0 & isFALSE(only.features)) {
      message("No signatures found, returning the Seurat object",call. = FALSE, immediate. = TRUE)
      return(seurat_object)
    }
    if (length(found.features.list) == 0 & isTRUE(only.features)) {
      message("No signatures found, nothing to return",call. = FALSE, immediate. = TRUE)
      return(NULL)
    }
    signatures.list = found.features.list
    if (length(signatures.names) == 1) {
      if (signatures.names == "name" | signatures.names == "id") {
        signatures.names = names(found.features.list)
      }
    }
    if (length(signatures.names) == length(found.features.list)) {
      names(signatures.list) = signatures.names
    }
    if (length(signatures.names) != length(found.features.list)) {
      stop("The length of signatures.names does not match the number of pathways found, this is likely due to the function finding multiple results from one or multiple of your queries or some of your pathways of interest has been removed due to too few features present in the Seurat object, please use only.features = TRUE and set signatures.names = 'name' or 'id' to get the number of signatures you need to provide names for")
    }
    if (only.features == TRUE) {
      if (isTRUE(verbose)) {
        cat("No signatures will be added, returning the Seurat object as well as features and pathways found.","\n",sep="")
      }
      return.list = list(seurat_object,features.from.pathways.list,found.features.list)
      names(return.list) = c("Seurat object","features from pathways","features present in Seurat object")
      return (return.list)
    }
    if (length(found.features.list) > fail.safe) {
      message("The function is attempting to add ",length(found.features.list)," signatures to the Seurat object which is higher than the fail-safe threshold, no signatures will be added, returning the Seurat object as well as features and pathways found. Increase fail.safe and run the function again to proceed with adding this many signatures",call. = FALSE, immediate. = TRUE)
      return.list = list(seurat_object,features.from.pathways.list,found.features.list)
      names(return.list) = c("Seurat object","features from pathways","features present in Seurat object")
      return (return.list)
    }
    if (length(found.features.list) == 1) {
      if (isTRUE(verbose)) {
        cat("The function will now add ",length(found.features.list)," signature to the Seurat object, this may take a while depending on the number of cells and pathways...","\n",sep="")
      }
    }
    else {
      if (isTRUE(verbose)) {
        cat("The function will now add ",length(found.features.list)," signatures to the Seurat object, this may take a while depending on the number of cells and pathways...","\n",sep="")
      }
    }
    if (isFALSE(any(Layers(seurat_object[[assay]]) %in% layer))) {
      if (isTRUE(any(Layers(seurat_object[[assay]]) %in% "data"))) {
        message("Layer '",layer,"' does not exist in the Seurat object's '",assay,"' assay, using 'data' instead",call. = FALSE, immediate. = TRUE)
        layer = "data"
      }
      else {
        message("Layer '",layer,"' does not exist in the Seurat object's '",assay,"' assay, using 'counts' instead",call. = FALSE, immediate. = TRUE)
        layer = "counts"
      }
    }
    if (method == "UCell") {
      if (isTRUE(verbose)) {
        cat("Using Module Score from UCell...","\n",sep="")
      }
      seurat_object = suppressWarnings(AddModuleScore_UCell(seurat_object, features = signatures.list, assay = assay, name = NULL, slot = layer, ...))
    }
    if (method == "Seurat") {
      if (isTRUE(verbose)) {
        cat("Using Module Score from Seurat...","\n",sep="")
      }
      seurat_object = AddModuleScore(seurat_object, features = signatures.list, assay = assay, name = "Pathway_", slot = layer, ...)
      modulescore.names = paste0("Pathway_",1:length(signatures.list))
      for (i in 1:length(modulescore.names)) {
        seurat_object[[signatures.names[i]]] = seurat_object[[modulescore.names[i]]]
        seurat_object[[modulescore.names[i]]] = NULL
      }
    }
    if (method != "UCell" & method != "Seurat") {
      stop("The method you have provided is not supported, please use 'UCell' or 'Seurat'")
    }
    if (isTRUE(verbose)) {
      cat("Done.","\n",sep="")
    }
    return.list = list(seurat_object,features.from.pathways.list,found.features.list,signatures.names)
    names(return.list) = c("Seurat object","features from pathways","features present in Seurat object","Names of signatures in Seurat object")
    return (return.list)
  }
  else {
    stop("The pathways you have provided are not in the correct format, please provide a character vector of the pathways names or ids")
  }
}
