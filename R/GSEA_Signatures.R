#' @title Gene signatures from GSEA-MSigDB pathways
#'
#' @description Following gene set enrichment analysis (GSEA), one often would like to explore the expression of genes comprised in enriched pathways. For this purpose, this function builds a pathway database from \href{https://www.gsea-msigdb.org/gsea/msigdb}{MSigDB} and creates signatures (module scores calculated from \pkg{UCell} or \pkg{Seurat}'s respective functions) from features found in a \pkg{Seurat} object and extracted from supplied pathways. It also returns the feature names which can be used to visualize their individual expression, using for example \code{\link[RightOmicsTools]{DotPlot_Heatmap}}.
#'
#' @param seurat_object A \pkg{Seurat} object.
#' @param assay Character. The name of an assay containing the \code{layer} with the expression matrix. If the \code{seurat_object} contains multiple 'RNA' assays, you may specify which one to use (for example, 'RNA2' if you have created a second 'RNA' assay you named 'RNA2'. See \href{https://satijalab.org/seurat/articles/seurat5_essential_commands.html#create-seurat-or-assay-objects}{Seurat v5 vignettes} for more information). You may also use another assay, such as 'SCT', to pull feature expression from.
#' @param layer Character. The name of a layer (formerly known as slot) which stores the expression matrix. If the \code{seurat_object} contains split layers, the function will always join them before searching features and adding the signatures.
#' @param species Character. The species name to be internally passed to \code{\link[msigdbr]{msigdbr}} to build the pathway database. Use \pkg{msigdbr}::\code{\link[msigdbr]{msigdbr_species}} for the names of available species.
#' @param category Character. The names of one or several categories to be internally passed to \code{\link[msigdbr]{msigdbr}} to build the pathway database. Use \pkg{msigdbr}::\code{\link[msigdbr]{msigdbr_collections}} for the names of available categories (gs_cat column). If \code{NULL}, all categories will be used.
#' @param subcategory Character. The names of one or several subcategories to be internally passed to \code{\link[msigdbr]{msigdbr}} to build the pathway database. Use \pkg{msigdbr}::\code{\link[msigdbr]{msigdbr_collections}} for the names of available subcategories (gs_subcat column). If \code{NULL}, all subcategories will be used.
#' @param pathways Character. The names of one or several pathways to be searched in the pathway database and added as signatures. You may provide either a pathway id (for example, 'GO:0006574') or a name matching the pattern found in the gs_name column (uppercase letters and underscores between words). Please note that you may also provide a partial match (for example, 'TYPE_I_INTERFERON') and the function will find all pathways containing this partial pattern. Beware that this may result in a large number of pathways to be added as signatures (using \code{only.features} = \code{TRUE} is highly recommended) but is very handy to explore all pathways of interest in a particular biological process.
#' @param min.features Numeric. The minimum number of features present in the \code{seurat_object} for a pathway to be added as a signature.
#' @param signatures.names Character. Either 'id', which will add the ids of the \code{pathways} as signature names (for example, 'GO:0004657', 'hsa05200' etc), or 'name', which will add the names of the \code{pathways} as signature names (for example, 'GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY'). You may also provide custom names to be used as signature names, whose length must match the length of \code{pathways} supplied. If multiple results are found for a pathway, the function will append a number to the corresponding custom signature name for each result.
#' @param method Character. The method used to calculate the module scores, either 'UCell' or 'Seurat'.
#' @param only.features Logical. If \code{TRUE}, the function will not add any signature to the \code{seurat_object} and will only return the feature names from the \code{pathways} found in the \code{seurat_object} and the feature names present in the \code{seurat_object}.
#' @param fail.safe Numeric. The maximum number of signatures the function will attempt to add to the \code{seurat_object}. If the number of signatures found is higher than this number, the function will not add any signature, and will instead return the \code{seurat_object} as well as the feature names from the \code{pathways} found in the \code{seurat_object} and the feature names present in the \code{seurat_object}. This prevents the function from adding a large number of signatures to the \code{seurat_object} by mistake.
#' @param verbose Logical. If \code{FALSE}, does not print progress messages and output, but warnings and errors will still be printed.
#' @param ... Additional arguments to be passed to \code{\link[UCell]{AddModuleScore_UCell}} or \code{\link[Seurat]{AddModuleScore}}, such as \code{nbin} or \code{maxRank}.
#'
#' @return A \code{list} containing the \code{seurat_object} with added signatures, all feature names from the \code{pathways} found in the \code{seurat_object}, the feature names present in the \code{seurat_object} and the signature names. If \code{only.features} = \code{TRUE}, the function will instead return a \code{list} containing the feature names from the \code{pathways} found in the \code{seurat_object} and the feature names present in the \code{seurat_object}.
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

  if (length(signatures.names) != length(pathways)) {
    if (signatures.names != "id" & signatures.names != "name") {
    stop("signatures.names needs to be a character vector of same length as the number of specified pathways")
    }
  }

  if (isFALSE(any(Assays(seurat_object) %in% assay))) {
    message("Assay '",assay,"' was not found in the Seurat object, using 'RNA' instead")
    assay = "RNA"
  }

  if (length(Layers(seurat_object, search = "counts")) > 1) {
    if (isTRUE(verbose)) {
      cat("Joining layers...","\n")
    }
    seurat_object[[assay]] = JoinLayers(seurat_object[[assay]])
  }

  if (isFALSE(any(Layers(seurat_object[[assay]]) %in% layer))) {
    if (isTRUE(any(Layers(seurat_object[[assay]]) %in% "data"))) {
      message("Layer '",layer,"' was not found in the Seurat object's '",assay,"' assay, using 'data' instead")
      layer = "data"
    }
    else {
      message("Layer '",layer,"' was not found in the Seurat object's '",assay,"' assay, using 'counts' instead")
      layer = "counts"
    }
  }

  if (isTRUE(verbose)) {
    cat("Building MSigDB's pathway database with provided parameters...","\n",sep="")
  }
  msigdbgmt <- msigdbr(species = species, category = category, subcategory = subcategory)
  msigdbgmt$gs_name_source=paste0(msigdbgmt$gs_name, "_", msigdbgmt$gs_exact_source)

  if (is.character(pathways)) {
    features.from.pathways.list = found.features.list = list()
    pathways.names = NULL
    for (i in 1:length(pathways)) {
      if (isTRUE(verbose)) {
        message("Searching for '",pathways[i],"'...")
      }
      features.from.pathways.list[[i]] = msigdbgmt[grep(pathways[i], msigdbgmt$gs_name_source),]
      if (length(features.from.pathways.list[[i]]$gs_name) == 0) {
        message("'",pathways[i],"' not found in MSigDB's pathway database")
        found.features.list[[i]] = features.from.pathways.list[[i]] = NULL
      }
      else {
        if (length(signatures.names) == 1) {
          if (signatures.names == "id") {
            pathways.found = unique(as.vector(features.from.pathways.list[[i]]$gs_exact_source))
            sig.name = pathways.found
          }
          else if (signatures.names == "name") {
            pathways.found = unique(as.vector(features.from.pathways.list[[i]]$gs_name))
            sig.name = pathways.found
          }
          else {
            pathways.found = unique(as.vector(features.from.pathways.list[[i]]$gs_name))
            sig.name = signatures.names
          }
        }
        else {
          pathways.found = unique(as.vector(features.from.pathways.list[[i]]$gs_name))
          sig.name = signatures.names[i]
        }
        if (length(pathways.found) > 1) {
          if (isTRUE(verbose)) {
            message("'",pathways[i],"' has returned multiple results")
          }
          if (length(signatures.names) > 1) {
            pathways.names = c(pathways.names,paste0(sig.name,"_",1:length(pathways.found)))
          }
          else {
            pathways.names = sig.name
          }
          features.multiple.pathways = found.multiple.pathways = list()
          for (j in 1:length(pathways.found)) {
            features.multiple.pathways[[j]] = features.from.pathways.list[[i]][grep(pathways.found[j], features.from.pathways.list[[i]]$gs_name_source),]
            features.multiple.pathways[[j]] = unique(as.vector(t(features.multiple.pathways[[j]]$human_gene_symbol)))
            found.multiple.pathways[[j]] = intersect(features.multiple.pathways[[j]],rownames(LayerData(seurat_object, assay = assay, layer = layer)))
            names(features.multiple.pathways)[j] = pathways.found[j]
            names(found.multiple.pathways)[j] = pathways.found[j]
            if (isTRUE(verbose)) {
              cat("Found ",pathways.found[j],"\n",length(found.multiple.pathways[[j]])," feature(s) present / ",length(features.multiple.pathways[[j]])," total features in pathway","\n",sep="")
            }
          }
          features.from.pathways.list[[i]] = features.multiple.pathways
          found.features.list[[i]] = found.multiple.pathways
        }
        else {
          pathways.names = c(pathways.names,sig.name)
          features.from.pathways.list[[i]] = list(unique(as.vector(t(features.from.pathways.list[[i]]$human_gene_symbol))))
          found.features.list[[i]] = list(intersect(features.from.pathways.list[[i]][[1]],rownames(LayerData(seurat_object, assay = assay, layer = layer))))
          names(features.from.pathways.list[[i]]) = pathways.found
          names(found.features.list[[i]]) = pathways.found
          if (isTRUE(verbose)) {
            cat("Found ",pathways.found,"\n",length(found.features.list[[i]][[1]])," feature(s) present / ",length(features.from.pathways.list[[i]][[1]])," total features in pathway","\n",sep="")
          }
        }
      }
    }
    tmp = tmp2 = list()
    for (i in 1:length(features.from.pathways.list)) {
      tmp = c(tmp, features.from.pathways.list[[i]])
      tmp2 = c(tmp2, found.features.list[[i]])
    }
    features.from.pathways.list = tmp
    found.features.list = tmp2
    features.from.pathways.list = features.from.pathways.list[lengths(features.from.pathways.list) != 0]
    found.features.list = found.features.list[lengths(found.features.list) != 0]
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
    if (isTRUE(verbose) & (length(features.from.pathways.list) - length(found.features.list)) > 0) {
      message("Removed ",length(features.from.pathways.list) - length(found.features.list)," pathway(s) with less than ",min.features," features present in Seurat object, ",length(found.features.list)," remaining","\n",sep="")
    }
    features.from.pathways.list = features.from.pathways.list[lengths(features.from.pathways.list) != 0]
    found.features.list = found.features.list[lengths(found.features.list) != 0]
    if (length(found.features.list) == 0) {
      message("No pathways found, returning the Seurat object")
      return(seurat_object)
    }
    signatures.list = found.features.list
    names(signatures.list) = pathways.names
    if (only.features == TRUE) {
      if (isTRUE(verbose)) {
        cat("No signatures will be added, returning features and pathways found","\n",sep="")
      }
      return.list = list(features.from.pathways.list,found.features.list)
      names(return.list) = c("Features from pathways","Features present in the Seurat object")
      return (return.list)
    }
    if (length(found.features.list) > fail.safe) {
      message("The function is attempting to add ",length(found.features.list)," signatures to the Seurat object which is higher than the fail-safe threshold.\nNo signatures will be added, returning the Seurat object as well as features and pathways found.\nIncrease fail.safe and run the function again to proceed with adding this many signatures")
      return.list = list(seurat_object,features.from.pathways.list,found.features.list)
      names(return.list) = c("Seurat object","Features from pathways","Features present in the Seurat object")
      return (return.list)
    }
    if (isTRUE(verbose)) {
      cat("The function will now add ",length(found.features.list)," signature(s) to the Seurat object\nThis may take a while depending on the number of cells and pathways...","\n",sep="")
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
        seurat_object[[pathways.names[i]]] = seurat_object[[modulescore.names[i]]]
        seurat_object[[modulescore.names[i]]] = NULL
      }
    }
    if (method != "UCell" & method != "Seurat") {
      stop("The method you have provided is not supported, please use 'UCell' or 'Seurat'")
    }
    if (isTRUE(verbose)) {
      cat("Done.","\n",sep="")
    }
    return.list = list(seurat_object,features.from.pathways.list,found.features.list,pathways.names)
    names(return.list) = c("Seurat object","Features from pathways","Features present in the Seurat object","Signature names in the Seurat object")
    return (return.list)
  }
  else {
    stop("Please provide a character vector of the pathways names or ids")
  }
}
