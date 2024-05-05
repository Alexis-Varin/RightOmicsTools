#' @title GSEA_Signatures
#'
#' @description This function creates signatures (module scores calculated from UCell or Seurat's respective functions) from the pathways's genes present in the Seurat object.
#'
#' @param seurat_object A Seurat object.
#' @param assay Character. If the Seurat object contains multiple RNA assays, you may specify which one to use (for example "RNA2" if you have created a second RNA assay you named "RNA2". See Seurat v5 vignettes for more information). You may also use another assay such as SCT to pull gene expression from.
#' @param layer Character. Formerly known as slot. If you have split layers the function will always join them before adding the signatures.
#' @param species Character. The species name to be passed to msigdbr to build the pathways database. Use msigdbr::msigdbr_species() for the list of available species.
#' @param category Character. The category or categories to be passed to msigdbr to build the pathways database. Use msigdbr::msigdbr_collections() for the list of available categories (gs_cat column). Leave NULL to use all categories.
#' @param subcategory Character. The subcategory or subcategories to be passed to msigdbr to build the pathways database. Use msigdbr::msigdbr_collections() for the list of available subcategories (gs_subcat column). Leave NULL to use all subcategories.
#' @param pathways Character. The names of the pathways to be used to create the signatures. You may provide either a pathway id (for example, "GO:0006574") or a name matching the pattern found in msigdbr$gs_name (all caps and underscores between words). Please note that you may provide a partial match (for example, "TYPE_I_INTERFERON") and the function will find all pathways containing this partial string. Beware that this may result in a large number of pathways to be added as signatures (using only.genes = TRUE is highly recommended) but is very handy to explore all pathways of interest in a particular field or biological process.
#' @param min.genes Numeric. The minimum number of genes present in the Seurat object for a pathway to be considered.
#' @param signatures.names Character. "id" will add pathways ids as signatures (e.g. GO:0004657, hsa05200 etc), "name" will add pathways names as signatures, which might be very long. You may also provide a vector of names to be used as signatures, whose length must match the number of pathways found and kept. It is recommended to use only.genes = TRUE and set signatures.names = "name" or "id" to get the number of signatures you need to provide names for.
#' @param method Character. The method you want to use to calculate the module scores, either "UCell" or "Seurat".
#' @param only.genes Logical. If TRUE, the function will not add any signature to the Seurat object and will only return the Seurat object as well as the genes from the pathways found in the Seurat object and the genes present in the Seurat object.
#' @param fail.safe Numeric. The maximum number of signatures the function will attempt to add to the Seurat object. If the number of signatures found is higher than this number, the function will not add any signature to the Seurat object and will only return the Seurat object as well as the genes from the pathways found in the Seurat object and the genes present in the Seurat object.
#' @param verbose Logical. If FALSE, does not print progress messages and output, but warnings and errors will still be printed.
#' @param ... Additional arguments to be passed to UCell's or Seurat's AddModuleScore() function.
#'
#' @return A list containing the Seurat object with the added signatures if only.genes = FALSE, the genes from the pathways, the genes present in the Seurat object and the names of the signatures in the Seurat object.
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
                           pathways = NULL,
                           min.genes = 2,
                           signatures.names = "name",
                           method = "UCell",
                           only.genes = FALSE,
                           fail.safe = 10,
                           verbose = TRUE,
                           ...) {

  if (isTRUE(verbose)) {
    cat("Starting...","\n",sep="")
  }
  if (is.null(pathways)) {
    stop("You have not provided any pathways to create signatures from")
  }

  if (isTRUE(verbose)) {
    cat("Building MSigDB's pathways database with provided parameters...","\n",sep="")
  }
  msigdbgmt <- msigdbr(species = species, category = category, subcategory = subcategory)
  msigdbgmt$gs_name_source=paste0(msigdbgmt$gs_name, "_", msigdbgmt$gs_exact_source)

  if (is.character(pathways)) {
    genes.from.pathways.list = list()
    found.genes.list = list()
    for (i in 1:length(pathways)) {
      if (isTRUE(verbose)) {
        cat("Searching for '",pathways[i],"'...","\n",sep="")
      }
      genes.from.pathways.list[[i]] = msigdbgmt[grep(pathways[i], msigdbgmt$gs_name_source),]
      if (length(genes.from.pathways.list[[i]]$gs_name) == 0) {
        warning("'",pathways[i],"' not found in MSigDB's pathways database, please check the spelling or the species/category/subcategory parameters",call. = FALSE, immediate. = TRUE)
        found.genes.list[[i]] = genes.from.pathways.list[[i]]
      }
      else {
        if (length(signatures.names) == 1) {
          if (signatures.names == "id") {
            pathways.found = unique(as.vector(genes.from.pathways.list[[i]]$gs_exact_source))
          }
          else {
            pathways.found = unique(as.vector(genes.from.pathways.list[[i]]$gs_name))
          }
        }
        else {
          pathways.found = unique(as.vector(genes.from.pathways.list[[i]]$gs_name))
        }
        if (length(pathways.found) > 1) {
          message("'",pathways[i],"' has returned multiple results")
          genes.multiple.pathways = list()
          found.multiple.pathways = list()
          for (j in 1:length(pathways.found)) {
            genes.multiple.pathways[[j]] = genes.from.pathways.list[[i]][grep(pathways.found[j], genes.from.pathways.list[[i]]$gs_name_source),]
            genes.multiple.pathways[[j]] = unique(as.vector(t(genes.multiple.pathways[[j]]$human_gene_symbol)))
            found.multiple.pathways[[j]] = intersect(genes.multiple.pathways[[j]],rownames(LayerData(seurat_object, assay = assay, layer = layer)))
            names(genes.multiple.pathways)[j] = pathways.found[j]
            names(found.multiple.pathways)[j] = pathways.found[j]
            if (isTRUE(verbose)) {
              cat("Found ",names(genes.multiple.pathways)[j],"\n",length(found.multiple.pathways[[j]])," gene(s) present / ",length(genes.multiple.pathways[[j]])," total genes in pathway","\n",sep="")
            }
          }
          genes.from.pathways.list[[i]] = genes.multiple.pathways
          found.genes.list[[i]] = found.multiple.pathways
        }
        else {
          genes.from.pathways.list[[i]] = unique(as.vector(t(genes.from.pathways.list[[i]]$human_gene_symbol)))
          found.genes.list[[i]] = intersect(genes.from.pathways.list[[i]],rownames(LayerData(seurat_object, assay = assay, layer = layer)))
          names(genes.from.pathways.list)[i] = pathways.found
          names(found.genes.list)[i] = pathways.found
          if (isTRUE(verbose)) {
            cat("Found ",names(genes.from.pathways.list)[i],"\n",length(found.genes.list[[i]])," gene(s) present / ",length(genes.from.pathways.list[[i]])," total genes in pathway","\n",sep="")
          }
        }
      }
    }
    if (isTRUE(verbose)) {
      cat("Removing pathways with less than ",min.genes," genes present in Seurat object...","\n",sep="")
    }
    k = length(genes.from.pathways.list)
    i = 1
    while (i <= k) {
      if (length(genes.from.pathways.list[[i]]) == 0) {
        genes.from.pathways.list[[i]] = NULL
        found.genes.list[[i]] = NULL
        i = i - 1
      }
      if (is.list(genes.from.pathways.list[[i]])) {
        genes.from.pathways.list = c(genes.from.pathways.list[1:(i-1)],genes.from.pathways.list[[i]],genes.from.pathways.list[(i+1):length(genes.from.pathways.list)])
        found.genes.list = c(found.genes.list[1:(i-1)],found.genes.list[[i]],found.genes.list[(i+1):length(found.genes.list)])
        genes.from.pathways.list[[i]] = NULL
        found.genes.list[[i]] = NULL
      }
      else {
        i = i + 1
      }
    }
    k = length(found.genes.list)
    i = 1
    while (i <= k) {
      if (length(found.genes.list[[i]]) < min.genes) {
        found.genes.list[[i]] = NULL
        k = k - 1
      }
      else {
        i = i + 1
      }
    }
    if (isTRUE(verbose)) {
      cat("Removed ",length(genes.from.pathways.list) - length(found.genes.list)," pathways, ",length(found.genes.list)," remaining","\n",sep="")
    }
      if (length(found.genes.list) == 0 & isFALSE(only.genes)) {
      warning("No signatures found, returning the Seurat object",call. = FALSE, immediate. = TRUE)
      return(seurat_object)
    }
    if (length(found.genes.list) == 0 & isTRUE(only.genes)) {
      warning("No signatures found, nothing to return",call. = FALSE, immediate. = TRUE)
      return(NULL)
    }
    signatures.list = found.genes.list
    if (length(signatures.names) == 1) {
      if (signatures.names == "name" | signatures.names == "id") {
        signatures.names = names(found.genes.list)
      }
    }
    if (length(signatures.names) == length(found.genes.list)) {
      names(signatures.list) = signatures.names
    }
    if (length(signatures.names) != length(found.genes.list)) {
      stop("The length of signatures.names does not match the number of pathways found, this is likely due to the function finding multiple results from one or multiple of your queries or some of your pathways of interest has been removed due to too few genes present in the Seurat object, please use only.genes = TRUE and set signatures.names = 'name' or 'id' to get the number of signatures you need to provide names for")
    }
    if (only.genes == TRUE) {
      if (isTRUE(verbose)) {
        cat("No signatures will be added, returning the Seurat object as well as genes and pathways found.","\n",sep="")
      }
      return.list = list(seurat_object,genes.from.pathways.list,found.genes.list)
      names(return.list) = c("Seurat object","Genes from pathways","Genes present in Seurat object")
      return (return.list)
    }
    if (length(found.genes.list) > fail.safe) {
      warning("The function is attempting to add ",length(found.genes.list)," signatures to the Seurat object which is higher than the fail-safe threshold, no signatures will be added, returning the Seurat object as well as genes and pathways found. Increase fail.safe and run the function again to proceed with adding this many signatures",call. = FALSE, immediate. = TRUE)
      return.list = list(seurat_object,genes.from.pathways.list,found.genes.list)
      names(return.list) = c("Seurat object","Genes from pathways","Genes present in Seurat object")
      return (return.list)
    }
    if (length(found.genes.list) == 1) {
      if (isTRUE(verbose)) {
        cat("The function will now add ",length(found.genes.list)," signature to the Seurat object, this may take a while depending on the number of cells and pathways...","\n",sep="")
      }
    }
    else {
      if (isTRUE(verbose)) {
        cat("The function will now add ",length(found.genes.list)," signatures to the Seurat object, this may take a while depending on the number of cells and pathways...","\n",sep="")
      }
    }
    if (isFALSE(any(Layers(seurat_object[[assay]]) %in% layer))) {
      if (isTRUE(any(Layers(seurat_object[[assay]]) %in% "data"))) {
        warning("Layer '",layer,"' does not exist in the Seurat object's '",assay,"' assay, using 'data' instead",call. = FALSE, immediate. = TRUE)
        layer = "data"
      }
      else {
        warning("Layer '",layer,"' does not exist in the Seurat object's '",assay,"' assay, using 'counts' instead",call. = FALSE, immediate. = TRUE)
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
    return.list = list(seurat_object,genes.from.pathways.list,found.genes.list,signatures.names)
    names(return.list) = c("Seurat object","Genes from pathways","Genes present in Seurat object","Names of signatures in Seurat object")
    return (return.list)
  }
  else {
    stop("The pathways you have provided are not in the correct format, please provide a character vector of the pathways names or ids")
  }
}
