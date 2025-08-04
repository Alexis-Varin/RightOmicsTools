#' @title Builder of CIBERSORTx Reference Matrix
#'
#' @description This function builds a Reference Matrix for \href{https://cibersortx.stanford.edu/index.php}{CIBERSORTx} from a \pkg{Seurat} object.
#'
#' @param seurat_object A \pkg{Seurat} object.
#' @param assay Character. The name of an assay containing the \code{layer} with the expression matrix. If the \code{seurat_object} contains multiple 'RNA' assays, you may specify which one to use (for example 'RNA2' if you have created a second 'RNA' assay you named 'RNA2'. See \href{https://satijalab.org/seurat/articles/seurat5_essential_commands.html#create-seurat-or-assay-objects}{Seurat v5 vignettes} for more information). You may also use another assay, such as 'SCT', if you want to extract the expression matrix for projects other than CIBERSORTx.
#' @param layer Character. The name of a layer (formerly known as slot) which stores the expression matrix. Must be 'counts' layer for CIBERSORTx. You may also specify another layer such as 'data' (log-normalized counts) or any other present in the specified \code{assay} if you want to extract the expression matrix for projects other than CIBERSORTx. If you have split layers the function will always join them before extracting the expression matrix. The raw counts matrix will be normalized to Counts Per Million (CPM) for CIBERSORTx.
#' @param ident.1 Character. The name of a metadata (for example, 'orig.ident', 'seurat_clusters', etc) from which to extract identities as column names in the Reference Matrix. If \code{NULL}, uses the active.ident metadata.
#' @param ident.2 Character. The name of a metadata from which to split \code{ident.1} identities by (for example if \code{ident.1} is 'seurat_clusters' and \code{ident.2} is 'orig.ident', identity 0 will be separated into several 0_identities, each corresponding to an orig.ident: 0_sample1, 0_sample2 etc).
#' @param double.ident Logical. If \code{TRUE}, column names of the Reference Matrix will be set with dual identities (ident.1_ident.2). If \code{FALSE}, column names will be set as \code{ident.1} identities. Ignored if \code{ident.2} = \code{NULL}.
#' @param double.ident.sep Character. The separator to use between the identity names of \code{ident.1} and \code{ident.2}. Ignored if \code{ident.2} = \code{NULL} or \code{double.ident} = \code{FALSE}.
#' @param reverse.double.ident Logical. If \code{TRUE}, the identity names of \code{ident.1} and \code{ident.2} will be reversed (ident.2_ident.1). Ignored if \code{ident.2} = \code{NULL} or \code{double.ident} = \code{FALSE}.
#' @param subset.ident.1 Character. The names of one or several \code{ident.1} identities to select. If \code{NULL}, all identities are used.
#' @param subset.ident.1.invert Logical. If \code{TRUE}, inverts the selection from \code{subset.ident.1} (for example, if \code{subset.ident.1} = c('0','1','2') and \code{subset.ident.1.invert} = \code{TRUE}, all identities except 0, 1 and 2 will be kept). Ignored if \code{subset.ident.1} = \code{NULL}.
#' @param subset.ident.2 Character. The names of one or several \code{ident.2} identities to select. If \code{NULL}, all identities are used. Ignored if \code{ident.2} = \code{NULL}.
#' @param subset.ident.2.invert Logical. If \code{TRUE}, inverts the selection from \code{subset.ident.2} (for example, if \code{subset.ident.2} = c('patient1','patient4') and \code{subset.ident.2.invert} = \code{TRUE}, all identities except patient 1 and patient 4 will be kept). Ignored if \code{subset.ident.2} = \code{NULL} or \code{ident.2} = \code{NULL}.
#' @param downsample.object.first Numeric. The number of cells to downsample from the entire \code{seurat_object} (not from each identity) before subsetting with \code{subset.ident.1} and/or \code{subset.ident.2}. Please note that the relative proportion of each identity will be kept (for example if an identity is 10% of total cells, it will remain 10% after downsampling). If \code{NULL}, all cells are used.
#' @param downsample.object.last Numeric. The number of cells to downsample from the entire \code{seurat_object} after subsetting with \code{subset.ident.1} and/or \code{subset.ident.2}. Please note that the relative proportion of each identity will be kept (for example if an identity is 10% of total cells, it will remain 10% after downsampling). If \code{NULL}, all cells are used.
#' @param downsample.cluster Numeric. The number of cells to downsample from each \code{ident.1} identity. Will be performed after \code{downsample.object.last}. If \code{NULL}, all cells are used.
#' @param automatic.downsample Logical. If \code{TRUE}, automatically downsamples the \code{seurat_object} (respecting the relative proportion of each \code{ident.1} identity, it is therefore similar to \code{downsample.object.first} and \code{downsample.object.last}) so that the Reference Matrix file written to disk would be just under the \code{max.matrix.size} limit (empirically estimated). Performed last. Ignored if \code{check.size} = \code{TRUE} or \code{write} = \code{FALSE}. Please \href{https://github.com/Alexis-Varin/RightOmicsTools/issues}{report an issue} if you see a significant difference between the file size written to disk and \code{max.matrix.size} (for example, \code{max.matrix.size} is set to 200 MB but the file is 400 MB or 100 MB).
#' @param check.size Logical. If \code{TRUE}, prints the estimated size of the Reference Matrix file that would be written to disk and the number of cells to downsample if need be.
#' @param max.matrix.size Numeric. The maximum size of the Reference Matrix file written to disk in MB. Will stop the function if the Reference Matrix file written to disk is estimated to be over this limit, or if \code{automatic.downsample} = \code{TRUE}, will downsample the \code{seurat_object} instead so that the Reference Matrix written to disk is under the size limit. Ignored if \code{write} = \code{FALSE}.
#' @param cell.barcodes Logical. Must be \code{FALSE} for CIBERSORTx. If \code{TRUE}, keeps the cell barcodes and does not rename with cell identities, if you want to extract the expression matrix for projects other than CIBERSORTx.
#' @param file.name Character. The name of the Reference Matrix file written to disk. Must not contain any space for CIBERSORTx, the function will automatically replace any space with an underscore. Ignored if \code{check.size} = \code{TRUE} or \code{write} = \code{FALSE}.
#' @param file.format Character. The format of the Reference Matrix file written to disk. Must be 'txt' or 'tsv' for CIBERSORTx, but you may also specify 'csv' for example, if you want to extract the expression matrix for projects other than CIBERSORTx. Accepts any format the \code{\link[data.table]{fwrite}} function would accept. Ignored if \code{check.size} = \code{TRUE} or \code{write} = \code{FALSE}.
#' @param file.sep Character. The separator to use in the Reference Matrix file written to disk. Must be tabulation for CIBERSORTx, but you may also specify a comma for example, if you want to extract the expression matrix for projects other than CIBERSORTx. Accepts any separator the \code{\link[data.table]{fwrite}} function would accept. Ignored if \code{check.size} = \code{TRUE} or \code{write} = \code{FALSE}.
#' @param write.path Character. The path to write the Reference Matrix into. If \code{NULL}, uses current working directory. Ignored if \code{check.size} = \code{TRUE} or \code{write} = \code{FALSE}.
#' @param write Logical. If \code{TRUE}, writes to disk the Reference Matrix file. Ignored if \code{check.size} = \code{TRUE}.
#' @param verbose Logical. If \code{FALSE}, does not print progress messages and output, but warnings and errors will still be printed.
#'
#' @return A \code{data.table} object containing the normalized counts to CPM from the \code{seurat_object} or any other specified \code{assay} and \code{layer}, with cell identities or barcodes as column names and feature names as first column. If \code{write} = \code{TRUE}, the \code{data.table} object is also written to disk. If \code{check.size} = \code{TRUE}, will instead print the estimated size of the Reference Matrix file that would be written to disk and the number of cells to downsample if need be.
#'
#' @examples
#' \dontshow{
#' suppressWarnings(suppressPackageStartupMessages(library(Seurat)))
#' }
#' # Prepare data
#' pbmc3k <- Right_Data("pbmc3k")
#'
#' # Example 1:
#' refmat <- Reference_Matrix_Builder(pbmc3k,
#'                                    ident.1 = "seurat_annotations",
#'                                    ident.2 = "orig.ident",
#'                                    downsample.object.first = 1500,
#'                                    subset.ident.1 = c("Naive CD4 T",
#'                                                       "Memory CD4 T",
#'                                                       "CD8 T",
#'                                                       "NK"),
#'                                    subset.ident.2 = "Donor_2",
#'                                    subset.ident.2.invert = TRUE,
#'                                    write = FALSE)
#' refmat[1:5, 1:5]
#'
#' # Example 2:
#' Reference_Matrix_Builder(pbmc3k,
#'                          check.size = TRUE,
#'                          max.matrix.size = 40)
#'
#' @import Seurat
#' @import SeuratObject
#' @import data.table
#' @export

Reference_Matrix_Builder = function(
    seurat_object,
    assay = "RNA",
    layer = "counts",
    ident.1 = NULL,
    ident.2 = NULL,
    double.ident = TRUE,
    double.ident.sep = "_",
    reverse.double.ident = FALSE,
    subset.ident.1 = NULL,
    subset.ident.1.invert = FALSE,
    subset.ident.2 = NULL,
    subset.ident.2.invert = FALSE,
    downsample.object.first = NULL,
    downsample.object.last = NULL,
    downsample.cluster = NULL,
    automatic.downsample = FALSE,
    check.size = FALSE,
    max.matrix.size = 500,
    cell.barcodes = FALSE,
    file.name = "Reference_Matrix",
    file.format = "txt",
    file.sep = "\t",
    write.path = NULL,
    write = TRUE,
    verbose = TRUE) {

  if (isTRUE(verbose)) {
    cat("Starting...","\n")
  }

  if (is.null(ident.2) & is.character(subset.ident.2)) {
    stop("You must provide an ident.2 metadata name to subset subset.ident.2 identities from")
  }

  to.put.back = Idents(seurat_object)

  if (is.character(ident.1)) {
    Idents(seurat_object) = ident.1
  }

  seurat_object@meta.data$ident.1 = Idents(seurat_object)

  if (is.character(ident.2)) {
    checkident.1 = Idents(seurat_object)
    Idents(seurat_object) = ident.2
    seurat_object@meta.data$ident.2 = Idents(seurat_object)
    checkident.2 = Idents(seurat_object)
    if (checkident.1[1] == checkident.2[1]) {
      stop("ident.1 and ident.2 or seurat_object@active.ident and ident.2 are the same, please choose different identities")
    }
    if (isTRUE(reverse.double.ident)) {
      seurat_object@meta.data$double.ident = paste(checkident.2,checkident.1,sep = double.ident.sep)
    }
    else {
      seurat_object@meta.data$double.ident = paste(checkident.1,checkident.2,sep = double.ident.sep)
    }
    Idents(seurat_object) = "double.ident"
  }

  gc(verbose = FALSE)

  if(is.numeric(downsample.object.first)) {
    if (isTRUE(verbose)) {
      cat("Downsampling...","\n")
    }
    if (downsample.object.first > ncol(seurat_object)) {
      warning("downsample.object.first is greater than the number of cells in seurat_object, no downsampling was done",immediate. = T)
    }
    props = as.data.frame(table(Idents(seurat_object)))
    cell.list = lapply(unique(Idents(seurat_object)), function(celltype) {
      WhichCells(seurat_object, idents = celltype, downsample = ceiling(downsample.object.first*props[props$Var1 == celltype,"Freq"]/ncol(seurat_object)))
    })
    seurat_object = seurat_object[,unlist(cell.list)]
  }

  if(is.character(subset.ident.1)) {
    if (isTRUE(verbose)) {
      cat("Subsetting...","\n")
    }
    Idents(seurat_object) = "ident.1"
    seurat_object = subset(seurat_object, idents = subset.ident.1, invert = subset.ident.1.invert)
    if (is.character(ident.2)) {
      Idents(seurat_object) = "double.ident"
    }
  }

  if (is.character(subset.ident.2)) {
    if (isTRUE(verbose)) {
      cat("Subsetting...","\n")
    }
    Idents(seurat_object) = "ident.2"
    seurat_object = subset(seurat_object, idents = subset.ident.2, invert = subset.ident.2.invert)
    Idents(seurat_object) = "double.ident"
  }

  if(is.numeric(downsample.object.last)) {
    if (isTRUE(verbose)) {
      cat("Downsampling...","\n")
    }
    if (downsample.object.last > ncol(seurat_object)) {
      warning("downsample.object.last is greater than the number of cells in seurat_object, no downsampling was done",immediate. = T)
    }
    props = as.data.frame(table(Idents(seurat_object)))
    cell.list = lapply(unique(Idents(seurat_object)), function(celltype) {
      WhichCells(seurat_object, idents = celltype, downsample = ceiling(downsample.object.last*props[props$Var1 == celltype,"Freq"]/ncol(seurat_object)))
    })
    seurat_object = seurat_object[,unlist(cell.list)]
  }

  if(is.numeric(downsample.cluster)) {
    if (isTRUE(verbose)) {
      cat("Downsampling...","\n")
    }
    Idents(seurat_object) = "ident.1"
    cell.list = WhichCells(seurat_object, downsample = downsample.cluster)
    seurat_object2 = seurat_object[,cell.list]
    if (length(colnames(seurat_object2)) == length(colnames(seurat_object))) {
      warning("downsample.cluster is greater than the number of cells in each identity, no downsampling was done",immediate. = T)
    }
    seurat_object = seurat_object2
    if (is.character(ident.2)) {
      Idents(seurat_object) = "double.ident"
    }
  }

  if (isFALSE(double.ident)) {
    Idents(seurat_object) = "ident.1"
  }

   if (length(Layers(seurat_object, search = layer)) > 1) {
     if (isTRUE(verbose)) {
      cat("Joining layers...","\n")
    }
    seurat_object2 = seurat_object
    seurat_object2[[assay]] = JoinLayers(seurat_object2[[assay]])
  }
  else {
    seurat_object2 = seurat_object
  }
  if (layer == "counts" & assay == "RNA") {
    if (isTRUE(verbose)) {
      cat("Normalizing and extracting the expression matrix...","\n")
    }
    refmat = NormalizeData(seurat_object2, assay = assay, normalization.method = "RC", scale.factor = 1e6, verbose = FALSE)
    refmat = LayerData(seurat_object2, assay = assay, layer = "data")
  }
  else {
    if (isTRUE(verbose)) {
      cat("Extracting the expression matrix...","\n")
    }
    refmat = LayerData(seurat_object2, assay = assay, layer = layer)
  }

  projected.cell.number = trunc(425000/1.02*max.matrix.size/length(rownames(refmat)))

  if (isTRUE(check.size)) {
    refmat.size = length(refmat)/425000/1.02

    gc(verbose = FALSE)

    if (projected.cell.number >= length(colnames(refmat))) {
      return(cat("Current estimated Reference Matrix size on CIBERSORTx web portal is between ",
          trunc(refmat.size/1.01),
          " and ",
          trunc(refmat.size*1.02),
          " MB :",
          "\n",
          "Matrix of ",
          length(colnames(refmat)),
          " cells by ",
          length(rownames(refmat)),
          " features",
          "\n",
          "You do not need to downsample your Seurat object","\n",sep=""))
    }

    else {
      return(cat("Current estimated Reference Matrix size on CIBERSORTx web portal is between ",
          trunc(refmat.size/1.01),
          " and ",
          trunc(refmat.size*1.02),
          " MB :",
          "\n",
          "Matrix of ",
          length(colnames(refmat)),
          " cells by ",
          length(rownames(refmat)),
          " features",
          "\n",
          "You may want to downsample your Seurat object to ",
          projected.cell.number,
          " cells for a Reference Matrix under ",
          max.matrix.size,
          " MB","\n",sep=""))
    }
  }

  if (length(refmat) > 425000*max.matrix.size/1.02 & isTRUE(write)) {
    if (isFALSE(automatic.downsample)) {
      stop(paste0("The Reference Matrix file is projected to be over the size limit of ",
               max.matrix.size,
               " MB on CIBERSORTx web portal :",
               "\n",
               " Matrix of ",
               ncol(refmat),
               " cells by ",
               nrow(refmat),
               " features",
               "\n",
               " Please subset identities, downsample the number of cells or set automatic.downsample = TRUE"))
    }
    else {
      warning(paste0("The Reference Matrix file is projected to be over the size limit of ",
      max.matrix.size,
      " MB on CIBERSORTx web portal :",
      "\n",
      " Matrix of ",
      ncol(refmat),
      " cells by ",
      nrow(refmat),
      " features",
      "\n",
      " Automatically downsampling to ",
      projected.cell.number,
      " cells to be under the size limit..."),immediate. = T)
      props = as.data.frame(table(Idents(seurat_object)))
      cell.list = lapply(unique(Idents(seurat_object)), function(celltype) {
        WhichCells(seurat_object, idents = celltype, downsample = ceiling(projected.cell.number*props[props$Var1 == celltype,"Freq"]/ncol(seurat_object)))
      })
      seurat_object = seurat_object[,unlist(cell.list)]
      if (length(Layers(seurat_object, search = layer)) > 1) {
        if (isTRUE(verbose)) {
          cat("Joining layers...","\n")
        }
        seurat_object[[assay]] = JoinLayers(seurat_object[[assay]])
      }
      if (layer == "counts" & assay == "RNA") {
        refmat = NormalizeData(seurat_object, assay = assay, normalization.method = "RC", scale.factor = 1e6, verbose = FALSE)
        refmat = LayerData(seurat_object, assay = assay, layer = "data")
      }
      else {
        refmat = LayerData(seurat_object, assay = assay, layer = layer)
      }
    }
  }

  if (isTRUE(verbose)) {
    cat("Building the data.table...","\n")
  }
  refmat = as.data.frame(refmat)
  refmat = cbind(Gene = rownames(refmat),refmat)
  refmat = setDT(refmat)
  if (isFALSE(cell.barcodes)) {
    colnames(refmat) = c("Gene", as.character(seurat_object@active.ident))
  }

  if (isTRUE(write)) {
    if(is.null(write.path)) {
      write.path = getwd()
    }
    if (isTRUE(grepl(" ",file.name))) {
      warning("The file name contains one or several spaces, renaming with underscores as CIBERSORTx will otherwise report an error with the Reference Matrix...",immediate. = T)
      file.name = gsub(" ","_",file.name)
    }
    if (isTRUE(verbose)) {
      cat("Writing to disk...","\n")
    }
    fwrite(refmat, file = paste0(write.path,"/",file.name,".",file.format), sep = file.sep, quote = F)
  }

  if (isTRUE(verbose)) {
    cat("Cleaning...","\n")
  }
  gc(verbose = FALSE)
  if (isTRUE(verbose)) {
    cat("Done.","\n")
  }
  return(refmat)
}
