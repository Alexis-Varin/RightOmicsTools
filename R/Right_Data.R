#' @title Data loader for vignettes
#'
#' @description This function downloads and prepares data used for vignettes, such as \code{Seurat} or \code{SingleCellExperiment} objects. The main advantage of using a data loader is that nothing is stored in the package, therefore it remains lightweight even if you use large objects in your tutorials. Anyone interested in adding a dataset to this function for use in their own package's vignettes (provided the raw data can be downloaded from a reputable source, \href{https://ftp.ncbi.nlm.nih.gov/geo/datasets/}{NCBI's GEO FTP site} for example) can open an \href{https://github.com/Alexis-Varin/RightSeuratTools/issues}{issue on GitHub}.
#'
#' @param dataset Character. The name of the dataset to load.
#' @param list.datasets Logical. If TRUE, prints the list of available datasets.
#'
#' @return An object of various classes, depending on the dataset asked, or a message containing the list of available datasets.
#'
#' @examples
#' Right_Data(list.datasets = TRUE)
#'
#' pbmc3k = Right_Data("pbmc3k")
#' pbmc3k
#' @import utils
#' @import Seurat
#' @import SeuratObject
#' @export

Right_Data = function(dataset = NULL,
                      list.datasets = FALSE) {

  if (isTRUE(list.datasets)) {
    message("Available datasets:")
    return(cat(paste0('\033[0;36m"pbmc3k"\033[0m',"   a Seurat v4 object of 2700 cells by 13714 genes","\n")))
  }

  if (dataset == "pbmc3k") {
    seurat.version = getOption("Seurat.object.assay.version")
    on.exit({
      try(suppressWarnings(suppressMessages(unlink(c("pbmc3k_filtered_gene_bc_matrices.tar.gz"), recursive = TRUE))), silent  = TRUE)
      try(suppressWarnings(suppressMessages(unlink(c("filtered_gene_bc_matrices"), recursive = TRUE))), silent  = TRUE)
      options(Seurat.object.assay.version = seurat.version)
    }, add = TRUE)
    suppressWarnings(suppressMessages(download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz", destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz", quiet = TRUE)))
    suppressWarnings(suppressMessages(untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")))
    pbmc3k.mat = suppressWarnings(suppressMessages(Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")))
    colnames(pbmc3k.mat) = suppressWarnings(suppressMessages(gsub("-1$","",colnames(pbmc3k.mat))))
    options(Seurat.object.assay.version = "v3")
    pbmc3k = suppressWarnings(suppressMessages(CreateSeuratObject(counts = pbmc3k.mat, project = "pbmc3k", min.cells = 3, min.features = 200)))
    pbmc3k = suppressWarnings(suppressMessages(NormalizeData(pbmc3k, verbose = FALSE)))
    pbmc3k@meta.data$seurat_annotations = pbmc3k.data$anno
    VariableFeatures(pbmc3k) = pbmc3k.data$hvg
    pbmc3k[["umap"]] = pbmc3k.data$umap
    Idents(pbmc3k) = "seurat_annotations"
    return(pbmc3k)
  }
}
