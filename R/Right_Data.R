#' @title Data loader for vignettes
#'
#' @description This function downloads and prepares data used for vignettes. The main advantage of using a data loader is that nothing is stored in the package, therefore it remains lightweight even if you use large objects in your tutorials. Anyone interested in adding a dataset to this function for use in their own package's vignettes (provided the raw data can be downloaded from a reputable source, \href{https://ftp.ncbi.nlm.nih.gov/geo/datasets/}{NCBI's GEO FTP site} for example) can open an \href{https://github.com/Alexis-Varin/RightOmicsTools/issues}{issue on GitHub}.
#'
#' @param dataset Character. The name of the dataset to load.
#' @param list.datasets Logical. If \code{TRUE}, prints the names and descriptions of available datasets.
#'
#' @return An object of various classes, depending on the dataset, or a message containing the names and descriptions of available datasets.
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
    return(cat(paste0('\033[0;36m"pbmc3k"\033[0m',"   a Seurat v4 object of 2,700 cells by 13,714 genes\n",
                      "           PBMCs of a healthy donor from 10XGenomics\n\n",
                      '\033[0;36m"monolps"\033[0m',"   a Seurat v5 object of 3,693 cells by 23,798 genes\n",
                      "            Monocytes stimulated or not with LPS from GSE226488\n")))
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
    pbmc3k = suppressWarnings(suppressMessages(CreateSeuratObject(counts = pbmc3k.mat, project = "pbmc3k", min.cells = 3, min.features = 200, meta.data = as.data.frame(all.data[["pbmc3k"]][1:3]))))
    pbmc3k = suppressWarnings(suppressMessages(NormalizeData(pbmc3k, verbose = FALSE)))
    VariableFeatures(pbmc3k) = all.data[["pbmc3k"]]$hvg
    pbmc3k[["umap"]] = all.data[["pbmc3k"]]$umap
    Idents(pbmc3k) = "seurat_annotations"
    return(pbmc3k)
  }

  if (dataset == "monolps") {
    on.exit({
      try(suppressWarnings(suppressMessages(unlink(c("pbmc_resting.tar.gz"), recursive = TRUE))), silent  = TRUE)
      try(suppressWarnings(suppressMessages(unlink(c("pbmc_stimulated.tar.gz"), recursive = TRUE))), silent  = TRUE)
      try(suppressWarnings(suppressMessages(unlink(c("pbmc_resting"), recursive = TRUE))), silent  = TRUE)
      try(suppressWarnings(suppressMessages(unlink(c("pbmc_stimulated"), recursive = TRUE))), silent  = TRUE)
    }, add = TRUE)
    suppressWarnings(suppressMessages(download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7077nnn/GSM7077865/suppl/GSM7077865%5FD1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.tar.gz", destfile = "pbmc_resting.tar.gz", quiet = TRUE)))
    suppressWarnings(suppressMessages(untar("pbmc_resting.tar.gz", exdir = "pbmc_resting")))
    pbmc.resting.mat = suppressWarnings(suppressMessages(Read10X(data.dir = "pbmc_resting")))
    pbmc.resting = suppressWarnings(suppressMessages(CreateSeuratObject(counts = pbmc.resting.mat, project = "resting", min.cells = 3, min.features = 200)))
    suppressWarnings(suppressMessages(download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7077nnn/GSM7077866/suppl/GSM7077866%5FG1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.tar.gz", destfile = "pbmc_stimulated.tar.gz", quiet = TRUE)))
    suppressWarnings(suppressMessages(untar("pbmc_stimulated.tar.gz", exdir = "pbmc_stimulated")))
    pbmc.stimulated.mat = suppressWarnings(suppressMessages(Read10X(data.dir = "pbmc_stimulated")))
    pbmc.stimulated = suppressWarnings(suppressMessages(CreateSeuratObject(counts = pbmc.stimulated.mat, project = "stimulated", min.cells = 3, min.features = 200)))
    pbmc.all = suppressWarnings(suppressMessages(merge(pbmc.resting, pbmc.stimulated)))
    pbmc.all[["RNA"]] = suppressWarnings(suppressMessages(JoinLayers(pbmc.all[["RNA"]])))
    pbmc.all@project.name = "Monocytes resting or LPS stimulated"
    mono = pbmc.all[ , all.data[["mono"]]$cells]
    mono = suppressWarnings(suppressMessages(NormalizeData(mono, verbose = FALSE)))
    VariableFeatures(mono) = all.data[["mono"]]$hvg
    mono[["umap"]] = all.data[["mono"]]$umap
    mono@meta.data$seurat_annotations = all.data[["mono"]]$seurat_annotations
    Idents(mono) = "seurat_annotations"
    return(mono)
  }
}
