#' @title Builder of CIBERSORTx Mixture File
#'
#' @description This function builds a Mixture File for \href{https://cibersortx.stanford.edu/index.php}{CIBERSORTx} from bulk RNA-seq and/or microarray files and/or R objects.
#'
#' @param objects A \code{data.frame}, \code{data.table} or \code{matrix} object or a mixed \code{list} of \code{data.frame}, \code{data.table} and/or \code{matrix} objects corresponding to bulk RNA-seq and/or microarray experiments, with genes as rows and gene names as row names or first column. A single object may contain multiple samples and/or experiments as columns as long as it contains a single gene names column as row names or first column (i.e. | gene | sample1 | sample2 | sample3 | etc). You may mix objects with different gene names order and number (for example, an object containing 13567 genes and another with 18458 genes).
#' @param files.path Character. The path to read txt, csv and/or tsv files corresponding to bulk RNA-seq or microarray experiments, with genes as rows and gene names as first column. A single file may contain multiple samples and/or experiments as columns as long as it contains a single gene names column as first column. You may mix files with different genes names order and number (for example, a file containing 21848 genes and another with 19457 genes).
#' @param file.name Character. The name of the Mixture File written to disk. Must not contain any space for CIBERSORTx, the function will automatically replace any space with an underscore. Ignored if \code{write} = \code{FALSE}.
#' @param file.format Character. The format of the Mixture File written to disk. Must be 'txt' or 'tsv' for CIBERSORTx, but you may also specify 'csv' for example, for projects other than CIBERSORTx. Accepts any format the \code{\link[data.table]{fwrite}} function would accept. Ignored if \code{write} = \code{FALSE}.
#' @param file.sep Character. The separator to use in the Mixture File written to disk. Must be tabulation for CIBERSORTx, but you may also specify a comma for example, for projects other than CIBERSORTx. Accepts any separator the \code{\link[data.table]{fwrite}} function would accept. Ignored if \code{write} = \code{FALSE}.
#' @param write.path Character. The path to write the Mixture File into. If \code{NULL}, uses current working directory. Ignored if \code{write} = \code{FALSE}.
#' @param write Logical. If \code{TRUE}, writes to disk the Mixture File.
#' @param verbose Logical. If \code{FALSE}, does not print progress messages and output, but warnings and errors will still be printed.
#'
#' @return A \code{data.table} object containing the counts of the objects and/or files provided, with each column being a sample. If \code{write} = \code{TRUE}, the \code{data.table} object is also written to disk.
#'
#' @import data.table
#' @export

Mixture_File_Builder = function(
    objects = NULL,
    files.path = NULL,
    file.name = "Mixture_File",
    file.format = "txt",
    file.sep = "\t",
    write.path = NULL,
    write = TRUE,
    verbose = TRUE) {

  if (is.null(objects) & is.null(files.path))
    stop("Please provide at least one object or a path to your files")

  if (is.character(files.path)) {
    if (isTRUE(verbose))
      cat("Reading files from ",files.path,"...","\n",sep="")
    dt.list = lapply(list.files(path=files.path, pattern="*.txt|*.csv|*.tsv", full.names = T),
                     function(x) fread(file = x, header = T))
  }

  if (!is.null(objects)) {
    if (isTRUE(verbose))
      cat("Adding objects...","\n")

    if (!inherits(class(objects), "list")) objects = list(objects)
    objects = lapply(objects, function(x) {
      if (is.matrix(x)) x = as.data.frame(x)
      x = setDT(x)
      })
    if (is.character(files.path))
      dt.list = c(dt.list,objects)
    else
      dt.list = objects
  }

  if (isTRUE(verbose))
    cat("Building the data.table...","\n")
  dt.list = lapply(dt.list, function(x) {names(x)[1] = "Gene"; x})
  merged.dt = Reduce(function(x, y) merge(x, y, all = TRUE), dt.list)
  merged.dt[is.na(merged.dt)] = 0

  if (isTRUE(write)) {
    if(is.null(write.path))
      write.path = getwd()
    if (isTRUE(grepl(" ",file.name))) {
      warning("The file name contains one or several spaces, renaming with underscores as CIBERSORTx will otherwise report an error...",immediate. = T)
      file.name = gsub(" ","_",file.name)
    }
    if (isTRUE(verbose))
      cat("Writing to disk...","\n")
    fwrite(x=merged.dt, file = paste0(write.path,"/",file.name,".",file.format), sep = file.sep, quote = F)
  }

  if (isTRUE(verbose))
    cat("Cleaning...","\n")
  gc(verbose = FALSE)
  if (isTRUE(verbose))
    cat("Done.","\n")
  return(merged.dt)
}
