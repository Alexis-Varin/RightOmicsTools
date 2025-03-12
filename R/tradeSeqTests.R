#' @title Perform \pkg{tradeSeq} tests in parallel
#'
#' @description This function is a wrapper around \pkg{tradeSeq} test functions, such as \code{\link[tradeSeq]{associationTest}} or \code{\link[tradeSeq]{conditionTest}}, and allows multiple tests to be performed simultaneously, with all necessary parameters specified internally. It also organizes the output of each test into separate \code{list} elements, making it easy to access results. Finally, it calculates an adjusted p-value and orders based on decreasing Waldstat.
#'
#' @param models A \pkg{SingleCellExperiment} object containing the fitted GAM smoothers, computed using \code{\link[tradeSeq]{fitGAM}}, with or without conditions provided.
#' @param tests Character. The names of one or several tests to perform. You may provide partial names, such as 'assoc' for \code{\link[tradeSeq]{associationTest}}, as long as a single match is found. The available tests are \code{\link[tradeSeq]{associationTest}}, \code{\link[tradeSeq]{patternTest}}, \code{\link[tradeSeq]{diffEndTest}}, \code{\link[tradeSeq]{earlyDETest}}, \code{\link[tradeSeq]{startVsEndTest}} and \code{\link[tradeSeq]{conditionTest}}. \code{\link[tradeSeq]{conditionTest}} is ignored if \code{models} was computed using \code{\link[tradeSeq]{fitGAM}} without conditions.
#' @param global Logical. (from \code{\link[tradeSeq]{earlyDETest}}/\code{\link[tradeSeq]{startVsEndTest}} documentation) If \code{TRUE}, test for all pairwise comparisons/lineages simultaneously.
#' @param pairwise Logical. (from \code{\link[tradeSeq]{earlyDETest}} documentation) If \code{TRUE}, test for all pairwise comparisons independently. Ignored in tests that do not have \code{pairwise}.
#' @param lineages Logical. (from \code{\link[tradeSeq]{startVsEndTest}} documentation) If \code{TRUE}, test for all lineages independently. Ignored in tests that do not have \code{lineages}.
#' @param l2fc Numeric. (from \code{\link[tradeSeq]{earlyDETest}} documentation) The log2 fold change threshold to test against. Note, that this will affect both the global test and the pairwise comparisons.
#' @param eigen.thresh Numeric. (from \code{\link[tradeSeq]{earlyDETest}} documentation) Eigenvalue threshold for inverting the variance-covariance matrix of the coefficients to use for calculating the Wald test statistics. Lower values are more lenient to adding more information but also decrease computational stability. This argument should in general not be changed by the user but is provided for back-compatability. Set to 1e-8 to reproduce results of older version of 'tradeSeq'. Ignored in tests that do not have \code{eigen.thresh}.
#' @param n.points Numeric. (from \code{\link[tradeSeq]{earlyDETest}} documentation) The number of points to be compared between lineages. Defaults to twice the number of knots. Ignored in tests that do not have \code{n.points}.
#' @param knots Numeric or List. (from \code{\link[tradeSeq]{earlyDETest}} documentation) A vector of length 2 specifying the knots at the start and end of the region of interest. You may also provide a \code{list} of multiple elements (for example, list(c(2,4), c(3,4))) to repeat the test for multiple sets of knots. Ignored in tests that do not have \code{knots}.
#' @param pseudotime.values Numeric or List. (from \code{\link[tradeSeq]{startVsEndTest}} documentation) A vector of length 2, specifying two pseudotime values to be compared against each other, for every lineage of the trajectory. @details Note that this test assumes that all lineages start at a pseudotime value of zero, which is the starting point against which the end point is compared. You may also provide a \code{list} of multiple elements (for example, list(c(8,12), c(6,14))) to repeat the test for multiple sets of pseudotime values. Ignored in tests that do not have \code{pseudotime.values}.
#' @param parallelized Logical. If \code{TRUE}, the \code{tests} will be parallelized using \pkg{BiocParallel}. Please note that parallelization is complex and depends on your operating system (Windows users might not see a gain or might even experience a slowdown).
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object to be used for parallelization. If \code{NULL} and \code{parallelized} = \code{TRUE}, the function will use a \code{\link[BiocParallel]{SerialParam}} object configured to use a single worker (core) and is therefore not parallelized, in order to prevent accidental use of large computation resources. Ignored if \code{parallelized} = \code{FALSE}.
#' @param tidy Logical. If \code{TRUE}, a \code{list} is returned, with \code{data.frame} objects corresponding to each test divided into each global, pairwise and/or lineage comparison results. An adjusted p-value (False Discovery Rate) is also calculated and each \code{data.frame} object is ordered based on decreasing Waldstat.
#' @param raw Logical. If \code{TRUE}, a \code{list} is returned, with \code{data.frame} objects corresponding to each test results.
#' @param verbose Logical. If \code{FALSE}, does not print progress messages and output, but warnings and errors will still be printed.
#'
#' @return A \pkg{list}, with \code{data.frame} objects corresponding to each test results and/or \code{data.frame} objects corresponding to each test and each global, pairwise and/or lineage comparison results.
#'
#' @import SingleCellExperiment
#' @import slingshot
#' @import tradeSeq
#' @import BiocParallel
#' @importFrom stats p.adjust
#' @importFrom SummarizedExperiment assays
#' @export

tradeSeqTests = function(models,
                         tests = c("assoc", "pattern", "diffend", "early", "startvsend", "condition"),
                         global = TRUE,
                         pairwise = TRUE,
                         lineages = TRUE,
                         l2fc = 0,
                         eigen.thresh = 0.01,
                         n.points = 2 * nknots(models),
                         knots = NULL,
                         pseudotime.values = NULL,
                         parallelized = FALSE,
                         BPPARAM = NULL,
                         tidy = TRUE,
                         raw = FALSE,
                         verbose = TRUE) {

  if (isTRUE("conditions" %in% colnames(colData(models)$tradeSeq))) {
    if (isTRUE(pairwise) & length(levels(factor(colData(models)$tradeSeq$conditions))) < 3) {
      if (isTRUE(verbose)) {
        message("Less than three conditions; skipping pairwise comparison in conditionTest")
      }
      pairwise2 = FALSE
    }
    else {
      pairwise2 = pairwise
    }
  }

  if(length(grep("pseudotime.Lineage", colnames(colData(models)$crv))) < 3) {
    if(isTRUE(pairwise)) {
      if (isTRUE(verbose)) {
        message("Less than three lineages; skipping pairwise comparison")
      }
      pairwise = pairwise2 = FALSE
    }
    global2 = TRUE
  }
  else {
    global2 = global
  }

  test.fun = function(test) {
    if (isTRUE(grepl("assoc", test))) {
      if (isTRUE(verbose)) {
        cat("Performing associationTest...", ifelse(isFALSE(parallelized),"\n",""))
      }
      return(list("associationTest" = as.data.frame(associationTest(models, global = global, lineages = lineages, l2fc = l2fc, nPoints = n.points))))
    }
    if (isTRUE(grepl("pattern", test))) {
      if (length(grep("pseudotime.Lineage", colnames(colData(models)$crv))) == 1) {
        if (isTRUE(verbose)) {
          message("Only one lineage; skipping patternTest")
        }
      }
      else {
        if (isTRUE(verbose)) {
          cat("Performing patternTest...", ifelse(isFALSE(parallelized),"\n",""))
        }
        return(list("patternTest" = as.data.frame(patternTest(models, global = global2, pairwise = pairwise, l2fc = l2fc, nPoints = n.points, eigenThresh = eigen.thresh))))
      }
    }
    if (isTRUE(grepl("early", test))) {
      if (length(grep("pseudotime.Lineage", colnames(colData(models)$crv))) == 1) {
        if (isTRUE(verbose)) {
          message("Only one lineage; skipping EarlyDETest")
        }
      }
      else {
        if(is.list(knots)) {
          return(unlist(lapply(seq_along(knots), function(n) {
            if (isTRUE(verbose)) {
              cat("Performing earlyDETest between knot ", knots[[n]][1], " and knot ", knots[[n]][2], "...", ifelse(n < length(knots),"\n",ifelse(isFALSE(parallelized),"\n","")), sep = "")
            }
            return(setNames(list(as.data.frame(earlyDETest(models, global = global2, pairwise = pairwise, l2fc = l2fc, nPoints = n.points, eigenThresh = eigen.thresh, knots = knots[[n]]))), paste0("earlyDETest", knots[[n]][1], "vs", knots[[n]][2])))
          }), recursive = FALSE))
        }
        else {
          if (isTRUE(verbose)) {
            cat("Performing earlyDETest...", ifelse(isFALSE(parallelized),"\n",""))
          }
          return(list("earlyDETest" = as.data.frame(earlyDETest(models, global = global2, pairwise = pairwise, l2fc = l2fc, nPoints = n.points, eigenThresh = eigen.thresh, knots = knots))))
        }
      }
    }
    if (isTRUE(grepl("diff", test))) {
      if (length(grep("pseudotime.Lineage", colnames(colData(models)$crv))) == 1) {
        if (isTRUE(verbose)) {
          message("Only one lineage; skipping diffEndTest")
        }
      }
      else {
        if (isTRUE(verbose)) {
          cat("Performing diffEndTest...", ifelse(isFALSE(parallelized),"\n",""))
        }
        return(list("diffEndTest" = as.data.frame(diffEndTest(models, global = global2, pairwise = pairwise, l2fc = l2fc))))
      }
    }
    if (isTRUE(grepl("start", test))) {
      if(is.list(pseudotime.values)) {
        return(unlist(lapply(seq_along(pseudotime.values), function(n) {
          if (isTRUE(verbose)) {
            cat("Performing startVsEndTest between the pseudotime value ", pseudotime.values[[n]][1], " and the pseudotime value ", pseudotime.values[[n]][2], "...", ifelse(n < length(pseudotime.values),"\n",ifelse(isFALSE(parallelized),"\n","")), sep = "")
          }
          return(setNames(list(as.data.frame(startVsEndTest(models, global = global, lineages = lineages, l2fc = l2fc, pseudotimeValues = pseudotime.values[[n]]))), paste0("startVsEndTest", pseudotime.values[[n]][1], "vs", pseudotime.values[[n]][2])))
        }), recursive = FALSE))
      }
      else {
        if (isTRUE(verbose)) {
          cat("Performing startVsEndTest...", ifelse(isFALSE(parallelized),"\n",""))
        }
        return(list("startVsEndTest" = as.data.frame(startVsEndTest(models, global = global, lineages = lineages, l2fc = l2fc, pseudotimeValues = pseudotime.values))))
      }
    }
    if(isTRUE(grepl("condition", test)) & isTRUE("conditions" %in% colnames(colData(models)$tradeSeq))) {
      if (length(grep("pseudotime.Lineage", colnames(colData(models)$crv))) == 1) {
        if (isTRUE(lineages)) {
          if (isTRUE(verbose)) {
            message("Only one lineage; skipping separate lineage comparison in conditionTest")
          }
          lineages = FALSE
        }
        global = TRUE
      }
      if(is.list(knots)) {
        return(unlist(lapply(seq_along(knots), function(n) {
          if (isTRUE(verbose)) {
            cat("Performing conditionTest between knot ", knots[[n]][1], " and knot ", knots[[n]][2], "...", ifelse(n < length(knots),"\n",ifelse(isFALSE(parallelized),"\n","")), sep = "")
          }
          return(setNames(list(as.data.frame(conditionTest(models, global = global, pairwise = pairwise2, lineages = lineages, l2fc = l2fc, eigenThresh = eigen.thresh, knots = knots[[n]]))), paste0("conditionTest", knots[[n]][1], "vs", knots[[n]][2])))
        }), recursive = FALSE))
      }
      else {
        if (isTRUE(verbose)) {
          cat("Performing conditionTest...", ifelse(isFALSE(parallelized),"\n",""))
        }
        return(list("conditionTest" = as.data.frame(conditionTest(models, global = global, pairwise = pairwise2, lineages = lineages, l2fc = l2fc, eigenThresh = eigen.thresh, knots = knots))))
      }
    }
  }

  if (isTRUE(parallelized)) {
    if (is.null(BPPARAM)) {
      warning("No BPPARAM parameter provided, using BiocParallel::SerialParam(), which is not parallelized", immediate. = TRUE)
      BPPARAM = SerialParam()
      parallelized = FALSE
    }
    tests.list = unlist(suppressWarnings(bplapply(tests, test.fun, BPPARAM = BPPARAM)), recursive = FALSE)
  }
  else {
    tests.list = unlist(lapply(tests, test.fun), recursive = FALSE)

  }

  if(isFALSE(tidy)) {
    if (isFALSE(raw)) {
      warning("Both tidy and raw parameters are FALSE; returning raw data")
    }
    return(tests.list)
  }

  if (isTRUE(verbose)) {
    cat("Tidying data...\n")
  }
  tests.list2 = list()
  if (length(grep("pseudotime.Lineage", colnames(colData(models)$crv))) > 1) {
    combinations = combn(1:length(grep("pseudotime.Lineage", colnames(colData(models)$crv))), 2, function(x) paste(x, collapse = "vs"))
  }
  else {
    combinations = NA
  }
  for(i in seq_along(tests.list)) {
    tests.list2[[i]] = tests.list[[i]]
    if(length(grep("pvalue", colnames(tests.list[[i]]))) > 1) {
      tmp.list = list()
      for(j in seq_along(grep("pvalue", colnames(tests.list[[i]])))) {
        tmp.list[[j]] = tests.list[[i]][ , ((j*3)-2):(j*3)]
        tmp.list[[j]]$padj = p.adjust(tmp.list[[j]][ , grep("pvalue",colnames(tmp.list[[j]])), drop = T], "fdr")
        if(isTRUE(grepl("startVsEnd", names(tests.list)[i]))) {
          if(isTRUE(global) & j > 1) {
            tmp.list[[j]] = cbind(tmp.list[[j]], tests.list[[i]][ , ncol(tests.list[[i]])-length(grep("pseudotime.Lineage", colnames(colData(models)$crv)))+j-1, drop = F])
            colnames(tmp.list[[j]]) = c("waldStat", "df", "pvalue", "padj", "LogFC")
          }
          else if(isFALSE(global)) {
            tmp.list[[j]] = cbind(tmp.list[[j]], tests.list[[i]][ , ncol(tests.list[[i]])-length(grep("pseudotime.Lineage", colnames(colData(models)$crv)))+j, drop = F])
            colnames(tmp.list[[j]]) = c("waldStat", "df", "pvalue", "padj", "LogFC")
          }
          else {
            colnames(tmp.list[[j]]) = c("waldStat", "df", "pvalue", "padj")
          }
        }
        else if(isTRUE(grepl("diffEnd", names(tests.list)[i]))) {
          if(isTRUE(global2) & j > 1) {
            tmp.list[[j]] = cbind(tmp.list[[j]], tests.list[[i]][ , ncol(tests.list[[i]])-length(combinations)+j-1, drop = F])
            colnames(tmp.list[[j]]) = c("waldStat", "df", "pvalue", "padj", "LogFC")
          }
          else if(isFALSE(global2)) {
            tmp.list[[j]] = cbind(tmp.list[[j]], tests.list[[i]][ , ncol(tests.list[[i]])-length(combinations)+j, drop = F])
            colnames(tmp.list[[j]]) = c("waldStat", "df", "pvalue", "padj", "LogFC")
          }
          else {
            colnames(tmp.list[[j]]) = c("waldStat", "df", "pvalue", "padj")
          }
        }
        else if(isTRUE(grepl("association", names(tests.list)[i]))) {
          tmp.list[[j]] = cbind(tmp.list[[j]], tests.list[[i]][ , ncol(tests.list[[i]]), drop = F])
          colnames(tmp.list[[j]]) = c("waldStat", "df", "pvalue", "padj", "meanLogFC")
        }
        else if(isTRUE(grepl("pattern|earlyDE", names(tests.list)[i]))) {
          tmp.list[[j]] = cbind(tmp.list[[j]], tests.list[[i]][ , ncol(tests.list[[i]]), drop = F])
          colnames(tmp.list[[j]]) = c("waldStat", "df", "pvalue", "padj", "medianLogFC")
        }
        else {
          colnames(tmp.list[[j]]) = c("waldStat", "df", "pvalue", "padj")
        }
        tmp.list[[j]] = tmp.list[[j]][order(-tmp.list[[j]][, 1]), ]
      }
      if(isTRUE(grepl("association|condition", names(tests.list)[i]))) {
        if (isTRUE("conditions" %in% colnames(colData(models)$tradeSeq)) && (isTRUE(grepl("association", names(tests.list)[i])) | (isTRUE(grepl("condition", names(tests.list)[i])) & length(levels(factor(colData(models)$tradeSeq$conditions))) > 2))) {
          if(length(tmp.list) > length(grep("pseudotime.Lineage", colnames(colData(models)$crv)))*length(levels(factor(colData(models)$tradeSeq$conditions)))) {
            if (length(grep("pseudotime.Lineage", colnames(colData(models)$crv))) == 1) {
              names(tmp.list) = c("Global", levels(factor(colData(models)$tradeSeq$conditions)))
            }
            else {
              names(tmp.list) = c("Global", apply(expand.grid(levels(factor(colData(models)$tradeSeq$conditions)), paste0("Lineage", 1:(length(grep("pseudotime.Lineage", colnames(colData(models)$crv)))))), 1, function(x) paste(x[2], x[1], sep = "_")))
            }
          }
          else {
            if (length(grep("pseudotime.Lineage", colnames(colData(models)$crv))) == 1) {
              names(tmp.list) = levels(factor(colData(models)$tradeSeq$conditions))
            }
            else {
              names(tmp.list) = apply(expand.grid(levels(factor(colData(models)$tradeSeq$conditions)), paste0("Lineage", 1:(length(grep("pseudotime.Lineage", colnames(colData(models)$crv)))))), 1, function(x) paste(x[2], x[1], sep = "_"))
            }
          }
        }
        else {
          if(length(tmp.list) > length(grep("pseudotime.Lineage", colnames(colData(models)$crv)))) {
            names(tmp.list) = c("Global", paste0("Lineage", 1:(length(grep("pseudotime.Lineage", colnames(colData(models)$crv))))))
          }
          else {
            names(tmp.list) = paste0("Lineage", 1:length(grep("pseudotime.Lineage", colnames(colData(models)$crv))))
          }
        }
      }
      else if(isTRUE(grepl("startVsEnd", names(tests.list)[i]))) {
        if(length(tmp.list) > length(grep("pseudotime.Lineage", colnames(colData(models)$crv)))) {
          names(tmp.list) = c("Global", paste0("Lineage", 1:(length(grep("pseudotime.Lineage", colnames(colData(models)$crv))))))
        }
        else {
          names(tmp.list) = paste0("Lineage", 1:length(grep("pseudotime.Lineage", colnames(colData(models)$crv))))
        }
      }
      else if(isTRUE(grepl("pattern|earlyDE|diffEnd", names(tests.list)[i]))) {
        if(length(tmp.list) > length(combinations)) {
          names(tmp.list) = c("Global", paste0("Lineage", combinations))
        }
        else {
          names(tmp.list) = paste0("Lineage", combinations)
        }
      }
      tests.list[[i]] = tmp.list
    }
    else {
      tests.list[[i]]$padj = p.adjust(tests.list[[i]]$pvalue, "fdr")
      tests.list[[i]] = tests.list[[i]][order(-tests.list[[i]]$waldStat), ]
      if(ncol(tests.list[[i]]) == 5) {
        tests.list[[i]] = tests.list[[i]][ , c(1:3, 5, 4)]
        if(isTRUE(grepl("startVsEnd|diffEnd", names(tests.list)[i]))) {
          colnames(tests.list[[i]]) = c("waldStat", "df", "pvalue", "padj", "LogFC")
        }
        else if(isTRUE(grepl("association", names(tests.list)[i]))) {
          colnames(tests.list[[i]]) = c("waldStat", "df", "pvalue", "padj", "meanLogFC")
        }
        else if(isTRUE(grepl("pattern|earlyDE", names(tests.list)[i]))) {
          colnames(tests.list[[i]]) = c("waldStat", "df", "pvalue", "padj", "medianLogFC")
        }
      }
      else {
        colnames(tests.list[[i]]) = c("waldStat", "df", "pvalue", "padj")
      }
    }
  }
  names(tests.list2) = names(tests.list)
  if (isFALSE(raw)) {
    return(tests.list)
  }
  else {
    return(list("Raw data" = tests.list2, "Tidy data" = tests.list))
  }
}
