#' @title Scatterplot of log-transformed counts, and average fitted GAM smoother curves along pseudotime
#'
#' @description This function generates a scatterplot of the log-transformed counts, and displays on top the average fitted GAM smoother curves, computed using \code{\link[tradeSeq]{fitGAM}}, along pseudotime, and is a reworked version of \code{link[tradeSeq]{plotSmoothers}}, it allows for the plotting of multiple genes simultaneously, as well as selected lineages and conditions. Additionally, the scatterplot is rasterized to reduce the size of the output file and branching points may be displayed on each curve.
#'
#' @param models A \pkg{SingleCellExperiment} object containing the fitted GAM smoothers, computed using \code{\link[tradeSeq]{fitGAM}}, with or without \code{conditions} provided.
#' @param predictSmooth.df A \code{data.frame} object containing the average fitted GAM smoothers, computed using \code{\link[tradeSeq]{predictSmooth}} with \code{nPoints}, \code{lineages} and \code{tidy} = \code{TRUE}. If \code{NULL}, average fitted GAM smoothers will be internally computed.
#' @param genes Character. The names of one or several genes to plot the log-transformed counts and the average fitted GAM smoothers from.
#' @param lineages Numeric. The indices of the lineages (for example, c(1, 5), 2:4 etc) to plot the log-transformed counts, the average fitted GAM smoothers and the pseudotime values from.
#' @param conditions Character. The names of one or several \code{conditions} identities to select. If \code{NULL}, all identities are used, and each unique condition will be plotted, for each of the \code{lineages} provided (for example, if you have four lineages and three conditions, twelve curves will be plotted). Ignored if the \code{models} object was computed using \code{\link[tradeSeq]{fitGAM}} without \code{conditions}. Please note that if the \code{models} object was computed using \code{\link[tradeSeq]{fitGAM}} with \code{conditions}, it is not possible to plot global lineages (without \code{conditions}), you would need to compute a new \code{models} object using \code{\link[tradeSeq]{fitGAM}} without \code{conditions}. This is due to a limitation in how \code{\link[tradeSeq]{predictSmooth}} returns average fitted GAM smoothers (if the \code{models} object was computed with \code{conditions}, the function will always mean the fitted GAM smoothers for each lineage and each condition independently).
#' @param facets Character. The names of one or two metadata among 'genes', 'lineages' and/or 'conditions' to facet the plot. If one metadata is provided, \code{\link[ggplot2]{facet_wrap}} will be used, and if two are provided, \code{\link[ggplot2]{facet_grid}} will be used, with the first metadata as rows and the second as columns. If the \code{models} object was computed using \code{\link[tradeSeq]{fitGAM}} without \code{conditions}, you cannot facet by 'conditions'. Please note that if the \code{models} object was computed using \code{\link[tradeSeq]{fitGAM}} with \code{conditions}, it is not possible to plot global lineages (without \code{conditions}), you would need to compute a new \code{models} object using \code{\link[tradeSeq]{fitGAM}} without \code{conditions}. This is due to a limitation in how \code{\link[tradeSeq]{predictSmooth}} returns average fitted GAM smoothers (if the \code{models} object was computed with \code{conditions}, the function will always mean the fitted GAM smoothers for each lineage and each condition independently).
#' @param nPoints Numeric. (from \code{\link[tradeSeq]{predictSmooth}} documentation) The number of points used to create the grid along the smoother for each lineage.
#' @param branch.points Numeric, Character or List. Branching points may be shown on the fitted GAM smoother curves, to partition the lineages and help visualize differences. Either one or more pseudotime values to plot the branching points at, and/or one or several values containing 'knot' followed by a number (for example, 'knot2', 'knot4' etc), which correspond to the knots (\code{k}) input in \code{\link[tradeSeq]{fitGAM}} and which divide each lineage into segments; the function will extract the pseudotime values from each knot number and plot the branching points. Mixing any of the two options is possible (for example, c(4.3, 7, 'knot4')). You may also provide a \code{list} containing any of the two options, and named after one or several identities of a single metadata among 'genes', 'lineages' or 'conditions' (for example, list('CD4' = c(2, 'knot5')) which will only plot branching points on CD4 curves). For \code{lineages}, the \code{list} names need to be provided as 'Lineage' (with capital L) followed by the index (for example, list('Lineage1' = c(4.3, 'knot4'), 'Lineage2' = c('knot2', 7, 9)), which will only plot branching points on lineage 1 and 2 curves). Finally, you may also provide identities corresponding to several metadata at the same time by pasting them together with '_' (for example, list('Lineage1_CD4' = c(1, 'knot1', 12)), which will only plot branching points on lineage 1 and CD4 curves).
#' @param points.size Numeric. The size of the log-transformed count points.
#' @param points.alpha Numeric. The transparency of the log-transformed count points.
#' @param curves.width Numeric. The width of the fitted GAM smoother curves.
#' @param curves.alpha Numeric. The transparency of the fitted GAM smoother curves.
#' @param facets.scales Character. (from \code{\link[ggplot2]{facet_grid}} documentation) Are scales shared across all facets (the default, "\code{fixed}"), or do they vary across rows ("\code{free_x}"), columns ("\code{free_y}"), or both rows and columns ("\code{free}")?
#' @param facets.axes Character. (from \code{\link[ggplot2]{facet_grid}} documentation) Determines which axes will be drawn. When "\code{margins}" (default), axes will be drawn at the exterior margins. "\code{all_x}" and "\code{all_y}" will draw the respective axes at the interior panels too, whereas "\code{all}" will draw all axes at all panels.
#' @param colors Character. The color names for each identity of the metadata excluded from \code{facets} among 'genes', 'lineages' and/or 'conditions' (for example, if \code{facets} = c('genes', 'conditions'), the \code{colors} will be for each of the \code{lineages} provided, or if \code{facets} = 'lineages', the \code{colors} will be for each of the \code{genes} times \code{conditions} provided). If \code{NULL}, uses \pkg{Seurat}'s default colors.
#' @param nrow Numeric. The number of rows to use for the facets. Ignored if two metadata are provided to \code{facets}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @import SingleCellExperiment
#' @import slingshot
#' @import tradeSeq
#' @import data.table
#' @import colorRamp2
#' @import ggplot2
#' @import ggborderline
#' @import scattermore
#' @import scales
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors metadata
#' @export

curveSmoothers = function(models,
                           predictSmooth.df = NULL,
                           genes = if (is.null(predictSmooth.df)) NA else unique(predictSmooth.df$gene),
                           lineages,
                           conditions = NULL,
                           facets = c("genes", "conditions"),
                           nPoints = 100,
                           branch.points = NULL,
                           points.size = 3,
                           points.alpha = 0.8,
                           curves.width = 1,
                           curves.alpha = 1,
                           facets.scales = "free",
                           facets.axes = "all",colors = NULL,
                           nrow = floor(sqrt(length(genes)))) {

  if (any(is.na(genes))) {
    stop("Please provide the genes names to plot the smoothed expression from")
  }

  yhat = time = NULL
  lineages = paste0("Lineage", lineages)
  knots = ceiling(unname(metadata(models)$tradeSeq$knots)*100)/100

  if (is.null(predictSmooth.df)) {
    predictSmooth.df = predictSmooth(models, gene = genes, nPoints = nPoints, tidy = TRUE)
  }
  predictSmooth.df$lineage = paste0("Lineage", predictSmooth.df$lineage)
  predictSmooth.df = predictSmooth.df[predictSmooth.df$lineage %in% lineages, , drop = FALSE]
  predictSmooth.df$yhat = log1p(predictSmooth.df$yhat)

  mat = assays(models)$counts[genes, , drop = FALSE]
  mat = as.data.frame(log1p(t(as.matrix(mat))))
  mat = cbind(mat, do.call(cbind, lapply(lineages, function(lin) {
    df = as.data.frame(colData(models)$crv[ , grep(lin, colnames(colData(models)$crv)), drop = FALSE])
    df[ , grep("pseudotime", colnames(df))] = ifelse(df[ , grep("cellWeights", colnames(df))] == 0, NA, df[ , grep("pseudotime", colnames(df))])
    return(setNames(data.frame(df[ , grep("pseudotime", colnames(df)), drop = TRUE]), lin))
  })))
  if (isTRUE("condition" %in% colnames(predictSmooth.df))) {
    mat$conditions = as.character(colData(models)$tradeSeq$conditions)
    if (is.character(conditions)) {
      mat = mat[mat$conditions %in% conditions, , drop = FALSE]
      colnames(predictSmooth.df) = c("lineages", "time", "conditions", "genes", "yhat")
      predictSmooth.df = predictSmooth.df[predictSmooth.df$conditions %in% conditions, , drop = FALSE]
    }
    else {
      colnames(predictSmooth.df) = c("lineages", "time", "conditions", "genes", "yhat")
      conditions = unique(predictSmooth.df$conditions)
    }
  }
  else {
    colnames(predictSmooth.df) = c("lineages", "time", "genes", "yhat")
    mat$conditions = predictSmooth.df$conditions = conditions = "therearenoconditionsinthisobject"
    if (length(setdiff(facets, "conditions")) == 0) {
      stop("The object does not contain conditions, facets need to be either 'genes', 'lineages' or both")
    }
    else {
      facets = setdiff(facets, "conditions")
    }
  }
  mat = na.omit(as.data.frame(melt(melt(setDT(mat),
                                        measure.vars = colnames(mat)[!grepl("Lineage|conditions", colnames(mat))],
                                        variable.name = "genes", value.name = "yhat"),
                                   measure.vars = colnames(mat)[grepl("Lineage", colnames(mat))],
                                   variable.name = "lineages", value.name = "time")))
  if (length(facets) == 2) {
    loop.var = ifelse(any(conditions == "therearenoconditionsinthisobject"),
                      ifelse(length(setdiff(c("genes", "lineages"), facets)) == 0,
                             "lineages", setdiff(c("genes", "lineages"), facets)),
                      setdiff(c("genes", "lineages", "conditions"), facets))
    predictSmooth.df[ , "loop.var"] = predictSmooth.df[ , loop.var]
    mat[ , "loop.var"] = mat[ , loop.var]
    loop.var = get(loop.var)
  }
  else {
    loop.var = gsub("_therearenoconditionsinthisobject", "", paste(rep(get(setdiff(c("genes", "lineages", "conditions"), facets)[1]),
                                                                       each = length(get(setdiff(c("genes", "lineages", "conditions"), facets)[2]))),
                                                                   get(setdiff(c("genes", "lineages", "conditions"), facets)[2]), sep = "_"))
    predictSmooth.df[ , "loop.var"] = gsub("_therearenoconditionsinthisobject", "", paste(predictSmooth.df[ , setdiff(c("genes", "lineages", "conditions"), facets)[1]],
                                                                                          predictSmooth.df[ , setdiff(c("genes", "lineages", "conditions"), facets)[2]], sep = "_"))
    mat[ , "loop.var"] = gsub("_therearenoconditionsinthisobject", "", paste(mat[ , setdiff(c("genes", "lineages", "conditions"), facets)[1]],
                                                                             mat[ , setdiff(c("genes", "lineages", "conditions"), facets)[2]], sep = "_"))
  }
  mat[c(facets, "loop.var")] = lapply(c(facets, "loop.var"),  function(col) factor(mat[[col]], levels = get(col)))
  predictSmooth.df[c(facets, "loop.var")] = lapply(c(facets, "loop.var"),  function(col) factor(predictSmooth.df[[col]], levels = get(col)))
  if (!is.null(branch.points)) {
    branch.points = do.call(rbind, lapply(seq_along(branch.points), function(bp) {
      lst = list(lineages, conditions, genes)
      if (!is.null(names(branch.points)[bp])) {
        lst = lapply(lst, function(vec) {
          vec = if (any(unlist(strsplit(names(branch.points)))[bp] %in% vec))
            vec[vec %in% unlist(strsplit(names(branch.points)))[bp]] else vec
        })
      }
      bp = apply(expand.grid(branch.points[[bp]], lst[[1]], lst[[2]], lst[[3]]), 1, function(vec) paste(vec, collapse = "_"))
      return(do.call(rbind, lapply(bp, function(vec) {
        vec = strsplit(vec, "_")[[1]]
        if (isTRUE(grepl("knot", vec[1]))) {
          vec[1] = knots[as.numeric(gsub("knot", "", vec[1]))]
        }
        df = predictSmooth.df[predictSmooth.df$lineages %in% vec[2] &
                                predictSmooth.df$conditions %in% vec[3] &
                                predictSmooth.df$genes %in% vec[4], , drop = FALSE]
        vec = as.numeric(ifelse(as.numeric(vec[1]) > max(df$time), nrow(df),
                                ifelse(as.numeric(vec[1]) < min(df$time), 1,
                                       Position(function(pos) pos >= as.numeric(vec[1]), df$time))))
        return(df[vec, , drop = FALSE])
      })))
    }))
  }

  if (is.null(colors)) {
    colors = hue_pal()(n = length(get(ifelse(any(conditions == "therearenoconditionsinthisobject"),
                                             ifelse(length(setdiff(c("lineages", "genes"), facets)) == 0,
                                                    "lineages", setdiff(c("lineages", "genes"), facets)),
                                             setdiff(c("lineages", "conditions", "genes"), facets))[1])))
  }
  if (is.null(names(colors))) {
    names(colors) = get(ifelse(any(conditions == "therearenoconditionsinthisobject"),
                               ifelse(length(setdiff(c("lineages", "genes"), facets)) == 0,
                                      "lineages", setdiff(c("lineages", "genes"), facets)),
                               setdiff(c("lineages", "conditions", "genes"), facets))[1])
  }
  if (length(colors) < length(loop.var)) {
    col.var = get(setdiff(c("genes", "conditions", "lineages"), facets)[1])
    if (length(colors)*length(col.var) == length(loop.var)) {
      colors = unlist(lapply(seq_along(colors), function(col) {
        col = setNames(rev(colorRamp2(c(0, 1), c("white", colors[col]))
                           (seq(0, 1, length.out = length(col.var) + 1))[-1]), paste(names(colors[col]), col.var, sep = "_"))
      }))
    }
    else {
      stop("Not enough colors supplied, please provide either", length(loop.var), "colors or", length(loop.var)/length(col.var), "colors")
    }
  }

  p = ggplot(NULL, aes(x = time, y = yhat,
                       col = factor(.data[["loop.var"]]))) +
    geom_scattermore(data = mat, pointsize = points.size + 0.1, alpha = points.alpha, show.legend = FALSE) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0)) +
    labs(x = "Pseudotime", y = "Log expression + 1", col = "") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)), breaks = pretty_breaks()) +
    scale_y_continuous(expand = c(0,0), breaks = pretty_breaks()) +
    scale_color_manual(values = colors)

  for (i in loop.var) {
    p = p +
      geom_borderline(data = predictSmooth.df[predictSmooth.df$loop.var %in% i, ],
                      linewidth = curves.width + 0.1, borderwidth = 0.3, alpha = curves.alpha)

    if (is.data.frame(branch.points)) {
      p = p +
        geom_point(data = branch.points[branch.points$loop.var %in% i, ],
                   size = curves.width*2 + 0.1, show.legend = FALSE)
    }
  }

  if (length(facets) == 2) {
    p = p +
      facet_grid(factor(.data[[facets[1]]])~
                   factor(.data[[facets[2]]]),
                 scales = facets.scales,
                 axes = facets.axes) +
      theme(panel.spacing = unit(12, "points"),
            aspect.ratio = 1)
  }
  else {
    p = p +
      facet_wrap(~factor(.data[[facets[1]]]),
                 scales = facets.scales,
                 axes = facets.axes,
                 nrow = nrow) +
      theme(panel.spacing = unit(12, "points"),
            aspect.ratio = 1)
  }
  return(p)
}
