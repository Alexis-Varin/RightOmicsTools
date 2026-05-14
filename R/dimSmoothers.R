#' @title Scatterplot of the dimensionality reduction and fitted GAM smoothers for each cell
#'
#' @description This function generates a scatterplot of the log-transformed counts, and displays on top the average fitted GAM smoother curves, computed using \code{\link[tradeSeq]{fitGAM}}, along pseudotime, and is a reworked version of \code{\link[tradeSeq]{plotSmoothers}}, it allows for the plotting of multiple genes simultaneously, as well as selected lineages and conditions. Additionally, the scatterplot is rasterized to reduce the size of the output file and branching points may be displayed on each curve. Finally, the log-transformed counts and average fitted GAM smoother curves are rescaled to the same range between genes to allow for a better comparison.
#'
#' @param sds A \pkg{SingleCellExperiment} object containing the pseudotime values of one or several lineages, computed using \code{\link[slingshot]{slingshot}} (usually, the input object to \code{\link[tradeSeq]{fitGAM}}).
#' @param models A \pkg{SingleCellExperiment} object containing the fitted GAM smoothers, computed using \code{\link[tradeSeq]{fitGAM}}, with or without \code{conditions} provided.
#' @param sds.conditions A list of \pkg{SingleCellExperiment} objects obtained using \code{\link[condiments]{slingshot_conditions}}. The conditions used to compute the pseudotime values must be the same as those used to fit the GAM smoothers in the \code{models} object, if provided.
#' @param genes Character. The names of one or several genes to plot the fitted GAM smoothers for each cell from. If the \code{models} object is not provided, the function will plot the log-transformed counts.
#' @param lineages Numeric. The indices of the lineages (for example, c(1, 5), 2:4 etc) to plot the pseudotime values from. Set to \code{NA} to remove all lineages.
#' @param clusters Character or Factor. Either the name of a metadata present in the \code{sds} object (for example, 'annotations', 'seurat_clusters', etc) to color the cells, or the identities, as character or factor, of length equal to the number of cells in the \code{sds} object. Set to 'pseudotime' to color the cells by their pseudotime values.
#' @param conditions Character. The names of one or several \code{conditions} identities to select. If \code{NULL}, all identities are used, and each unique condition will be plotted, for each of the \code{lineages} provided (for example, if you have four lineages and three conditions, twelve curves will be plotted). Ignored if the \code{models} object was computed using \code{\link[tradeSeq]{fitGAM}} without \code{conditions}. Please note that if the \code{models} object was computed using \code{\link[tradeSeq]{fitGAM}} with \code{conditions}, it is not possible to plot global lineages (without \code{conditions}), you would need to compute a new \code{models} object using \code{\link[tradeSeq]{fitGAM}} without \code{conditions}. This is due to a limitation in how \code{\link[tradeSeq]{predictSmooth}} returns average fitted GAM smoothers (if the \code{models} object was computed with \code{conditions}, the function will always mean the fitted GAM smoothers for each lineage and each condition independently).
#' @param knots Numeric. The indices of the knots to display on the \code{lineages}. Set to \code{NA} to remove all knots. Ignored if the \code{models} object is not provided.
#' @param show.knots.labels Logical. If \code{TRUE}, the knot numbers will be displayed using \code{\link[ggrepel]{geom_label_repel}}. Ignored if the \code{models} object is not provided.
#' @param rescale Logical. If \code{TRUE}, fitted GAM smoothers will be adjusted using \code{\link[scales]{rescale}} between the first numerical value of \code{rescale.range} (lowest value) and the second numerical value (highest value). This is different than \code{\link[base]{scale}} as this doesn't compare values to any mean or standard deviation and is therefore not a Z-score, it only refits each value (independently for each gene) in order to visualize all \code{genes} in the same dimension regardless of their differences in fitted GAM smoothers.
#' @param rescale.range Numeric. The minimum and maximum values to resize the fitted GAM smoothers and internally passed to \code{\link[scales]{rescale}}. These values are arbitrary and will not change the visualization, only the values in the legend. Ignored if \code{rescale} = \code{FALSE}.
#' @param pseudotime.type Character. Determines the \code{cells.colors} range scale when \code{clusters} = 'pseudotime' and not all \code{lineages} are displayed: either 'absolute', where the \code{cells.colors} range is based on the highest pseudotime value among the displayed \code{lineages}, or 'relative', where the \code{cells.colors} range is based on the highest pseudotime value of all the \code{lineages}, including the non-displayed ones.
#' @param points.size Numeric. The size of the cells.
#' @param lineages.width Numeric. The width of the lineage curves.
#' @param lineages.arrow An \code{\link[grid]{arrow}} object added at the end of the lineage curves. Use the function's parameters to modify its appearance. Set to \code{NULL} to remove the arrow.
#' @param knots.size Numeric. The size of the knots on the lineage curves. Ignored if the \code{models} object is not provided.
#' @param cells.colors Character. If \code{genes} is provided, either two color names to use for the cells, corresponding to the lowest and highest values in the fitted GAM smoothers and internally passed to \code{\link[ggplot2]{scale_color_gradient}}, or the name of a palette and internally passed to \code{\link[ggplot2]{scale_color_viridis_c}} (such as 'turbo', 'inferno', 'magma', 'viridis' etc, check the function for all palettes available). If \code{clusters} is provided, the color names for each identity in \code{clusters}. If \code{NULL}, uses \pkg{Seurat}'s default colors. If \code{clusters} = 'pseudotime', the name of a palette and internally passed to \code{\link[ggplot2]{scale_color_viridis_c}}. If neither \code{genes} nor \code{clusters} are provided, the name of a color to use for all cells. If \code{NULL}, defaults to 'lightgrey'.
#' @param na.color Character. The name of a color to use for cells not participating in the \code{lineages} (for example, if you have four lineages and you only plot lineages 1 and 3, the cells participating in lineages 2 and 4 will be colored with \code{na.color}). Set to 'transparent' to not display these cells.
#' @param show.na Logical. If \code{FALSE}, ignores \code{na.color} and colors all cells with \code{cells.colors}, even those not participating in the \code{lineages}.
#' @param lineages.colors Character. either a single color name to use for all lineages, or the color names for each of the \code{lineages} displayed, or for the number of \code{lineages} times the number of \code{conditions} if the \code{models} object was computed with \code{conditions}.
#' @param knots.colors Character. either a single color name to use for all knots, or a vector of color names of length equal to the number of \code{lineages}. Ignored if the \code{models} object is not provided.
#' @param knots.labels.text.colors Character. either a single color name to use for all knot labels, or a vector of color names of length equal to the number of \code{lineages}. Ignored if the \code{models} object is not provided or if \code{show.knots.labels} = \code{FALSE}.
#' @param knots.labels.fill.colors Character. either a single color name to use for all knot label backgrounds, or a vector of color names of length equal to the number of \code{lineages}. Ignored if the \code{models} object is not provided or if \code{show.knots.labels} = \code{FALSE}.
#' @param axis.text.size Numeric. The font size of the dimensionality reduction values.
#' @param axis.title.size Numeric. The font size of the axis title.
#' @param plot.title.size Numeric. The font size of the plot title.
#' @param lineages.legend.names Character. You may provide custom names for the lineages displayed in the legend. Ignored if \code{lineages.colors} is a single color.
#' @param legend.text.size Numeric. The font size of the legend text.
#' @param raster Logical. If \code{TRUE}, the cells will be rasterized using \code{\link[ggrastr]{rasterize}} to reduce the size of the output file. The points may appear slightly blurry.
#'
#' @return A \pkg{ggplot2} object, or a \code{patchwork} object containing \pkg{ggplot2} objects.
#'
#' @import SingleCellExperiment
#' @import slingshot
#' @import tradeSeq
#' @import data.table
#' @import ggplot2
#' @import ggborderline
#' @import ggnewscale
#' @import ggrastr
#' @import ggrepel
#' @import patchwork
#' @import scales
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors metadata
#' @export

dimSmoothers = function(sds,
                        models = NULL,
                        sds.conditions = NULL,
                        genes = NULL,
                        lineages = NULL,
                        clusters = NULL,
                        conditions = NULL,
                        knots = NULL,
                        show.knots.labels = FALSE,
                        rescale = FALSE,
                        rescale.range = c(0, 3),
                        pseudotime.type = "relative",
                        points.size = 1,
                        lineages.width = 1,
                        lineages.arrow = arrow(),
                        knots.size = lineages.width*3,
                        cells.colors = if (is.character(genes)) "turbo" else if (any(clusters == "pseudotime")) "viridis" else NULL,
                        na.color = "lightgrey",
                        show.na = TRUE,
                        lineages.colors = "grey50",
                        knots.colors = lineages.colors,
                        knots.labels.text.colors = "white",
                        knots.labels.fill.colors = lineages.colors,
                        axis.text.size = 9,
                        axis.title.size = 11,
                        plot.title.size = 16,
                        lineages.legend.names = NULL,
                        legend.text.size = 9,
                        raster = FALSE) {

  lincols = NULL

  if (isTRUE("SingleCellExperiment" %in% class(models))) {
    if (is.numeric(knots)) {
      knots2 = knots
      knots = unname(metadata(models)$tradeSeq$knots)[knots]
    }
    else if (!any(is.na(knots))) {
      knots = unname(metadata(models)$tradeSeq$knots)
      knots2 = seq_len(length(knots))
    }
  }

  if (isTRUE("SingleCellExperiment" %in% class(models)) && isTRUE("conditions" %in% colnames(colData(models)$tradeSeq))) {
    if (is.null(conditions)) {
      conditions = levels(as.factor(colData(models)$tradeSeq$conditions))
    }
    conditions = setNames(factor(colData(models)$tradeSeq$conditions[colData(models)$tradeSeq$conditions %in% conditions], levels = conditions), colnames(models))
    conditions = split(seq_along(conditions), conditions)
    if (!is.list(sds.conditions)) {
      message("The models object was computed with conditions, splitting the sds object accordingly\nFor more accurate trajectories within conditions, consider using slingshot_conditions from condiments")
    }
    else {
      for (i in seq_along(sds.conditions)) {
        if (length(unique(conditions[colnames(sds.conditions[[i]])])) > 1) {
          stop("The sds.conditions and the models objects were not computed with the same conditions")
        }
      }
    }
  }
  else conditions = split(seq_along(colnames(sds)), factor(rep("all", ncol(sds))))

  pseudo = data.frame(slingPseudotime(sds))
  pseudo = ceiling(pseudo*100)/100
  if (is.null(lineages)) {
    lineages = seq_len(ncol(pseudo))
  }
  if (!any(is.na(lineages))) {
    lineages = as.numeric(lineages)
    pseudo = pseudo[ , lineages, drop = F]
    if (is.null(names(lineages.colors)) & length(lineages.colors) > 1) {
      if (length(conditions) > 1) {
        names(lineages.colors) = paste(rep(colnames(pseudo), each = length(conditions)), names(conditions), sep = "_")
      }
      else {
        names(lineages.colors) = colnames(pseudo)
      }
      if (!is.character(lineages.legend.names)) {
        lineages.legend.names = names(lineages.colors)
      }
    }
    else if (is.null(names(lineages.colors))) {
      names(lineages.colors) = lineages.colors
    }
  }
  missing.cells = which(apply(pseudo, 1, function(row) all(is.na(row))))

  data = data.frame(slingReducedDim(sds))
  if (is.character(genes) & isTRUE("SingleCellExperiment" %in% class(models))) {
    mat = log1p(t(predictCells(models, gene = genes)))
    if (isTRUE(rescale)) {
      mat = apply(mat, 2, rescale, to = rescale.range)
    }
    data = cbind(data, data.frame(mat))
  }
  else if (is.character(genes)) {
    mat = log1p(t(as.matrix(assays(sds)$counts[genes, , drop = F])))
    if (isTRUE(rescale)) {
      mat = apply(mat, 2, rescale, to = rescale.range)
    }
    data = cbind(data, data.frame(mat))
  }

  if (!any(is.na(lineages))) {
    if (length(conditions) > 1) {
      curves = list()
      knots.list = list()
      for (i in 1:length(conditions)) {
        if (is.list(sds.conditions)) {
          curves[[i]] = slingCurves(sds.conditions[[i]])[lineages]
        }
        else {
          curves[[i]] = slingCurves(sds)[lineages]
        }
        curves[[i]] = lapply(seq_along(curves[[i]]), function(cu) {
          df = as.data.frame(curves[[i]][[cu]]$s[curves[[i]][[cu]]$ord, , drop = F])
          df$condition = names(conditions)[i]
          df$lincols = ifelse(length(lineages.colors) == 1, names(lineages.colors), names(lineages.colors)[cu+length(cu)*(i-1)])
          return(df)
        })
        if (!any(is.na(knots)) & isTRUE("SingleCellExperiment" %in% class(models))) {
          knots.list[[i]] = lapply(seq_along(curves[[i]]), function(pts) {
            cu = curves[[i]][[pts]]
            time = pseudo[conditions[[i]], grep(paste0("Lineage", lineages[pts]), colnames(pseudo)) , drop = F]
            cu$time = seq(min(time, na.rm = T), max(time, na.rm = T), length.out = nrow(cu))
            cu$display = "yes"
            if (!is.null(lineages.arrow)) {
              cu$display[nrow(cu)] = "no"
            }
            cu$kshape = 21
            cu$kshape[1] = 23
            cu$kshape[nrow(cu)] = 22
            coords = sapply(knots, function(pt) {
              as.numeric(ifelse(as.numeric(pt) > max(cu$time), nrow(cu),
                                ifelse(as.numeric(pt) < min(cu$time), 1,
                                       Position(function(pos) pos >= as.numeric(pt), cu$time))))
            })
            cu = cu[coords, , drop = F]
            cu$kcols = ifelse(length(knots.colors) == 1, knots.colors, knots.colors[pts+length(pts)*(i-1)])
            cu$labtxtcols = ifelse(length(knots.labels.text.colors) == 1, knots.labels.text.colors, knots.labels.text.colors[pts])
            cu$labflcols = ifelse(length(knots.labels.fill.colors) == 1, knots.labels.fill.colors, knots.labels.fill.colors[pts])
            cu$label = knots2
            if (length(coords) > 1) {
              cu = cu[c(TRUE, sapply(seq_along(coords)[2:length(coords)], function(pt) if (coords[pt] == coords[pt-1]) FALSE else TRUE)), , drop = F]
            }
            return(cu)
          })
        }
      }
      curves = lapply(seq_along(curves[[1]]), function(cu) {
        do.call(rbind, lapply(curves, function(cond) cu = cond[[cu]]))
      })
      knots = lapply(seq_along(knots.list[[1]]), function(cu) {
        do.call(rbind, lapply(knots.list, function(cond) cond[[cu]]))
      })
    }
    else {
      curves = slingCurves(sds)[lineages]
      curves = lapply(seq_along(curves), function(cu) {
        df = as.data.frame(curves[[cu]]$s[curves[[cu]]$ord, , drop = F])
        df$condition = "all"
        df$lincols = ifelse(length(lineages.colors) == 1, names(lineages.colors), names(lineages.colors)[cu])
        return(df)
      })
      if (!any(is.na(knots)) & isTRUE("SingleCellExperiment" %in% class(models))) {
        knots = lapply(seq_along(curves), function(pts) {
          cu = curves[[pts]]
          time = pseudo[ , grep(paste0("Lineage", lineages[pts]), colnames(pseudo)) , drop = F]
          cu$time = seq(min(time, na.rm = T), max(time, na.rm = T), length.out = nrow(cu))
          cu$display = "yes"
          if (!is.null(lineages.arrow)) {
            cu$display[nrow(cu)] = "no"
          }
          cu$kshape = 21
          cu$kshape[1] = 23
          cu$kshape[nrow(cu)] = 22
          coords = sapply(knots, function(pt) {
            as.numeric(ifelse(as.numeric(pt) > max(cu$time), nrow(cu),
                              ifelse(as.numeric(pt) < min(cu$time), 1,
                                     Position(function(pos) pos >= as.numeric(pt), cu$time))))
          })
          cu = cu[coords, , drop = F]
          cu$kcols = ifelse(length(knots.colors) == 1, knots.colors, knots.colors[pts])
          cu$labtxtcols = ifelse(length(knots.labels.text.colors) == 1, knots.labels.text.colors, knots.labels.text.colors[pts])
          cu$labflcols = ifelse(length(knots.labels.fill.colors) == 1, knots.labels.fill.colors, knots.labels.fill.colors[pts])
          cu$label = knots2
          if (length(coords) > 1) {
            cu = cu[c(TRUE, sapply(seq_along(coords)[2:length(coords)], function(pt) if (coords[pt] == coords[pt-1]) FALSE else TRUE)), , drop = F]
          }
          return(cu)
        })
      }
    }
  }

  if ((is.character(clusters) | is.factor(clusters)) & !is.character(genes)) {
    if (clusters == "pseudotime") {
      if (any(is.na(lineages))) stop("lineages should be provided if clusters is set to 'pseudotime'")
      pseudo[] = lapply(pseudo, function(col) {
        col[missing.cells] = 0
        return(col)
      })
      if (pseudotime.type == "relative") {
        pseudo = data.frame(slingPseudotime(sds))
        pseudo = ceiling(pseudo*100)/100
      }
      data$clusters = apply(pseudo, 1, max, na.rm = TRUE)
      gglimits = c(min(data$clusters, na.rm = T), max(data$clusters, na.rm = T))
      if (length(cells.colors) == 1) {
        ggcolors = scale_color_viridis_c(option = cells.colors, na.value = na.color,
                                         limits = gglimits, guide = guide_colourbar(title.hjust = 0.5))
      }
      else {
        ggcolors = scale_color_gradient(low = cells.colors[1], high = cells.colors[2],
                                        na.value = na.color, limits = gglimits, guide = guide_colourbar(title.hjust = 0.5))
      }
    }
    else if (length(clusters) == 1) {
      data$clusters = colData(sds)[ , clusters]
      if (is.null(cells.colors)) {
        cells.colors = hue_pal()(length(unique(data$clusters)))
      }
      if (is.null(names(cells.colors))) {
        if (is.factor(data$clusters)) {
          names(cells.colors) = levels(data$clusters)
        }
        else {
          names(cells.colors) = unique(data$clusters)
        }
      }
      ggcolors = scale_color_manual(values = cells.colors, breaks = names(cells.colors),
                                    na.value = na.color, guide = guide_legend(title.hjust = 0.5,
                                                                              override.aes = list(size = points.size*2)))
    }
    else if (length(clusters) == ncol(sds)) {
      data$clusters = clusters
      if (is.null(cells.colors)) {
        cells.colors = hue_pal()(length(unique(data$clusters)))
      }
      if (is.null(names(cells.colors))) {
        if (is.factor(data$clusters)) {
          names(cells.colors) = levels(data$clusters)
        }
        else {
          names(cells.colors) = unique(data$clusters)
        }
      }
      ggcolors = scale_color_manual(values = cells.colors, breaks = names(cells.colors),
                                    na.value = na.color, guide = guide_legend(title.hjust = 0.5,
                                                                              override.aes = list(size = points.size*2)))
    }
    else stop("clusters should be either the name of a metadata present in the sds object, or a vector of identities of length equal to the number of cells in the sds object")
  }
  else if (!is.character(genes)) {
    data$clusters = "All"
    ggcolors = scale_color_manual(values = ifelse(is.null(cells.colors), "lightgrey", cells.colors),
                                  na.value = na.color)
  }

  gg = lapply(seq_len(ncol(data))[3:ncol(data)], function(col) {
    dat = data[, c(1, 2, col)]
    if (is.character(genes)) {
      gglimits = c(min(dat[[3]], na.rm = T), max(dat[[3]], na.rm = T))
      if (length(cells.colors) == 1) {
        ggcolors = scale_color_viridis_c(option = cells.colors, na.value = na.color,
                                         limits = gglimits, guide = guide_colourbar(title.hjust = 0.5))
      }
      else {
        ggcolors = scale_color_gradient(low = cells.colors[1], high = cells.colors[2],
                                        na.value = na.color, limits = gglimits, guide = guide_colourbar(title.hjust = 0.5))
      }
    }
    if (isTRUE(show.na)) {
      dat[[3]][missing.cells] = NA
    }
    if (isTRUE("SingleCellExperiment" %in% class(models)) && isTRUE("conditions" %in% colnames(colData(models)$tradeSeq))) {
      dat$condition = colData(models)$tradeSeq$conditions
    }
    else dat$condition = "all"
    g = ggplot(dat, aes(x = .data[[colnames(dat)[1]]], y = .data[[colnames(dat)[2]]], col = .data[[colnames(dat)[3]]])) +
      labs(col = ifelse(any(clusters == "pseudotime"), "Pseudotime",
                        ifelse(is.character(clusters) | is.factor(clusters), "Clusters",
                               ifelse(isTRUE("SingleCellExperiment" %in% class(models)),
                                      "Smoothed\nExpression", "Logged\nExpression")))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = plot.title.size))

    if (is.character(genes)) {
      g = g + ggtitle(colnames(dat)[3])
    }

    if (isTRUE(raster)) {
      g = g +
        ggcolors +
        rasterize(geom_point(data = dat[is.na(dat[[3]]), ], size = points.size/2, show.legend = ifelse(is.null(genes) & is.null(clusters), FALSE, TRUE))) +
        rasterize(geom_point(data = dat[!is.na(dat[[3]]), ], size = points.size/2, show.legend = ifelse(is.null(genes) & is.null(clusters), FALSE, TRUE)))
    }
    else {
      g = g +
        ggcolors +
        geom_point(data = dat[is.na(dat[[3]]), ], size = points.size, show.legend = ifelse(is.null(genes) & is.null(clusters), FALSE, TRUE)) +
        geom_point(data = dat[!is.na(dat[[3]]), ], size = points.size, show.legend = ifelse(is.null(genes) & is.null(clusters), FALSE, TRUE))
    }
    if (!any(is.na(lineages))) {
      if (!any(is.na(knots)) & isTRUE("SingleCellExperiment" %in% class(models))) {
        for (cu in seq_along(curves)) {
          if (length(lineages.colors) == 1) {
            g = g + geom_borderpath(data = curves[[cu]], col = lineages.colors,
                                    linewidth = lineages.width, arrow = lineages.arrow) +
              geom_point(data = knots[[cu]][knots[[cu]]$display == "yes", ], stroke = lineages.width,
                         fill = knots[[cu]][knots[[cu]]$display == "yes", "kcols"], col = "white",
                         size = knots.size, shape = knots[[cu]][knots[[cu]]$display == "yes", "kshape"])
          }
          else {
            if (cu == 1) {
              g = g + new_scale_color()
            }
            g = g +
              geom_path(data = curves[[cu]], aes(x = .data[[colnames(curves[[cu]])[1]]],
                                                 y = .data[[colnames(curves[[cu]])[2]]],
                                                 col = lincols),
                        linewidth = lineages.width, inherit.aes = FALSE) +
              geom_borderpath(data = curves[[cu]], aes(x = .data[[colnames(curves[[cu]])[1]]],
                                                       y = .data[[colnames(curves[[cu]])[2]]],
                                                       col = lincols),
                              linewidth = lineages.width, inherit.aes = FALSE,
                              arrow = lineages.arrow, show.legend = FALSE) +
              geom_point(data = knots[[cu]][knots[[cu]]$display == "yes", ], aes(x = .data[[colnames(knots[[cu]])[1]]],
                                                                                 y = .data[[colnames(knots[[cu]])[2]]]),
                         fill = knots[[cu]][knots[[cu]]$display == "yes", "kcols"], col = "white", stroke = lineages.width,
                         size = knots.size, shape = knots[[cu]][knots[[cu]]$display == "yes", "kshape"], inherit.aes = FALSE)
            if (cu == length(curves)) {
              g = g +
                scale_color_manual(values = lineages.colors, labels = lineages.legend.names,
                                   name = "", guide = guide_legend(override.aes = list(linewidth = 2)))
            }
          }
        }
        if (isTRUE(show.knots.labels)) {
          knots = do.call(rbind, knots)
          g = g +
            geom_label_repel(data = knots, aes(x = .data[[colnames(knots)[1]]],
                                               y = .data[[colnames(knots)[2]]]),
                             label = knots$label, col = knots$labtxtcols,
                             fill = knots$labflcols, max.overlaps = Inf,
                             force = 3, segment.linetype = 0, inherit.aes = FALSE)
        }
      }
      else {
        for (cu in seq_along(curves)) {
          if (length(lineages.colors) == 1) {
            g = g + geom_borderpath(data = curves[[cu]], col = lineages.colors,
                                    linewidth = lineages.width, arrow = lineages.arrow)
          }
          else {
            if (cu == 1) {
              g = g + new_scale_color()
            }
            g = g + geom_path(data = curves[[cu]], aes(x = .data[[colnames(curves[[cu]])[1]]],
                                                       y = .data[[colnames(curves[[cu]])[2]]],
                                                       col = lincols),
                              linewidth = lineages.width, inherit.aes = FALSE) +
              geom_borderpath(data = curves[[cu]], aes(x = .data[[colnames(curves[[cu]])[1]]],
                                                       y = .data[[colnames(curves[[cu]])[2]]],
                                                       col = lincols),
                              linewidth = lineages.width, inherit.aes = FALSE,
                              arrow = lineages.arrow, show.legend = FALSE)
            if (cu == length(curves)) {
              g = g +
                scale_color_manual(values = lineages.colors, labels = lineages.legend.names,
                                   name = "", guide = guide_legend(override.aes = list(linewidth = 2)))
            }
          }
        }
      }
    }
    if (length(conditions) > 1) {
      g = g + facet_wrap(~factor(condition, levels = names(conditions)), axes = "all", nrow = 1)
    }
    return(g)
  })

  if (length(gg) > 1) {
    gg = wrap_plots(gg, ncol = 1, guides = "collect")
  }
  else gg = gg[[1]]

  return(gg)
}
