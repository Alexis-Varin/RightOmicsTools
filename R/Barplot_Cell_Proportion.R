#' @title Bar plot of the cell proportions in each identity
#'
#' @description This function generates a bar plot of the cell proportions from a \pkg{Seurat} object.
#'
#' @param seurat_object A \pkg{Seurat} object.
#' @param group.by Character. The name of a metadata (for example, 'orig.ident', 'seurat_clusters', etc) to group the identities of the active.ident metadata by. Returns a stacked bar plot.
#' @param split.by Character. The name of a metadata to split the identities of the active.ident metadata by. Returns separate \code{ggplot2} objects.
#' @param idents Character. The names of one or several identities in the active.ident metadata to select. If \code{NULL}, all identities are used.
#' @param group.idents Character. The names of one or several \code{group.by} identities to select. If \code{NULL}, all identities are used. Ignored if \code{group.by} = \code{NULL}.
#' @param split.idents Character. The names of one or several \code{split.by} identities to select. If \code{NULL}, all identities are used. Ignored if \code{split.by} = \code{NULL}.
#' @param colors Character. The color names for each identity of the active.ident metadata or in \code{idents}. If \code{NULL}, uses \pkg{Seurat}'s default colors.
#' @param order.prop Character or Numeric. Either 'reverse', or the identities (as names or as numeric values corresponding to the indices) of the active.ident metadata or in \code{idents} to order the cells.
#' @param order.group Character. Either 'reverse', or the names of the \code{group.by} identities or in \code{group.idents} to order the cells. Ignored if \code{group.by} = \code{NULL}.
#' @param order.split Character. Either 'reverse', or the names of the \code{split.by} identities or in \code{split.idents} to order the cells. Ignored if \code{split.by} = \code{NULL}.
#' @param order.colors Logical. If \code{TRUE}, the \code{colors} will automatically be ordered according to \code{order.prop}. Ignored if \code{order.prop} = \code{NULL}.
#' @param alpha Numeric. The transparency of the bar colors. A value between 0 and 1.
#' @param show.cellsum.label Logical. If \code{TRUE}, the cell sum will be shown at the top of each bar.
#' @param cellsum.label.size Numeric. The font size of the cell sum label. Ignored if \code{show.cellsum.label} = \code{FALSE}.
#' @param axis.text.size Numeric. The font size of the identity names and cell percent or count.
#' @param x.axis.angle Numeric. The rotation angle of the identity names.
#' @param x.axis.hjust Numeric. The horizontal justification of the identity names.
#' @param y.axis.title.size Numeric. The font size of the y axis title.
#' @param legend.text.size Numeric. The font size of the legend text. Ignored if \code{show.legend} = \code{FALSE}.
#' @param legend.side Character. The side where the legend will be displayed, either 'left', 'right', 'top' or 'bottom'. Ignored if \code{show.legend} = \code{FALSE}.
#' @param split.plot.title.size Numeric. The font size of the split plot titles. Ignored if \code{split.by} = \code{NULL}.
#' @param show.legend Logical. If \code{TRUE}, shows the legend.
#' @param prop.percent Logical. If \code{TRUE}, the cell proportions will be shown as a percentage of total cells for each identity. Ignored if \code{group.by} = \code{NULL}.
#' @param nrow Numeric. The number of rows in the \code{patchwork} object. Ignored if \code{group.by} = \code{NULL} and \code{split.by} = \code{NULL}.
#' @param unique.group.plot Logical. If \code{TRUE}, the stacked bar plots will be gathered in a single \code{ggplot} object. Ignored if \code{group.by} = \code{NULL}.
#' @param unique.split.plot Logical. If \code{TRUE}, the \code{ggplot} objects will be gathered in a single \code{patchwork}. Ignored if \code{split.by} = \code{NULL}.
#' @param output.data Logical. If \code{TRUE}, the function will return a \code{data.frame} object or a \code{list} of \code{data.frame} objects of the cell proportions, instead of displaying anything.
#'
#' @return A \code{ggplot} object, a \code{list} of \code{ggplot} objects, a \code{patchwork} object containing \code{ggplot} objects or a \code{list} of \code{patchwork} objects containing \code{ggplot} objects.
#'
#' @examples
#' \dontshow{
#' suppressWarnings(suppressPackageStartupMessages(library(Seurat)))
#' }
#' # Prepare data
#' pbmc3k <- Right_Data("pbmc3k")
#'
#' # Example 1: default parameters
#' Barplot_Cell_Proportion(pbmc3k)
#' @import Seurat
#' @import SeuratObject
#' @import ggplot2
#' @import patchwork
#' @import scales
#' @import grDevices
#' @export

Barplot_Cell_Proportion = function(seurat_object,
                                   group.by = NULL,
                                   split.by = NULL,
                                   idents = NULL,
                                   group.idents = NULL,
                                   split.idents = NULL,
                                   colors = NULL,
                                   order.prop = NULL,
                                   order.group = NULL,
                                   order.split = NULL,
                                   order.colors = TRUE,
                                   alpha = 1,
                                   show.cellsum.label = TRUE,
                                   cellsum.label.size = 3,
                                   axis.text.size = 9,
                                   x.axis.angle = 60,
                                   x.axis.hjust = 1,
                                   y.axis.title.size = 11,
                                   legend.text.size = 9,
                                   legend.side = "bottom",
                                   show.legend = TRUE,
                                   split.plot.title.size = 24,
                                   prop.percent = TRUE,
                                   nrow = 1,
                                   unique.group.plot = TRUE,
                                   unique.split.plot = FALSE,
                                   output.data = FALSE) {

  ident1 = ident2 = sumpercent = nbcells = NULL

  if (is.null(group.by)) {
    prop.percent = FALSE
  }

  if (is.character(group.idents) & !is.character(group.by)) {
    group.idents = NULL
  }
  if (is.character(split.idents) & !is.character(split.by)) {
    split.idents = NULL
  }

  if (is.character(idents)) {
    if (is.character(group.idents)) {
      if (is.character(split.idents)) {
        activeident = Idents(seurat_object)
        seurat_object@meta.data$id_group_split = paste(as.character(seurat_object@active.ident),
                                                       as.character(seurat_object[[group.by]][ , 1, drop = TRUE]),
                                                       as.character(seurat_object[[split.by]][ , 1, drop = TRUE]), sep = "_")
        idents = paste(rep(idents, each = length(group.idents)), group.idents, sep = "_")
        idents = paste(rep(idents, each = length(split.idents)), split.idents, sep = "_")
        Idents(seurat_object) = "id_group_split"
        seurat_object = subset(seurat_object, idents = idents)
        Idents(seurat_object) = activeident
      }
      else {
        activeident = Idents(seurat_object)
        seurat_object@meta.data$id_group = paste(as.character(seurat_object@active.ident),
                                                 as.character(seurat_object[[group.by]][ , 1, drop = TRUE]), sep = "_")
        idents = paste(rep(idents, each = length(group.idents)), group.idents, sep = "_")
        Idents(seurat_object) = "id_group"
        seurat_object = subset(seurat_object, idents = idents)
        Idents(seurat_object) = activeident
      }
    }
    else if (is.character(split.idents)) {
      activeident = Idents(seurat_object)
      seurat_object@meta.data$id_split = paste(as.character(seurat_object@active.ident),
                                               as.character(seurat_object[[split.by]][ , 1, drop = TRUE]), sep = "_")
      idents = paste(rep(idents, each = length(split.idents)), split.idents, sep = "_")
      Idents(seurat_object) = "id_split"
      seurat_object = subset(seurat_object, idents = idents)
      Idents(seurat_object) = activeident
    }
    else {
      seurat_object = subset(seurat_object, idents = idents)
    }
  }
  else if (is.character(group.idents)) {
    if (is.character(split.idents)) {
      activeident = Idents(seurat_object)
      seurat_object@meta.data$group_split = paste(as.character(seurat_object[[group.by]][ , 1, drop = TRUE]),
                                                  as.character(seurat_object[[split.by]][ , 1, drop = TRUE]), sep = "_")
      group.idents = paste(rep(group.idents, each = length(split.idents)), split.idents, sep = "_")
      Idents(seurat_object) = "group_split"
      seurat_object = subset(seurat_object, idents = group.idents)
      Idents(seurat_object) = activeident
    }
    else {
      activeident = Idents(seurat_object)
      Idents(seurat_object) = as.character(seurat_object[[group.by]][ , 1, drop = TRUE])
      seurat_object = subset(seurat_object, idents = group.idents)
      Idents(seurat_object) = activeident
    }
  }
  else if (is.character(split.idents)) {
    activeident = Idents(seurat_object)
    Idents(seurat_object) = as.character(seurat_object[[split.by]][ , 1, drop = TRUE])
    seurat_object = subset(seurat_object, idents = split.idents)
    Idents(seurat_object) = activeident
  }

  if (is.null(colors)) {
    colors = hue_pal()(n = length(levels(Idents(seurat_object))))
  }

  if (isTRUE(order.colors)) {
    if (!is.null(names(colors))) {
      colors = colors[levels(Idents(seurat_object))]
    }
    if (is.character(order.prop)) {
      if (length(order.prop) > 1) {
        names(colors) = levels(Idents(seurat_object))
        colors = colors[order.prop]
        }
      else {
        if (order.prop == "reverse") {
          colors = rev(colors)
        }
      }
    }
  }

  idents.df = data.frame("ident1" = Idents(seurat_object))

  if (is.character(group.by)) {
    Idents(seurat_object) = group.by
    idents.df$ident2 = Idents(seurat_object)
  }

  table.list = proportion.plot = list()

  if (is.character(split.by)) {
    Idents(seurat_object) = split.by
    idents.df$ident3 = Idents(seurat_object)
    if (is.null(order.split)) {
      order.split = levels(Idents(seurat_object))
    }
    else {
      if (length(order.split) == length(levels(Idents(seurat_object)))) {
        seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.split)
      }
      else {
        if (length(order.split) == 1 & any(order.split == "reverse")) {
          order.split = rev(levels(Idents(seurat_object)))
          seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.split)
        }
        else {
          stop("order.split needs to be either 'reverse' or a character or numeric vector of same length as the number of identities")
        }
      }
    }
    levels.split.by = levels(Idents(seurat_object))
  for (i in levels.split.by) {
    Idents(seurat_object) = split.by
    split1.df = idents.df[which(idents.df$ident3 == i),]
    if (is.character(group.by)) {
    Idents(seurat_object) = group.by
    if (is.null(order.group)) {
      order.group = levels(Idents(seurat_object))
    }
    if (!is.null(order.group)) {
      if (length(order.group) > 1) {
        seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.group)
      }
      else {
        if (order.group == "reverse") {
          order.group = rev(levels(Idents(seurat_object)))
          seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.group)
        }
        else {
          stop("order.group needs to be either 'reverse' or a character vector")
        }
      }
    }
    k = l = 1
    table.df = sum.df = data.frame()
    for (j in levels(Idents(seurat_object))) {
      split2.df = split1.df[which(split1.df$ident2 == j),]
      if (!is.null(order.prop)) {
      if (is.character(order.prop)) {
        if (length(order.prop) > 1) {
          table.list[[i]][[j]] = as.data.frame(table(split2.df$ident1))
          colnames(table.list[[i]][[j]]) = c("ident1","nbcells")
          table.list[[i]][[j]]$percent = table.list[[i]][[j]]$nbcells/sum(table.list[[i]][[j]]$nbcells)
          table.list[[i]][[j]]$ident1 = factor(table.list[[i]][[j]]$ident1, levels = order.prop)
        }
        else {
          if (order.prop == "reverse") {
            table.list[[i]][[j]] = as.data.frame(rev(table(split2.df$ident1)))
            colnames(table.list[[i]][[j]]) = c("ident1","nbcells")
            table.list[[i]][[j]]$percent = table.list[[i]][[j]]$nbcells/sum(table.list[[i]][[j]]$nbcells)
          }
          else {
            stop("order.prop needs to be either 'reverse' or a character vector")
          }
        }
      }
        else {
          stop("order.prop needs to be either 'reverse' or a character vector")
        }
      }
      else {
        table.list[[i]][[j]] = as.data.frame(table(split2.df$ident1))
        colnames(table.list[[i]][[j]]) = c("ident1","nbcells")
        table.list[[i]][[j]]$percent = table.list[[i]][[j]]$nbcells/sum(table.list[[i]][[j]]$nbcells)
      }
      table.list[[i]][[j]]$ident2 = j
      l = l+1
      if (isTRUE(unique.group.plot & sum(table.list[[i]][[j]]$nbcells) > 0)) {
        table.df = rbind(table.df, table.list[[i]][[j]])
        k = k+1
        sum.df = rbind(sum.df,data.frame("ident1" = i,
                                         "sum" = sum(table.list[[i]][[j]]$nbcells),
                                         "ident2" = j,
                                         "sumpercent" = sum(table.list[[i]][[j]]$percent)))

      }
      else {
        if (isFALSE(unique.group.plot)) {
        if (isTRUE(prop.percent)) {
          if (isTRUE(j == levels(Idents(seurat_object))[1])) {
            if (isTRUE(show.cellsum.label)) {
                proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=percent, x=ident2, fill = ident1)) +
                  geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                     x = j,
                                                                                     y="Relative number of cells") +
                  theme_bw() +
                  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                        axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                        axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                        axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                        plot.margin = margin(5.5,5.5,5.5,20.5))+
                  scale_fill_manual(values=colors)+
                  geom_text(data = data.frame("ident1" = i,
                                              "sum" = sum(table.list[[i]][[j]]$nbcells),
                                              "ident2" = j,
                                              "sumpercent" = sum(table.list[[i]][[j]]$percent)), aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                            fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                            inherit.aes = F)+
                  scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)
            }
            else {
            proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=percent, x=ident2, fill = ident1)) +
              geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                 x = j,
                                                                                 y="Relative number of cells") +
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                    plot.margin = margin(5.5,5.5,5.5,20.5))+
              scale_y_continuous(expand= c(0,0), labels = percent)+
              scale_fill_manual(values=colors)
            }
          }
          else {
            if (isTRUE(show.cellsum.label)) {
              proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                   x = j,
                                                                                   y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_fill_manual(values=colors)+
                geom_text(data = data.frame("ident1" = i,
                                            "sum" = sum(table.list[[i]][[j]]$nbcells),
                                            "ident2" = j,
                                            "sumpercent" = sum(table.list[[i]][[j]]$percent)), aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                          fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                          inherit.aes = F)+
                scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), breaks = NULL)
            }
            else {
              proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                   x = j,
                                                                                   y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_y_continuous(expand= c(0,0), breaks = NULL)+
                scale_fill_manual(values=colors)
            }
          }
        }
        else {
          if (isTRUE(show.cellsum.label)) {
          proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=nbcells, x=ident2, fill = ident1)) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                               x = j,
                                                                               y="Number of cells") +
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side)+
            scale_fill_manual(values=colors)+
            geom_text(data = data.frame("ident1" = i,
                                        "sum" = sum(table.list[[i]][[j]]$nbcells),
                                        "ident2" = j,
                                        "sumpercent" = sum(table.list[[i]][[j]]$percent)), aes(label = sum, x = ident2, y = sum), stat = "summary",
                      fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                      inherit.aes = F)+
            scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))
        }
          else {
            proportion.plot[[i]][[j]] = ggplot(table.list[[i]][[j]], aes(y=nbcells, x=ident2, fill = ident1)) +
              geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                 x = j,
                                                                                 y="Number of cells") +
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side)+
              scale_y_continuous(expand= c(0,0))+
              scale_fill_manual(values=colors)
          }
        }
        }
      }
    }
    if (isTRUE(unique.group.plot)) {
      if (isTRUE(prop.percent)) {
        if (isTRUE(show.cellsum.label)) {
            proportion.plot[[i]] = ggplot(table.df, aes(fill=ident1, y=percent, x=factor(ident2, levels = order.group))) +
              geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                 y="Relative number of cells") +
              geom_text(data = sum.df, aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                        fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                        inherit.aes = F)+
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                    legend.position = legend.side, aspect.ratio = l/k, plot.margin = margin(5.5,65.5,5.5,5.5))+
              scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)+
              scale_fill_manual(values=colors)+
              ggtitle(i)
        }
        else {
        proportion.plot[[i]] = ggplot(table.df, aes(fill=ident1, y=percent, x=factor(ident2, levels = order.group))) +
          geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                             y="Relative number of cells") +
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                legend.position = legend.side, aspect.ratio = l/k, plot.margin = margin(5.5,65.5,5.5,5.5))+
          scale_y_continuous(expand= c(0,0), labels = percent)+
          scale_fill_manual(values=colors)+
          ggtitle(i)
        }
      }
      else {
        if (isTRUE(show.cellsum.label)) {
          proportion.plot[[i]] = ggplot(table.df, aes(fill=ident1, y=nbcells, x=factor(ident2, levels = order.group))) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                               y="Number of cells") +
            geom_text(data = sum.df, aes(label = sum, x = ident2, y = sum), stat = "summary",
                      fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                      inherit.aes = F)+
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side, aspect.ratio = l/k, plot.margin = margin(5.5,65.5,5.5,5.5))+
            scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))+
            scale_fill_manual(values=colors)+
            ggtitle(i)
        }
        else {
        proportion.plot[[i]] = ggplot(table.df, aes(fill=ident1, y=nbcells, x=factor(ident2, levels = order.group))) +
          geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                             y="Number of cells") +
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                legend.position = legend.side, aspect.ratio = l/k, plot.margin = margin(5.5,65.5,5.5,5.5))+
          scale_y_continuous(expand= c(0,0))+
          scale_fill_manual(values=colors)+
          ggtitle(i)
        }
      }
    }
    }
    else {
      if (!is.null(order.prop)) {
        if (is.character(order.prop)) {
          if (length(order.prop) > 1) {
            table.list[[i]] = as.data.frame(table(split1.df$ident1))
            colnames(table.list[[i]]) = c("ident1","nbcells")
            table.list[[i]]$percent = table.list[[i]]$nbcells/sum(table.list[[i]]$nbcells)
            table.list[[i]]$ident1 = factor(table.list[[i]]$ident1, levels = order.prop)
          }
          else {
            if (order.prop == "reverse") {
              table.list[[i]] = as.data.frame(rev(table(split1.df$ident1)))
              colnames(table.list[[i]]) = c("ident1","nbcells")
              table.list[[i]]$percent = table.list[[i]]$nbcells/sum(table.list[[i]]$nbcells)
            }
            else {
              stop("order.prop needs to be either 'reverse' or a character vector")
            }
          }
        }
        else {
          stop("order.prop needs to be either 'reverse' or a character vector")
        }
      }
      else {
        table.list[[i]] = as.data.frame(table(split1.df$ident1))
        colnames(table.list[[i]]) = c("ident1","nbcells")
        table.list[[i]]$percent = table.list[[i]]$nbcells/sum(table.list[[i]]$nbcells)
      }
      table.list[[i]]$ident3 = i

      if (isTRUE(prop.percent)) {
        if (isTRUE(show.cellsum.label)) {
          proportion.plot[[i]] = ggplot(table.list[[i]], aes(y=percent, x=ident1, fill = ident1)) +
            geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Relative number of cells") +
            geom_text(aes(label = nbcells, x = ident1, y = percent), stat = "identity",
                      vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size)+
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side, plot.margin = margin(5.5,65.5,5.5,5.5))+
            scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)+
            scale_fill_manual(values=colors)+
            NoLegend()+
            ggtitle(i)
        }
        else {
          proportion.plot[[i]] = ggplot(table.list[[i]], aes(y=percent, x=ident1, fill = ident1)) +
            geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Relative number of cells") +
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side, plot.margin = margin(5.5,65.5,5.5,5.5))+
            scale_y_continuous(expand= c(0,0), labels = percent)+
            scale_fill_manual(values=colors)+
            NoLegend()+
            ggtitle(i)
        }
      }
      else {
        if (isTRUE(show.cellsum.label)) {
          proportion.plot[[i]] = ggplot(table.list[[i]], aes(y=nbcells, x=ident1, fill = ident1)) +
            geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Number of cells") +
            geom_text(aes(label = nbcells, x = ident1, y = nbcells), stat = "identity",
                      vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size)+
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side, plot.margin = margin(5.5,65.5,5.5,5.5))+
            scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))+
            scale_fill_manual(values=colors)+
            NoLegend()+
            ggtitle(i)
        }
        else {
          proportion.plot[[i]] = ggplot(table.list[[i]], aes(y=nbcells, x=ident1, fill = ident1)) +
            geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Number of cells") +
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side, plot.margin = margin(5.5,65.5,5.5,5.5))+
            scale_y_continuous(expand= c(0,0))+
            scale_fill_manual(values=colors)+
            NoLegend()+
            ggtitle(i)
        }
      }
    }
  }

    if (isFALSE(unique.group.plot) & is.character(group.by)) {
  for (i in 1:length(table.list)) {
    k = 1
    table.tmp = table.list[[i]]
    for (j in 1:length(table.tmp)) {
      if (sum(table.tmp[[j]]$nbcells) == 0) {
        proportion.plot[[i]][[k]] = NULL
        table.list[[i]][[k]] = NULL
      }
      else
        k = k+1
    }
  }

  for (i in 1:length(proportion.plot)) {
    for (j in 2:length(proportion.plot[[i]])) {
      proportion.plot[[i]][[j]] = proportion.plot[[i]][[j]]+
        theme(axis.title.y = element_blank())
    }
  }
  if (isTRUE(unique.split.plot)) {
    Idents(seurat_object) = split.by
    temp.plot = list()
    temp.title = ""
    for (i in levels(Idents(seurat_object))) {
      temp.plot = c(temp.plot,proportion.plot[[i]])
      temp.title = paste0(temp.title,i)
      if (i != levels(Idents(seurat_object))[length(levels(Idents(seurat_object)))]) {
        temp.title = paste0(temp.title," vs ")
        temp.plot = c(temp.plot,list(plot_spacer()))
      }
    }
    temp.plot = wrap_plots(temp.plot, nrow = nrow, guides = "collect")+
      plot_annotation(title = temp.title, theme = theme(plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                                    legend.position = legend.side))
    proportion.plot = temp.plot
    if (isFALSE(show.legend)) {
      proportion.plot = proportion.plot+
        plot_annotation(theme = theme(legend.position = "none"))
    }
  }
  else {
  Idents(seurat_object) = split.by
  for (i in levels(Idents(seurat_object))) {
    proportion.plot[[i]] = wrap_plots(proportion.plot[[i]], nrow = nrow, guides = "collect")+
      plot_annotation(title = i, theme = theme(plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                                    legend.position = legend.side))
    if (isFALSE(show.legend)) {
      proportion.plot[[i]] = proportion.plot[[i]]+
        plot_annotation(theme = theme(legend.position = "none"))
    }
  }
  }
      if(isTRUE(output.data)) {
        return(table.list)
      }
      else {
        return(proportion.plot)
      }
  }

    else {
    if (isTRUE(unique.split.plot)) {
    Idents(seurat_object) = split.by
    proportion.plot[[length(levels(Idents(seurat_object)))]] = proportion.plot[[length(levels(Idents(seurat_object)))]]+
      theme(plot.margin = margin(5.5,5.5,5.5,5.5))
    proportion.plot = wrap_plots(proportion.plot, nrow = nrow, guides = "collect")+
      plot_annotation(theme = theme(legend.position = legend.side))
    if (isFALSE(show.legend)) {
      proportion.plot = proportion.plot+
        plot_annotation(theme = theme(legend.position = "none"))
    }
    }
    else {
      for (i in 1:length(proportion.plot)) {
        proportion.plot[[i]] = proportion.plot[[i]]+
          theme(plot.margin = margin(5.5,5.5,5.5,5.5))
        if (isFALSE(show.legend)) {
          proportion.plot[[i]] = proportion.plot[[i]]+
            NoLegend()
        }
      }
    }
      if(isTRUE(output.data)) {
        for (i in 1:length(table.list)) {
          k = 1
          table.tmp = table.list[[i]]
          for (j in 1:length(table.tmp)) {
            if (sum(table.tmp[[j]]$nbcells) == 0) {
              table.list[[i]][[k]] = NULL
            }
            else
              k = k+1
          }
        }
        return(table.list)
      }
      else {
        return(proportion.plot)
      }
  }
  }

  if (is.character(group.by)) {
    Idents(seurat_object) = group.by
    if (is.null(order.group)) {
      order.group = levels(Idents(seurat_object))
    }
    if (!is.null(order.group)) {
      if (length(order.group) > 1) {
        seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.group)
      }
      else {
        if (order.group == "reverse") {
          order.group = rev(levels(Idents(seurat_object)))
          seurat_object@active.ident = factor(seurat_object@active.ident, levels = order.group)
        }
        else {
          stop("order.group needs to be either 'reverse' or a character vector")
        }
      }
    }
    table.df = sum.df = data.frame()
    for (j in levels(Idents(seurat_object))) {
      split1.df = idents.df[which(idents.df$ident2 == j),]
      if (!is.null(order.prop)) {
        if (is.character(order.prop)) {
          if (length(order.prop) > 1) {
            table.list[[j]] = as.data.frame(table(split1.df$ident1))
            colnames(table.list[[j]]) = c("ident1","nbcells")
            table.list[[j]]$percent = table.list[[j]]$nbcells/sum(table.list[[j]]$nbcells)
            table.list[[j]]$ident1 = factor(table.list[[j]]$ident1, levels = order.prop)
          }
          else {
            if (order.prop == "reverse") {
              table.list[[j]] = as.data.frame(rev(table(split1.df$ident1)))
              colnames(table.list[[j]]) = c("ident1","nbcells")
              table.list[[j]]$percent = table.list[[j]]$nbcells/sum(table.list[[j]]$nbcells)
            }
            else {
              stop("order.prop needs to be either 'reverse' or a character vector")
            }
          }
        }
        else {
          stop("order.prop needs to be either 'reverse' or a character vector")
        }
      }
      else {
        table.list[[j]] = as.data.frame(table(split1.df$ident1))
        colnames(table.list[[j]]) = c("ident1","nbcells")
        table.list[[j]]$percent = table.list[[j]]$nbcells/sum(table.list[[j]]$nbcells)
      }
      table.list[[j]]$ident2 = j
      if (isTRUE(unique.group.plot & sum(table.list[[j]]$nbcells) > 0)) {
        table.df = rbind(table.df, table.list[[j]])
        sum.df = rbind(sum.df,data.frame("sum" = sum(table.list[[j]]$nbcells),
                                         "ident2" = j,
                                         "sumpercent" = sum(table.list[[j]]$percent)))

      }
      else {
        if (isTRUE(prop.percent)) {
          if (isTRUE(j == levels(Idents(seurat_object))[1])) {
            if (isTRUE(show.cellsum.label)) {
              proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                   y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_fill_manual(values=colors)+
                geom_text(data = data.frame("sum" = sum(table.list[[j]]$nbcells),
                                            "ident2" = j,
                                            "sumpercent" = sum(table.list[[j]]$percent)), aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                          fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                          inherit.aes = F)+
                scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)
            }
            else {
              proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                   y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_y_continuous(expand= c(0,0), labels = percent)+
                scale_fill_manual(values=colors)
            }
          }
          else {
            if (isTRUE(show.cellsum.label)) {
              proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                   y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_fill_manual(values=colors)+
                geom_text(data = data.frame("sum" = sum(table.list[[j]]$nbcells),
                                            "ident2" = j,
                                            "sumpercent" = sum(table.list[[j]]$percent)), aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                          fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                          inherit.aes = F)+
                scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), breaks = NULL)
            }
            else {
              proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=percent, x=ident2, fill = ident1)) +
                geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                   y="Relative number of cells") +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                      axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                      axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side,
                      plot.margin = margin(5.5,5.5,5.5,20.5))+
                scale_y_continuous(expand= c(0,0), breaks = NULL)+
                scale_fill_manual(values=colors)
            }
          }
        }
        else {
          if (isTRUE(show.cellsum.label)) {
            proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=nbcells, x=ident2, fill = ident1)) +
              geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                 y="Number of cells") +
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side)+
              scale_fill_manual(values=colors)+
              geom_text(data = data.frame("sum" = sum(table.list[[j]]$nbcells),
                                          "ident2" = j,
                                          "sumpercent" = sum(table.list[[j]]$percent)), aes(label = sum, x = ident2, y = sum), stat = "summary",
                        fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                        inherit.aes = F)+
              scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))
          }
          else {
            proportion.plot[[j]] = ggplot(table.list[[j]], aes(y=nbcells, x=ident2, fill = ident1)) +
              geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                                 y="Number of cells") +
              theme_bw() +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                    axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                    axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), legend.position = legend.side)+
              scale_y_continuous(expand= c(0,0))+
              scale_fill_manual(values=colors)
          }
        }
      }
    }
    if (isFALSE(unique.group.plot)) {
      for (i in 2:length(proportion.plot)) {
          proportion.plot[[i]] = proportion.plot[[i]]+
            theme(axis.title.y = element_blank())
      }
      proportion.plot = wrap_plots(proportion.plot, nrow = nrow, guides = "collect")+
        plot_annotation(theme = theme(legend.position = legend.side))
      if (isFALSE(show.legend)) {
        proportion.plot = proportion.plot+
          plot_annotation(theme = theme(legend.position = "none"))
      }
      if(isTRUE(output.data)) {
          k = 1
          table.tmp = table.list
          for (j in 1:length(table.tmp)) {
            if (sum(table.tmp[[j]]$nbcells) == 0) {
              table.list[[k]] = NULL
            }
            else
              k = k+1
          }
        return(table.list)
      }
      else {
        return(proportion.plot)
      }
    }
    else {
      if (isTRUE(prop.percent)) {
        if (isTRUE(show.cellsum.label)) {
          proportion.plot = ggplot(table.df, aes(fill=ident1, y=percent, x=factor(ident2, levels = order.group))) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                               y="Relative number of cells") +
            geom_text(data = sum.df, aes(label = sum, x = ident2, y = sumpercent), stat = "summary",
                      fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                      inherit.aes = F)+
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side)+
            scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)+
            scale_fill_manual(values=colors)
        }
        else {
          proportion.plot = ggplot(table.df, aes(fill=ident1, y=percent, x=factor(ident2, levels = order.group))) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                               y="Relative number of cells") +
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side)+
            scale_y_continuous(expand= c(0,0), labels = percent)+
            scale_fill_manual(values=colors)
        }
      }
      else {
        if (isTRUE(show.cellsum.label)) {
          proportion.plot = ggplot(table.df, aes(fill=ident1, y=nbcells, x=factor(ident2, levels = order.group))) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                               y="Number of cells") +
            geom_text(data = sum.df, aes(label = sum, x = ident2, y = sum), stat = "summary",
                      fun = sum, vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size,
                      inherit.aes = F)+
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side)+
            scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))+
            scale_fill_manual(values=colors)
        }
        else {
          proportion.plot = ggplot(table.df, aes(fill=ident1, y=nbcells, x=factor(ident2, levels = order.group))) +
            geom_bar(alpha = alpha, position="stack", stat="identity", color = "black")+ labs(fill="",
                                                                               y="Number of cells") +
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                  axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                  axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size), plot.title = element_text(hjust = 0.5, size = split.plot.title.size),
                  legend.position = legend.side)+
            scale_y_continuous(expand= c(0,0))+
            scale_fill_manual(values=colors)
        }
      }
      if (isFALSE(show.legend)) {
        proportion.plot = proportion.plot+
          NoLegend()
      }
      if(isTRUE(output.data)) {
        return(table.df)
      }
      else {
        return(proportion.plot)
      }
    }
  }

  else {
    if (!is.null(order.prop)) {
      if (is.character(order.prop)) {
        if (length(order.prop) > 1) {
          table.df = as.data.frame(table(idents.df$ident1))
          colnames(table.df) = c("ident1","nbcells")
          table.df$percent = table.df$nbcells/sum(table.df$nbcells)
          table.df$ident1 = factor(table.df$ident1, levels = order.prop)
        }
        else {
          if (order.prop == "reverse") {
            table.df = as.data.frame(rev(table(idents.df$ident1)))
            colnames(table.df) = c("ident1","nbcells")
            table.df$percent = table.df$nbcells/sum(table.df$nbcells)
          }
          else {
            stop("order.prop needs to be either 'reverse' or a character vector")
          }
        }
      }
      else {
        stop("order.prop needs to be either 'reverse' or a character vector")
      }
    }
    else {
      table.df = as.data.frame(table(idents.df$ident1))
      colnames(table.df) = c("ident1","nbcells")
      table.df$percent = table.df$nbcells/sum(table.df$nbcells)
    }
    sum.df = data.frame("sum" = sum(table.df$nbcells),
                        "ident1" = levels(Idents(seurat_object)),
                        "sumpercent" = sum(table.df$percent))
    if (isTRUE(prop.percent)) {
      if (isTRUE(show.cellsum.label)) {
        proportion.plot = ggplot(table.df, aes(y=percent, x=ident1, fill = ident1)) +
          geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Relative number of cells") +
          geom_text(aes(label = nbcells, x = ident1, y = percent), stat = "identity",
                    vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size)+
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size),
                legend.position = legend.side)+
          scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)), labels = percent)+
          scale_fill_manual(values=colors)+
          NoLegend()
      }
      else {
        proportion.plot = ggplot(table.df, aes(y=percent, x=ident1, fill = ident1)) +
          geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Relative number of cells") +
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size),
                legend.position = legend.side)+
          scale_y_continuous(expand= c(0,0), labels = percent)+
          scale_fill_manual(values=colors)+
          NoLegend()
      }
    }
    else {
      if (isTRUE(show.cellsum.label)) {
        proportion.plot = ggplot(table.df, aes(y=nbcells, x=ident1, fill = ident1)) +
          geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Number of cells") +
          geom_text(aes(label = nbcells, x = ident1, y = nbcells), stat = "identity",
                    vjust = -0.5, position = position_dodge(width = 1), size = cellsum.label.size)+
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size),
                legend.position = legend.side)+
          scale_y_continuous(expand= expansion(mult = c(0,cellsum.label.size/50)))+
          scale_fill_manual(values=colors)+
          NoLegend()
      }
      else {
        proportion.plot = ggplot(table.df, aes(y=nbcells, x=ident1, fill = ident1)) +
          geom_bar(alpha = alpha, position= "stack", stat = "identity", color = "black")+ labs(y="Number of cells") +
          theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x = element_blank(), axis.text.x = element_text(angle = x.axis.angle, hjust = x.axis.hjust),
                axis.title.x = element_blank(), axis.title.y = element_text(size = y.axis.title.size),
                axis.text = element_text(size = axis.text.size), legend.text = element_text(size = legend.text.size),
                legend.position = legend.side)+
          scale_y_continuous(expand= c(0,0))+
          scale_fill_manual(values=colors)+
          NoLegend()
      }
    }
    if (isFALSE(show.legend)) {
      proportion.plot = proportion.plot+
        NoLegend()
    }
    if(isTRUE(output.data)) {
      return(table.df)
    }
    else {
      return(proportion.plot)
    }
  }
}
