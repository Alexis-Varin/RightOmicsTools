#' @title Discrete palette
#'
#' @description This function generates a discrete color palette from \code{\link[grDevices]{colors}}, starting from a given color. With base parameters, this function generates 200 colors, based on their luminance and saturation, chosen to be maximally distinct.
#'
#' @param n Numeric. The number of colors to return. If \code{NULL}, returns all the colors.
#' @param starting.color Character. The name of the color to start the palette from. The subsequent maximally distinct colors will then be generated, meaning choosing different starting colors produce different palettes.
#' @param exclude.hues Character. The names of the shades of colors to exclude from the palette (for example, 'green' will exclude all shades of green).
#' @param offset Numeric. The number of colors to skip from the start of the palette.
#' @param lum.range Numeric. The range of luminance values to consider for the palette. Contrasted colors have low luminance while pastel colors have high luminance.
#' @param sat.range Numeric. The range of saturation values to consider for the palette. Washed-out colors have low saturation while vivid colors have high saturation.
#'
#' @return A character vector of color names.
#'
#' @examples
#' \dontshow{
#' suppressWarnings(suppressPackageStartupMessages(library(scales)))
#' }
#' # Example 1: default parameters
#' show_col(Right_Palette(), labels = FALSE)
#'
#' # Example 2: a subset of colors
#' show_col(Right_Palette(100), labels = FALSE)
#'
#' # Example 3: with offset, skipping the first 10 colors
#' show_col(Right_Palette(100, offset = 10), labels = FALSE)
#'
#' # Example 4: different starting color
#' show_col(Right_Palette(100, starting.color = "green"), labels = FALSE)
#'
#' # Example 5: pastel colors palette, high luminance and low saturation
#' show_col(Right_Palette(30, lum.range = c(0.7, 1), sat.range = c(0, 0.3)), labels = FALSE)
#'
#' # Example 6: contrasted colors palette, low luminance and high saturation
#' show_col(Right_Palette(30, lum.range = c(0, 0.4), sat.range = c(0.6, 1)), labels = FALSE)
#'
#' # Example 7: washed-out colors palette, low luminance and low saturation
#' show_col(Right_Palette(30, lum.range = c(0, 0.4), sat.range = c(0, 0.4)), labels = FALSE)
#'
#' # Example 8: vivid colors palette, high luminance and high saturation
#' show_col(Right_Palette(30, lum.range = c(0.5, 1), sat.range = c(0.5, 1)), labels = FALSE)
#'
#' # Example 9: excluding some shades of colors
#' show_col(Right_Palette(100, exclude.hues = "blue"), labels = FALSE)
#' show_col(Right_Palette(100, exclude.hues = "red"), labels = FALSE)
#' show_col(Right_Palette(100, exclude.hues = "green"), labels = FALSE)
#' show_col(Right_Palette(20, exclude.hues = c("red", "orange", "yellow", "cyan",
#'                                             "blue", "purple", "pink")), labels = FALSE)
#'
#' @importFrom grDevices colors col2rgb rgb2hsv
#'
#' @export

Right_Palette = function(n = NULL,
                         starting.color = "orange",
                         exclude.hues = c("white", "grey", "black"),
                         offset = 0,
                         lum.range = c(0, 0.8),
                         sat.range = c(0.5, 1)) {

  cols = colors(distinct = TRUE)

  rgb.mat = col2rgb(cols)/255

  lumfct = function(c) ifelse(c <= 0.04045, c/12.92, ((c+0.055)/1.055)^2.4)

  lum = 0.2126*lumfct(rgb.mat[1 , ])+0.7152*lumfct(rgb.mat[2 , ])+0.0722*lumfct(rgb.mat[3 , ])

  hsv.mat = rgb2hsv(rgb.mat, maxColorValue = 1)
  sat = hsv.mat["s", ]

  cols = cols[lum >= lum.range[1] & lum <= lum.range[2] & sat >= sat.range[1] & sat <= sat.range[2]]
  lum = lum[lum >= lum.range[1] & lum <= lum.range[2] & sat >= sat.range[1] & sat <= sat.range[2]]

  rgb.mat = t(col2rgb(cols))/255

  lum.mat = ifelse(rgb.mat <= 0.04045, rgb.mat/12.92, ((rgb.mat+0.055)/1.055)^2.4)

  M = matrix(c(0.4124564, 0.3575761, 0.1804375,
               0.2126729, 0.7151522, 0.0721750,
               0.0193339, 0.1191920, 0.9503041),
             nrow = 3, byrow = TRUE)

  XYZ = lum.mat%*%t(M)

  XYZ[ , 1] = XYZ[ , 1]/0.95047
  XYZ[ , 3] = XYZ[ , 3]/1.08883

  labfct = function(col) ifelse(col > (6/29)^3, col^(1/3), col/(3*(6/29)^2)+4/29)

  fX = labfct(XYZ[ , 1])
  fY = labfct(XYZ[ , 2])
  fZ = labfct(XYZ[ , 3])

  lab = cbind(L = 116*fY-16,
              a = 500*(fX-fY),
              b = 200*(fY-fZ))

  rownames(lab) = cols

  lab.hue = (atan2(lab[ , "b"], lab[ , "a"])%%(2*pi))/(2*pi)
  lab.chroma = sqrt(lab[ , "a"]^2+lab[ , "b"]^2)

  hue.bins = list(red = c(0, 45)/360,
                  orange = c(45, 75)/360,
                  yellow = c(75, 115)/360,
                  green = c(115, 165)/360,
                  cyan = c(165, 230)/360,
                  blue = c(230, 310)/360,
                  purple = c(310, 340)/360,
                  pink = c(340, 360)/360)

  achromatic.bins = list(white = c(0.8, 1),
                         grey  = c(0.05, 0.8),
                         black = c(0, 0.05))

  hue.excluded = rep(FALSE, nrow(lab))
  if (is.character(exclude.hues)) {
    for (hues in exclude.hues[exclude.hues %in% names(hue.bins)]) {
      hue.excluded = hue.excluded | (lab.chroma > 5 & lab.hue >= hue.bins[[hues]][1] & lab.hue <= hue.bins[[hues]][2])
    }
    for (hues in exclude.hues[exclude.hues %in% names(achromatic.bins)]) {
      hue.excluded = hue.excluded | (lab.chroma <= 5 & lum >= achromatic.bins[[hues]][1] & lum <= achromatic.bins[[hues]][2])
    }
  }

  lab = lab[!hue.excluded, ]

  if (isFALSE(starting.color %in% rownames(lab))) {
    message("'", starting.color, "' is not present with the current parameters, using '", rownames(lab)[1], "' instead", call. = FALSE)
    starting.color = rownames(lab)[1]
  }

  chosen = starting.color
  remaining = setdiff(rownames(lab), chosen)

  dist.min = sqrt(rowSums((lab[remaining, , drop = FALSE]-matrix(lab[starting.color, ], length(remaining), 3, byrow = TRUE))^2))
  names(dist.min) = remaining

  while (length(remaining) > 1) {
    nxt = names(which.max(dist.min))
    chosen = c(chosen, nxt)
    remaining = setdiff(remaining, nxt)
    dist.new = sqrt(rowSums((lab[remaining, , drop = FALSE]-matrix(lab[nxt, ], length(remaining), 3, byrow = TRUE))^2))
    dist.min = pmin(dist.min[remaining], dist.new)
  }

  if (is.null(n)) n = length(chosen)

  if (n > length(chosen)) {
    warning("Only ", length(chosen), " distinct colors available, returning all", call. = FALSE)
    return(chosen)
  }
  else if (n+offset > length(chosen)) {
    warning("The number of colors asked with offset is higher than the number of distinct colors available,\ndecreasing offset to ", length(chosen)-n,  call. = FALSE)
    return(chosen[(1+length(chosen)-n):(length(chosen))])
  }
  else return(chosen[(1+offset):(n+offset)])
}
