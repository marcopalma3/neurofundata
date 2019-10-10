#' Colour palette
#'
#' Create colour palette with equidistant hues.
#'
#' @param n number of colours
#'
#' @author Marco Palma, \email{M.Palma@@warwick.ac.uk}
#' @keywords gg_color_hue
#'
#' @examples
#' gg_color_hue(3)
#'
#' @export


gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
