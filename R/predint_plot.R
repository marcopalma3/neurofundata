#' Illustration of crayon colors
#'
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param method2order method to order colors (\code{"hsv"} or \code{"cluster"})
#' @param cex character expansion for the text
#' @param mar margin parameters; vector of length 4 (see \code{\link[graphics]{par}})
#'
#' @return None
#'
#' @author Karl W Broman, \email{broman@@wisc.edu}
#' @references \url{http://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#'
#' @examples
#' plot_crayons()
#'
#' @export
#' @importFrom grDevices rgb2hsv
#' @importFrom graphics par plot rect text
#'
predint_plot <- function(pred_interval,
                         mycol = c("#009E73", "gold3", "#D55E00"),
                         xlab_input = "Difference from chronological age",
                         ylab_input = "Subjects") {

  data_forest_center <- pred_interval$data_forest_center
  #mycol <- mycol[as.numeric(data_forest_center$Dx)]
  excesspoints <- pred_interval$excesspoints
  #pal <- c("#009E73", "gold3", "#D55E00")

  ggplot(data_forest_center,
         aes(y = id,
             x = AgePredMed,
             xmin = AgePredLower,
             xmax = AgePredUpper,
             colour = Dx))+
    geom_point(cex = 0.5)+
    scale_colour_manual(values = mycol) +
    ggExtra::removeGridY() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_line(linetype = "blank"),
          strip.text = element_text(size = 12))+
    geom_errorbarh(height = 0.01) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 12)) +
    geom_point(data = excesspoints,
               mapping = aes_string(x = "x",y = "id"),
               inherit.aes = F, size = 1, shape = 18) +
    facet_wrap(. ~ Dx, scales = "free_y", shrink = TRUE) +
    xlab(xlab_input) +
    ylab(ylab_input)
}
