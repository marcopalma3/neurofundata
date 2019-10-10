#' Prediction intervals for functional quantile regression
#'
#' Return point predictions and prediction intervals for functional 
#' scalar-on-image quantile regression and outputs from FPCA.
#'
#' @param pred_interval rownames of \code{data_projected_name} used to fit the FPCA  
#' @param col rownames of \code{data_projected_name} for which to get the FPCA scores
#' @param xlab_input text file with the smoothing projections for each statistical unit
#' @param ylab_input demographic data with the scalar outcome of interest
#'
#' @author Marco Palma, \email{M.Palma@@warwick.ac.uk}
#' @keywords FPCA
#'
#' @export
#' @importFrom data.table fread
#' @importFrom Matrix crossprod chol
#' @importFrom spam kronecker


predint_plot <- function(pred_interval,
                         col = c("#009E73", "gold3", "#D55E00"),
                         xlab_input = "Difference from chronological age",
                         ylab_input = "Subjects") {

  data_forest_center <- pred_interval$data_forest_center
  #col <- mycol[as.numeric(data_forest_center$Dx)]
  excesspoints <- pred_interval$excesspoints
  #pal <- c("#009E73", "gold3", "#D55E00")

  ggplot(data_forest_center,
         aes(y = id,
             x = AgePredMed,
             xmin = AgePredLower,
             xmax = AgePredUpper,
             colour = Dx))+
    geom_point(cex = 0.5)+
    scale_colour_manual(values = col) +
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
