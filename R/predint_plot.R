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
