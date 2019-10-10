slices_plot <- function(image_vec,
                        dims = dims_mask,
                        voxels = results$voxel_grid_nonzero_mask,
                        col_threshold = 0,
                        legend_range = range(image_vec),
                        slice_number = FALSE, ...){

  img <- rep(NA, prod(dims))
  img[voxels] <- image_vec

  ybr <- seq(legend.range[1], legend.range[2], length.out = 101)
  rc1 <- grDevices::colorRampPalette(colors = c("blue", "white"), space="Lab")(length(which(ybr < col_threshold)))
  rc2 <- grDevices::colorRampPalette(colors = c("white", "red"), space="Lab")(length(which(ybr > col_threshold))-1)
  rampcols <- c(rc1, rc2)

  if (col_threshold == 0) {
    leg_text <- c(format(min(ybr), scientific = TRUE, digits = 3),
                  col.threshold,
                  format(max(ybr), scientific = TRUE, digits = 3))
  } else {
    leg_text <- c(floor(min(ybr)), col_threshold, ceiling(max(ybr)))
  }

  overlay(
    x = nifti(resize_image(mask, img_template)$array),
    y = nifti(array(img, dim = dims)),
    plot.type = "single",
    z = seq(16, 136, by = 5),
    col.y = scales::alpha(rampcols, 0.45),
    useRaster = TRUE,
    oma = c(4,0,0,0), ...
  )

  fields::image.plot(
    legend.only = TRUE,
    zlim = range(ybr),
    col = rampcols,
    horizontal = TRUE,
    legend.mar = 1.5,
    legend.cex = 0.5,
    legend.width = 0.8,
    legend.args = list(text = leg_text,
                       col = "white", cex = 0.7, side = 1,
                       at = c(min(ybr), col_threshold, max(ybr))
    )
  )

  if (slice_number == TRUE) {
    par(xpd = TRUE)
    lab_img <- 34 + seq(16, 136, by = 5)
    annot_grid <- expand.grid(0.025 + seq(0, 0.95, by = 0.2), seq(0.97, 0.1, by = -0.175))
    grid::grid.text(lab_img,
                    x = annot_grid[,1],
                    y = annot_grid[,2],
                    gp = gpar(fontsize=10, cex = 0.9, col = "white")
    )
  }

}
