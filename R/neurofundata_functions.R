resize_image <- function(mask, img = mask){
  mask <- drop(mask)
  img <- drop(img)
  dims <- dim(mask)

  if(!all.equal(dim(img), dims)) stop("Mask and img must have the same dimensions!")

  nonzero_coord <- matrix(NA, nrow = length(dims), ncol = 2)
  rownames(nonzero_coord) <- paste0("Dim", 1:length(dims))
  colnames(nonzero_coord) <- c("first","last")


  for(i in 1:length(dims)) nonzero_coord[i,] <- range(which(apply(mask,i,sum) != 0))
  ###select first and last coordinate for which the sum is non zero

  resized_array <- (img[,,] * (mask[,,]>0))[nonzero_coord[1,"first"]:nonzero_coord[1,"last"],
                                            nonzero_coord[2,"first"]:nonzero_coord[2,"last"],
                                            nonzero_coord[3,"first"]:nonzero_coord[3,"last"]]

  list(array = resized_array, original_coord = nonzero_coord)
}



calc_projection_Bspl <- function(knot_space, mask_fname){

  start_time <- Sys.time()
  mask <- oro.nifti::readNIfTI(fname = mask_fname)  %>%
    magrittr::is_greater_than(0.5) * 1      ### to avoid fuzzy boundaries
  mask_subset <- resize_image(mask)
  dims_mask <- dim(mask_subset$array)
  voxel_grid_nonzero_mask<- which(as.vector(mask_subset$array) != 0)

  basismat_dim1 <- fda::bsplineS(x = 1:dims_mask[1],
                                 breaks = c(seq(1,dims_mask[1], by = knot_space),
                                            dims_mask[1]),
                                 norder=3, nderiv=0, returnMatrix=TRUE)
  basismat_dim2 <- fda::bsplineS(x = 1:dims_mask[2],
                                 breaks = c(seq(1,dims_mask[2], by = knot_space),
                                            dims_mask[2]),
                                 norder=3, nderiv=0, returnMatrix=TRUE)
  basismat_dim3 <- fda::bsplineS(x = 1:dims_mask[3],
                                 breaks = c(seq(1,dims_mask[3], by = knot_space),
                                            dims_mask[3]),
                                 norder=3, nderiv=0, returnMatrix=TRUE)

  design_mat <- spam::kronecker(basismat_dim3, basismat_dim2, make.dimnames = TRUE) %>%
    spam::kronecker(., basismat_dim1, make.dimnames = TRUE) %>%
    .[voxel_grid_nonzero_mask, ]  %>%
    .[, colSums(.) != 0]

  #design_mat@x[design_mat@x<0.001] <- 0
  #design_mat <- design_mat[, colSums(design_mat) != 0]

  cat("The number of basis functions is ", ncol(design_mat),".\n", sep = "")

  list("basis_mat" = design_mat,
       "voxel_grid_nonzero_mask" = voxel_grid_nonzero_mask,
       "dims_mask" = dims_mask,
       "mask_subset" = mask_subset)
}





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



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}




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


FPCA <- function(basis_mat,
                 knot_space,
                 data_projected_name = paste0("data_projected_", knot_space, "mm.dat"),
                 train_sub,
                 test_sub) {

  data_projected <- data_projected_name %>%
    data.table::fread(.) %>%
    data.frame(., row.names = 1) %>%
    as.matrix()

  int_mat <- Matrix::crossprod(basis_mat)  ###matrix W
  basis_W_half <- try(Matrix::chol(int_mat))
  if("try-error" %in% class(basis_W_half)) basis_W_half <- Matrix::chol(int_mat, pivot = TRUE)

  data_projected_train <- data_projected[train_sub, ]
  data_projected_test <- data_projected[test_sub, ]

  coef_data <- data_projected_train %>% scale(scale = F)
  FPCA_mat <- Matrix::tcrossprod(coef_data, basis_W_half) / sqrt(nrow(coef_data) - 1)

  tmpSVD <- svd(FPCA_mat, nu = 0)  ###, nv = 50
  values <- tmpSVD$d ^ 2

  ncomp <- which(cumsum(values) / sum(values) >= pve_threshold)[1]
  vectors <- Matrix::solve(basis_W_half, tmpSVD$v[, 1:ncomp])
  eigenFuns <- basis_mat %*% -vectors

  scores_mat <- data_projected %>%
    data.matrix(., rownames.force = TRUE) %>%
    sweep(., 2, colMeans(data_projected_train)) %>%
    magrittr::multiply_by_matrix(., int_mat) %*% vectors %>%
    as.matrix()

  return(scores = scores_mat,
         eFuns = eigenFuns,
         eVals = values)
}



fqr_prediction <-  function(train_sub,
                            test_sub,
                            data_demo,
                            pred_table,
                            qr_pen = "LASSO",
                            qr_postLASSO = FALSE,
                            lambda = NULL,
                            tau_levels,
                            return_func_coef = FALSE){

  data_projected_train <- data_projected[train_sub, ]
  data_projected_test <- data_projected[test_sub, ]

  coef_data <- data_projected_train %>% scale(scale = F)
  FPCA_mat <- Matrix::tcrossprod(coef_data, basis_W_half) / sqrt(nrow(coef_data) - 1)

  tmpSVD <- svd(FPCA_mat, nu = 0)  ###, nv = 50
  values <- tmpSVD$d ^ 2

  ncomp <- which(cumsum(values) / sum(values) >= pve_threshold)[1]
  vectors <- Matrix::solve(basis_W_half, tmpSVD$v[, 1:ncomp])

  scores_mat <- data_projected %>%
    data.matrix(., rownames.force = TRUE) %>%
    sweep(., 2, colMeans(data_projected_train)) %>%
    magrittr::multiply_by_matrix(., int_mat) %*% vectors %>%
    as.matrix()

  regr_data <- scores_mat %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(., var = "PTID") %>%
    right_join(data_demo, ., by =  "PTID") %>%
    select(PTID, AGE, starts_with('V', ignore.case = FALSE)) %>%
    tibble::column_to_rownames(., "PTID")

  set.seed(657)
  model_lower <- rqPen::cv.rq.pen(x = as.matrix(regr_data[train_sub, -1]),
                                  y = regr_data[train_sub, ]$AGE,
                                  tau = tau_levels[1],
                                  method = "br",
                                  penalty = qr_pen
  )
  model_med <- rqPen::cv.rq.pen(x = as.matrix(regr_data[train_sub, -1]),
                                y = regr_data[train_sub, ]$AGE,
                                tau = tau_levels[2],
                                method = "br",
                                penalty = qr_pen
  )
  model_upper <- rqPen::cv.rq.pen(x = as.matrix(regr_data[train_sub, -1]),
                                  y = regr_data[train_sub,]$AGE,
                                  tau = tau_levels[3],
                                  method = "br",
                                  penalty = qr_pen)

  model_lower_coef <- with(model_lower, models[[which(cv[, 1] == lambda.min)]]$coefficients[-1])
  model_med_coef <- with(model_med, models[[which(cv[, 1] == lambda.min)]]$coefficients[-1])
  model_upper_coef <- with(model_upper, models[[which(cv[, 1] == lambda.min)]]$coefficients[-1])

  lambda_min_med <- model_med$lambda.min

  if (return_func_coef == TRUE) print(cv_plots(model_med, logLambda = TRUE, loi = NULL))

  if (qr_postLASSO == TRUE) {
    model_lower <- model_lower_coef %>%
      subset(., !not(.)) %>%
      purrr::when(length(names(.)) < 1  ~ 1, ~ names(.)) %>%
      paste(., collapse = "+") %>%
      paste("AGE ~ ", .) %>%
      formula(.) %>%
      rq(., tau = tau_lev[1], data = regr_data[train_sub, ])

    model_med <- model_med_coef %>%
      subset(., !not(.)) %>%
      when(length(names(.)) < 1  ~ 1, ~ names(.))  %>%
      paste(., collapse = "+") %>%
      paste("AGE ~ ", .) %>%
      formula(.) %>%
      rq(., tau = tau_lev[2], data = regr_data[train_sub, ])

    model_upper <- model_upper_coef %>%
      subset(., !not(.)) %>%
      when(length(names(.)) < 1  ~ 1, ~ names(.))  %>%
      paste(., collapse = "+") %>%
      paste("AGE ~ ", .) %>%
      formula(.) %>%
      rq(., tau = tau_lev[3], data = regr_data[train_sub, ])

    pred_table[test_sub, 1] <- model_lower %>%
      predict(newdata = regr_data[test_sub, ])
    pred_table[test_sub, 2] <- model_med %>%
      predict(newdata = regr_data[test_sub, ])
    pred_table[test_sub, 3] <- model_upper %>%
      predict(newdata = regr_data[test_sub, ])
  } else {
    pred_table[test_sub, 1] <- model_lower %>%
      predict(newx = as.matrix(regr_data[test_sub, -1]))
    pred_table[test_sub, 2] <- model_med %>%
      predict(newx = as.matrix(regr_data[test_sub, -1]))
    pred_table[test_sub, 3] <- model_upper %>%
      predict(newx = as.matrix(regr_data[test_sub, -1]))
  }
  if (return_func_coef == FALSE)
    return(pred_table)
  else
    return(
      list(
        pred_table = pred_table,
        no_fpc = ncomp,
        evec = vectors,
        model_med_coef = as(model_med_coef, "sparseMatrix"),
        model_lower_coef = as(model_lower_coef, "sparseMatrix"),
        model_upper_coef = as(model_upper_coef, "sparseMatrix"),
        lambda_min = lambda_min_med,
        no_comp = ncomp
      )
    )
}
