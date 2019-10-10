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
