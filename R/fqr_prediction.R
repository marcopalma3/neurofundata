#' Prediction intervals for functional quantile regression
#'
#' Return point predictions and prediction intervals for functional 
#' scalar-on-image quantile regression and outputs from FPCA.
#'
#' @param train_sub rownames of \code{data_projected_name} used to fit the FPCA  
#' @param test_sub rownames of \code{data_projected_name} for which to get the FPCA scores
#' @param data_projected_name text file with the smoothing projections for each statistical unit
#' @param data_demo demographic data with the scalar outcome of interest
#' @param pred_table table with predictions for each statistical unit
#' @param qr_pen penalty type for quantile regression
#' @param qr_postLASSO refit quantile regression with LASSO selected variables? 
#' @param lambda lambda parameter for penalised quantile regression
#' @param tau_levels quantiles for which quantile regression is fitted
#' @param return_func_coef return the functional coefficient?
#'
#' @return A list with the following elements
#' \describe{
#'   \item{pred_table}{table with predictions for each statistical unit}
#'   \item{no_fpc}{number of functional principal components selected}
#'   \item{evec}{Eigenimages}
#'   \item{model_med_coef}{Regression coefficient from median regression}
#'   \item{model_lower_coef}{Regression coefficient from regression with lower quantile}
#'   \item{model_upper_coef}{Regression coefficient from regression with lower quantile}
#'   \item{lambda_min}{Value of \code{lambda} for which \code{model_med_coef} is extracted}
#' }
#'
#' @author Marco Palma, \email{M.Palma@@warwick.ac.uk}
#' @keywords FPCA
#'
#' @export
#' @importFrom data.table fread
#' @importFrom Matrix crossprod chol
#' @importFrom spam kronecker


fqr_prediction <-  function(data_projected_name,
                            train_sub,
                            test_sub,
                            data_demo,
                            pve_threshold,
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
    dplyr::right_join(data_demo, ., by =  "PTID") %>%
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

  if (return_func_coef == TRUE) print(rqPen::cv_plots(model_med, logLambda = TRUE, loi = NULL))

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
        lambda_min = lambda_min_med
      )
    )
}
