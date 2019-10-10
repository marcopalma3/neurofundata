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
