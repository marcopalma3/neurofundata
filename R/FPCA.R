#' Functional principal component analysis
#'
#' Return eigenimages, eigenvectors and scores.
#'
#' @param basis_mat matrix of basis functions (each basis is vectorised in one column of the matrix)
#' @param knot_space space between consecutive knots (in mm, equal for each dimension)
#' @param pve_threshold proportion of variance explained
#' @param data_projected_name text file with the smoothing projections for each statistical unit
#' @param train_sub rownames of \code{data_projected_name} used to fit the FPCA  
#' @param test_sub rownames of \code{data_projected_name} for which to get the FPCA scores
#'
#' @return A list with the following elements
#' \describe{
#'   \item{scores}{Scores matrix (each row is for one statistical unit)}
#'   \item{eFuns}{Eigenimages (vectorised on each column)}
#'   \item{eVals}{Eigenvalues}
#' }
#'
#' @author Marco Palma, \email{M.Palma@@warwick.ac.uk}
#' @keywords FPCA
#'
#' @export
#' @importFrom data.table fread
#' @importFrom Matrix crossprod chol
#' @importFrom spam kronecker
#' @importFrom magrittr %>%



FPCA <- function(basis_mat,
                 knot_space,
                 pve_threshold,
                 data_projected_name = paste0("data_projected_", knot_space, "mm.dat"),
                 train_sub,
                 test_sub = train_sub) {

  data_projected <- data_projected_name %>%
    data.table::fread(.) %>%
    data.frame(., row.names = 1) %>%
    as.matrix()

  int_mat <- Matrix::crossprod(basis_mat)  ###matrix W
  basis_W_half <- try(Matrix::chol(int_mat))
  if("try-error" %in% class(basis_W_half)) basis_W_half <- Matrix::chol(int_mat, pivot = TRUE)

  data_projected_train <- data_projected[train_sub, ]
  #data_projected_test <- data_projected[test_sub, ]

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
