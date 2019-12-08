#' Bspline projection of 3D images
#'
#' Create a design matrix from an isotropic tensor product of Bsplines only 
#' for voxels within a mask. 
#'
#' @param knot_space space between consecutive knots (in mm, equal for each dimension)
#' @param mask_fname path of the mask file
#'
#' @return A list with the following elements
#' \describe{
#'   \item{basis_mat}{Matrix of basis functions (each basis is vectorised in one column of the matrix)}
#'   \item{voxel_grid_nonzero_mask}{indices of voxels within the mask}
#'   \item{dims_mask}{dimensions of the mask}
#'   \item{mask_subset}{resized mask}
#' }
#'
#' @author Marco Palma, \email{M.Palma@@warwick.ac.uk}
#' @keywords calc_projection_Bspl
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom oro.nifti readNIfTI
#' @importFrom fda bsplineS
#' @importFrom spam kronecker
#' @importFrom Matrix colSums



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
    .[, Matrix::colSums(.) != 0]

  #design_mat@x[design_mat@x<0.0001] <- 0
  #design_mat <- design_mat[, colSums(design_mat) != 0]

  cat("The number of basis functions is ", ncol(design_mat),".\n", sep = "")

  list("basis_mat" = design_mat,
       "voxel_grid_nonzero_mask" = voxel_grid_nonzero_mask,
       "dims_mask" = dims_mask,
       "mask_subset" = mask_subset)
}
