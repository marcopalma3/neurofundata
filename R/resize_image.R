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
