#' Eigen Decomposition of Bray-Curtis Dissimilarity Matrix
#'
#' Computes the Bray-Curtis dissimilarity matrix, performs double centering,
#' and returns eigenvalues and eigenvectors via Principal Coordinates Analysis (PCoA).
#'
#' @param X A numeric matrix of species abundances where rows are samples and
#'   columns are taxa.
#' @param pos Logical. If \code{TRUE}, return only eigenvectors and eigenvalues
#'   corresponding to positive eigenvalues. If \code{FALSE}, return all
#'   eigenvectors and eigenvalues.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{eigenvectors}{A matrix of eigenvectors (columns correspond to
#'       principal coordinates).}
#'     \item{eigenvalues}{A numeric vector of eigenvalues.}
#'   }
#'
#'
#' @importFrom vegan vegdist
#' @export

eigen_bray <- function(X, pos){
  # Step 1. Calculate bray-curtis dissimilarity
  D <- vegdist(X, method = "bray") # bray-curtis distance
  D <- as.matrix(D)

  # Step 2: center and square D
  S <- -0.5 * (D^2)

  # Step 3: Double center S
  N <- nrow(S)
  one_vector <- matrix(1, nrow = N, ncol = 1)
  k <- one_vector / N
  I <- diag(N)
  A_normalized <- (I - k %*% t(one_vector)) %*% S %*% (I - one_vector %*% t(k))

  # Step 4: Eigen Decomposition
  eigen_decomp <- eigen(A_normalized)

  if(pos==TRUE){
    # ACPCOA used this epsilon. To be consistent when adapting their function, we set the same epsilon.
    # https://github.com/YuWang28/acPCoA/blob/master/R/acPCoA.R
    epsilon <- sqrt(.Machine$double.eps)*10000
    positive_indices <- which(eigen_decomp$values > epsilon)
    if (length(positive_indices) == 0){
      stop("No positive eigenvalues found.")
    }
    # Step 5: Obtain reduced dimension X_hat with dimension of positive eigenvalues

    eigenvalues <- eigen_decomp$values[positive_indices]
    eigenvectors <- eigen_decomp$vectors[, positive_indices]
    print(paste("Number of features in X_hat:", dim(eigenvectors)[2]))

    return(list(eigenvectors = eigenvectors, eigenvalues = eigenvalues))
  }else{
    # return all eigenvectors and eigenvalues
    eigenvalues <- eigen_decomp$values
    eigenvectors <- eigen_decomp$vectors
    print(paste("Number of features in X_hat:", dim(eigenvectors)[2]))

    return(list(eigenvectors = eigenvectors, eigenvalues = eigenvalues))
  }
}
