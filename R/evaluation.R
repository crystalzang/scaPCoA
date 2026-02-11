
#' Calculate Feature Importance from PCoA Eigenvectors
#'
#' Computes feature importance scores by calculating the covariance between
#' original features and principal coordinate axes, weighted by eigenvalues.
#'
#' @param eigen_vectors A matrix of eigenvectors from PCoA decomposition.
#' @param eigen_values A numeric vector of corresponding eigenvalues.
#' @param X_train A numeric matrix of training data (rows = samples,
#'   columns = features).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{feature}{Feature name from column names of \code{X_train}.}
#'     \item{sum}{Absolute sum of standardized covariances across all axes.}
#'     \item{sum_weighted}{Absolute eigenvalue-weighted sum of covariances.}
#'     \item{rank_sum}{Rank based on \code{sum} (1 = most important).}
#'     \item{rank_sum_weighted}{Rank based on \code{sum_weighted} (1 = most important).}
#'   }
#'   Additional columns correspond to the standardized covariance with each
#'   principal coordinate axis.
#'
#' @details
#' Standardization of covariances follows the approach described in
#' Legendre & Legendre (1998). Each covariance is scaled by the inverse
#' square root of the eigenvalue divided by \code{n - 1}.
#'
#' @references
#' Legendre, P. & Legendre, L. (1998). \emph{Numerical Ecology}. Elsevier.
#'
#' @importFrom stats cov
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @export
calc_feature_importance <- function(eigen_vectors, eigen_values, X_train){
  eigen_vectors_standard <- scale(eigen_vectors)
  # Compute covariance of variables with all axes
  S <- cov(X_train, eigen_vectors_standard)
  sum_weighted <- as.vector(t(eigen_values %*% t(S[,1:ncol(S)]))) # weighted sum (CZ)

  # Standardize value of covariance (see Legendre & Legendre 1998)
  U <- S %*% diag((eigen_values/(nrow(X_train)- 1))^(-0.5))
  colnames(U) <- colnames(eigen_vectors)


  U_df <- U %>%
    data.frame()%>%
    mutate(sum = abs(rowSums(.)),
           sum_weighted = abs(sum_weighted))%>%
    rownames_to_column(var = "feature")%>%
    mutate(
      rank_sum = rank(-sum, ties.method = "min"),  # Rank 1 = highest sum
      rank_sum_weighted = rank(-sum_weighted, ties.method = "min")  # Rank 1 = highest sum_weighted
    )
  return(U_df)
}
