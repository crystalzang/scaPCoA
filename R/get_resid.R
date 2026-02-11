#' Compute Residual Matrices Adjusting for Outcome and Covariates
#'
#' Regresses each column of X on Y and A, then computes two residual matrices:
#' X.tilde (removing the effect of Y) and X.star (removing the effect of A).
#'
#' @param X A numeric matrix of features (rows = samples, columns = features).
#' @param Y A numeric vector or matrix of outcome values.
#' @param A A numeric matrix of covariates to adjust for.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{X.tilde}{Residual matrix after removing the effect of Y from X.}
#'     \item{X.star}{Residual matrix after removing the effect of A from X.}
#'   }
#'
#' @details
#' For each feature j in X, a linear model \code{X[,j] ~ Y + A} is fit.
#' X.tilde is obtained by subtracting the Y component, and X.star is obtained
#' by subtracting the A component. These residual matrices are used in
#' the scaPCoA-R (residualization) variant.
#'
#' @importFrom stats lm coef
#' @export
get.resid = function(X, Y, A){
  # Obtain residual X: X.tidle and X.star
  X = as.matrix(X) # make sure it's a matrix
  A = as.matrix(A)

  X.tilde = c()
  X.star = c()
  for (j in 1:ncol(X)){
    lm.X.j = lm(X[,j] ~ Y + A)
    X.tilde = cbind(X.tilde, (X[,j] - lm.X.j$coefficients["Y"]*Y))
    X.star = cbind(X.star, (X[,j] - A%*%lm.X.j$coefficients[startsWith(names(lm.X.j$coefficients),"A")]))
  }

  return(list(X.tilde = X.tilde,
              X.star = X.star))
}
