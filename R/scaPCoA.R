
#' PLS Maximization Function
#'
#' Computes the matrix-vector product for partial least squares eigenvalue
#' problem. Called internally by \code{\link{pedecure}}.
#'
#' @param v A numeric vector.
#' @param args A list containing \code{X} (feature matrix) and \code{Y}
#'   (outcome vector).
#'
#' @return A numeric vector resulting from \code{X'YY'X * v}.
#'
#' @keywords internal
plsMax = function(v,args){
  X = args$X
  Y = args$Y

  out = tcrossprod(crossprod(X,Y))%*%matrix(v,ncol=1)
  # Replace NA values with 0
  out[is.na(out)] <- 0

  return(out)
}

#' PeDecURe Maximization Function
#'
#' Computes the matrix-vector product for the penalized decomposition
#' eigenvalue problem. Called internally by \code{\link{pedecure}}.
#'
#' @param v A numeric vector.
#' @param args A list containing \code{X} (feature matrix),
#'   \code{X.penalize} (penalization matrix), \code{Y} (outcome),
#'   \code{A} (covariates), and \code{lambda} (tuning parameter).
#'
#' @return A numeric vector resulting from the penalized objective function
#'
#' @keywords internal
pedecureMax = function(v,args){
  # maximization function for PeDecURe (called within function pedecure()).
  # output: v which maximizes objective function
  X = args$X
  X.penalize = args$X.penalize
  Y = args$Y
  A = args$A
  lambda = args$lambda # tuning parameter lambda

  if (is.null(A)) stop("Error: A is NULL")

  XTY=crossprod(X,Y)
  K=A%*%t(A)


  out = XTY%*%t(XTY)%*%v - lambda*t(X.penalize)%*%K%*%X.penalize%*%v # based on output from calAv function in acPCA R package (Lin et al., 2016)

  # Replace NA values with 0
  out[is.na(out)] <- 0


  return(out)
}

#' Penalized Decomposition for Confounding Removal (PeDecURe)
#'
#' Performs penalized eigen decomposition that maximizes association with the
#' outcome Y while penalizing association with covariates A.
#'
#' @param X A numeric matrix of features (rows = samples, columns = features).
#' @param X.penalize A numeric matrix used for the penalization term. Must have
#'   the same number of columns as \code{X}. Can be \code{NULL} if
#'   \code{lambda = 0}.
#' @param A A numeric matrix of covariates to adjust for.
#' @param Y A numeric vector or matrix of outcome values.
#' @param lambda Numeric. Penalization tuning parameter. When \code{lambda = 0},
#'   no penalization is applied (equivalent to PLS).
#' @param nPC Integer. Number of principal components to compute.
#' @param centerX Logical. Whether to center X. Default is \code{FALSE}.
#' @param centerA Logical. Whether to center A. Default is \code{FALSE}.
#' @param centerY Logical. Whether to center Y. Default is \code{FALSE}.
#' @param scaleX Logical. Whether to scale X. Default is \code{FALSE}.
#' @param scaleA Logical. Whether to scale A. Default is \code{FALSE}.
#' @param scaleY Logical. Whether to scale Y. Default is \code{FALSE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{values}{A numeric vector of eigenvalues.}
#'     \item{vectors}{A matrix of eigenvectors.}
#'   }
#'   Returns \code{NULL} if the penalized matrix is not symmetric.
#'
#'
#' @importFrom RSpectra eigs_sym
#' @export
pedecure <- function(X, X.penalize = NULL, A, Y, lambda, nPC, centerX = FALSE, centerA = FALSE, centerY = FALSE, scaleX = FALSE, scaleA = FALSE, scaleY = FALSE) {
  # Center and scale input matrices
  A <- as.matrix(A)
  X <- as.matrix(X)
  X <- scale(X, center = centerX, scale = scaleX)


  if (!is.null(X.penalize) & lambda> 0) { # check for dimension
    if (ncol(X) != ncol(X.penalize)) {
      stop("X and X.penalize should have the same number of columns")
    }


    X.penalize <- as.matrix(X.penalize)
    X.penalize <- scale(X.penalize, center = centerX, scale = scaleX)


    args <- list(X = X, X.penalize = X.penalize, Y = Y, A = A, lambda = lambda)


    X = args$X
    X.penalize = args$X.penalize
    Y = args$Y
    A = args$A
    A <- as.matrix(A)
    lambda= args$lambda

    K <- A %*% t(A)  # Kernel matrix for penalization
    XTY <- crossprod(X, Y)  # Compute X'Y

    penalized_matrix <- XTY %*% t(XTY) - lambda* t(X.penalize) %*% K %*% X.penalize
    penalized_matrix[is.na(penalized_matrix)] <- 0  # Replace NAs with 0

    is_symmetric <- function(mat) all(abs(mat - t(mat)) < 1e-10)
    #print(paste("Penalized matrix is symmetric:",is_symmetric(penalized_matrix)))

    if (!is_symmetric(penalized_matrix)) {
      message(paste("Skipping lambda =", lambda, "because penalized_matrix is not symmetric"))
      return(NULL)
    }

    # Handle eigenvalue computation
    if (nPC == ncol(X)) {
      print("Full spectrum computation using eigen()")
      # Full spectrum computation using eigen()

      eigen_out <- eigen(penalized_matrix, symmetric = TRUE)

      eig_out <- list(values = eigen_out$values, vectors = eigen_out$vectors)
    } else {
      # Use eigs_sym for partial spectrum computation
      eig_out <- eigs_sym(
        pedecureMax,
        k = nPC,
        which = "LA",
        n = ncol(X),
        args = args,
        opts = list(tol = 1e-10, maxitr = 1e6)
      )
      #niter = eig_out$niter
      #nops = eig_out$nops
      if(all(eig_out$values == 0)){
        message("Warning: All eigenvalues are zero, returning NA")
      }
    }
  } else {
    # No penalty term (maximize X'YY'X)
    args <- list(X = X, Y = Y)

    if (nPC == ncol(X)) {
      # Full spectrum computation using eigen()
      X = args$X
      Y = args$Y

      XTY <- crossprod(X, Y)  # Compute X'Y
      pls_matrix <- XTY %*% t(XTY)  # Compute full matrix

      eigen_out <- eigen(pls_matrix, symmetric = TRUE)
      eig_out <- list(values = eigen_out$values, vectors = eigen_out$vectors)
    }  else {
      # Use eigs_sym for partial spectrum computation
      eig_out <- eigs_sym(
        plsMax,
        k = nPC,
        which = "LA", #largest eigenvalues (positive)
        n = ncol(X),
        args = args,
        opts = list(tol = 1e-10, maxitr = 1e6)
      )

      #niter = eig_out$niter
      #nops = eig_out$nops

    }
  }
  return(eig_out)
  # return(list(eig_out=eig_out,
  #             niter = niter,
  #             nops = nops))
}


#' Tune Lambda for PeDecURe
#'
#' Selects the optimal penalization parameter lambda by maximizing the
#' eigenvalue-weighted difference between partial correlations with the
#' outcome Y and partial correlations with covariates A.
#'
#' @param X.orig A numeric matrix of original features.
#' @param X.max A numeric matrix used for the maximization term.
#' @param X.penalize A numeric matrix used for the penalization term.
#' @param lambdas A numeric vector of candidate lambda values to evaluate.
#' @param A A numeric matrix of covariates.
#' @param Y A numeric vector or matrix of outcome values.
#' @param nPC Integer. Number of principal components to compute.
#' @param centerX Logical. Whether to center X. Default is \code{FALSE}.
#' @param scaleX Logical. Whether to scale X. Default is \code{FALSE}.
#' @param centerY Logical. Whether to center Y. Default is \code{FALSE}.
#' @param scaleY Logical. Whether to scale Y. Default is \code{FALSE}.
#' @param centerA Logical. Whether to center A. Default is \code{FALSE}.
#' @param scaleA Logical. Whether to scale A. Default is \code{FALSE}.
#' @param n.cores Integer. Number of cores for parallel computation.
#'   Default is 1.
#' @param plot Logical. Currently unused. Default is \code{FALSE}.
#' @param thresh Numeric. Currently unused. Default is 0.001.
#'
#' @return A list with:
#'   \describe{
#'     \item{weighted_sum_abs_diff}{Numeric vector of tuning criterion values
#'       for each valid lambda.}
#'     \item{lambdas}{The candidate lambda values evaluated.}
#'     \item{lambda_tune}{The optimal lambda that maximizes the tuning
#'       criterion. \code{NA} if no valid lambda was found.}
#'   }
#'
#' @details
#' For each candidate lambda, the function runs \code{\link{pedecure}},
#' computes in-sample PC scores, and evaluates partial correlations via
#' \code{\link{partial.cor}}. The tuning criterion is the mean across
#' covariates of the eigenvalue-weighted difference:
#' \code{|pcor(Y, PC)| - |pcor(A, PC)|}.
#'
#' @importFrom parallel mclapply
#' @export
pedecure.tune = function(X.orig, X.max, X.penalize, lambdas, A, Y, nPC, centerX = F, scaleX = F, centerY = F, scaleY = F, centerA = F, scaleA = F, n.cores = 1, plot = FALSE, thresh = 0.001){

  sum_abs_diff = numeric(length = length(lambdas))
  names_A = colnames(A)

  #lambdas = sort(lambdas)  # sort in ascending order

  best_lambda = best_weighted_sum_abs_diff = numeric(length = length(lambdas))
  all_weighted_sum_abs_diff = c()
  valid_lambdas = c()  # Track only lambdas that contribute to all_weighted_sum_abs_diff

  # Iterate through all lambdas without splitting into groups
  for (lambda_index in 1:length(lambdas)) {
    lambda_i = lambdas[lambda_index]
    print(paste("lambda=", lambda_i))

    weighted_sum_abs_diff = unlist(mclapply(lambda_i, FUN = function(lambda_i) {
      tryCatch({

        # Penalized decomposition using residuals for a given lambda:
        temp_out = pedecure(X = X.max, X.penalize = X.penalize, A = A, Y = Y, lambda= lambda_i, nPC = nPC,
                            centerX = centerX, scaleX = scaleX, centerA = centerA, scaleA = scaleA,
                            centerY = centerY, scaleY = scaleY)

        if (is.null(temp_out)){
          return(NA)
        }
        # In-sample scores for lambda = l:
        temp_scores = X.orig %*% temp_out$vectors

        # Partial correlations:
        temp_pcor = partial.cor(scores = temp_scores, A = A, Y = Y)$partial$estimates

        # Weighted sum of difference of absolute values of partial correlations
        # weighted_sum_abs_diff_pcor = sapply(1:ncol(A), FUN = function(a) {
        #   (temp_out$values / sum(temp_out$values)) %*% (abs(temp_pcor["Y", ]) - abs(temp_pcor[names_A[a], ]))
        # })

        weighted_sum_abs_diff_pcor = sapply(1:ncol(A), FUN = function(a) {
          if (sum(temp_out$values) == 0) {
            return(0)  # Avoid division by zero, return NA
          } else {
            return((temp_out$values / sum(temp_out$values)) %*%
                     (abs(temp_pcor["Y", ]) - abs(temp_pcor[names_A[a], ])))
          }
        })

        return(mean(weighted_sum_abs_diff_pcor))
      }, error = function(e){
        message(paste("Error for lambda =", lambda_i, ":", e$message))
        return(NA)  # Return NA if there's an error
      })

    },mc.cores = n.cores))

    # Store results for each lambda, only if values are valid
    if (!all(is.na(weighted_sum_abs_diff))) {
      all_weighted_sum_abs_diff = c(all_weighted_sum_abs_diff, weighted_sum_abs_diff)
      valid_lambdas = c(valid_lambdas, lambda_i)
    }
  }
  df <- data.frame(all_weighted_sum_abs_diff,valid_lambdas)
  # Check if df is empty before proceeding
  if (nrow(df) > 0) {
    max_index <- which.max(df$all_weighted_sum_abs_diff)
    lambda_max <- df$valid_lambdas[max_index]
  } else {
    lambda_max <- NA
    message("No valid lambda tuning values were found.")
  }

  return(list(weighted_sum_abs_diff = all_weighted_sum_abs_diff,
              lambdas = lambdas[1:length(all_weighted_sum_abs_diff)],
              lambda_tune = lambda_max))
}


#' Check Matrix Symmetry
#'
#' @param mat A numeric matrix.
#'
#' @return Logical. \code{TRUE} if the matrix is symmetric.
#'
#' @keywords internal
is_symmetric <- function(mat) {
  all.equal(mat, t(mat)) == TRUE
}

#' Estimate Covariate-Adjusted Feature Matrix (X Tilde)
#'
#' Estimates a reduced-dimension feature matrix that maximizes association with
#' covariates A while penalizing association with outcome Y. This is used to
#' obtain X tilde in the scaPCoA framework.
#'
#' @param X.orig A numeric matrix of original features.
#' @param A A numeric matrix of covariates.
#' @param Y A numeric vector or matrix of outcome values.
#' @param npc Integer. Number of principal components to compute.
#' @param lambda_ls A numeric vector of candidate lambda values for tuning.
#'
#' @return A list with:
#'   \describe{
#'     \item{X}{The projected feature matrix (X tilde).}
#'     \item{d}{Number of positive eigenvalues.}
#'     \item{lambda}{The optimal lambda selected by tuning.}
#'   }
#'   Returns \code{NULL} on error.
#'
#' @export
estimate_x_tilde <- function(X.orig, A, Y, npc,lambda_ls) {
  tryCatch({
    # Scale the input matrix
    X.orig <- scale(X.orig)

    # Tune lambda
    lambda.tune <- pedecure.tune( #PeDecURe::pedecure.tune
      X.orig = X.orig,
      X.max = X.orig,
      X.penalize = X.orig,
      lambdas = lambda_ls,
      A = Y,
      Y = A,
      nPC = npc,
      centerX=T,
      scaleX=T
    )

    best.lambda <- lambda.tune$lambda_tune
    print(best.lambda)

    # Run pedecure
    pedecure.out <- pedecure( #PeDecURe::pedecure
      X = X.orig,
      X.penalize = X.orig,
      A = Y,
      Y = A,
      lambda = best.lambda,
      nPC = npc,
      centerX=T,
      scaleX=T
    )
    #pedecure.out <-pedecure.out_ls$eig_out

    v_tilde <- pedecure.out$vectors # Eigenvectors
    eigen_values <- pedecure.out$values # Eigenvalues

    epsilon <- sqrt(.Machine$double.eps)*10000
    positive_indices <- which(eigen_values > epsilon) # Identify indices of positive eigenvalues

    d <- length(positive_indices)
    print(paste("number of positive eigenvalues:", d))

    if (d == length(eigen_values)) {
      print("Eigenvalues are all positive.")
    }

    # Project the data onto the eigenvectors
    PC.scores <- X.orig %*% v_tilde

    X_tilde <- PC.scores

    print(paste("Dim:", dim(X_tilde)))

    return(list(X = X_tilde, d = d, lambda = best.lambda))
  }, error = function(e) {
    # Handle errors and retry with reduced npc
    message("Error encountered: ", e$message)
    # if (npc > 1) {
    #  new_npc <- max(3,round(npc/ 4 *3))
    #   message("Retrying with npc = ", new_npc)
    #   return(estimate_x_tilde(X.orig, A, Y, new_npc))
    # } else {
    #   stop("All retries failed. Unable to estimate X_tilde.")
    # }
  })
}


#' Estimate Outcome-Associated Feature Matrix (X Star)
#'
#' Estimates a reduced-dimension feature matrix that maximizes association with
#' outcome Y while penalizing association with covariates A. This is used to
#' obtain X star in the scaPCoA framework.
#'
#' @param X.orig A numeric matrix of original features.
#' @param Y A numeric vector or matrix of outcome values.
#' @param A A numeric matrix of covariates.
#' @param npc Integer. Number of principal components to compute.
#' @param lambda_ls A numeric vector of candidate lambda values for tuning.
#'
#' @return A list with:
#'   \describe{
#'     \item{X}{The projected feature matrix (X star).}
#'     \item{d}{Number of positive eigenvalues.}
#'     \item{lambda}{The optimal lambda selected by tuning.}
#'   }
#'   Returns \code{NULL} on error.
#'
#' @export
estimate_x_star <- function(X.orig, Y, A, npc, lambda_ls) {
  tryCatch({
    # Scale the input matrix
    X.orig <- scale(X.orig)

    # Tune lambda
    lambda.tune <- pedecure.tune(
      X.orig = X.orig,
      X.max = X.orig,
      X.penalize = X.orig,
      lambdas = lambda_ls,
      A = A,
      Y = Y,
      nPC = npc,
      centerX=T,
      scaleX=T
    )

    best.lambda <- lambda.tune$lambda_tune
    print(best.lambda)

    # Run pedecure
    pedecure.out <- pedecure( #PeDecURe::pedecure
      X = X.orig,
      X.penalize = X.orig,
      A = A,
      Y = Y,
      lambda = best.lambda,
      nPC = npc,
      centerX=T,
      scaleX=T
    )
    #pedecure.out = pedecure.out_ls$eig_out


    v_star <- pedecure.out$vectors # Eigenvectors
    eigen_values <- pedecure.out$values # Eigenvalues

    epsilon <- sqrt(.Machine$double.eps)*10000
    positive_indices <- which(eigen_values > epsilon) # Identify indices of positive eigenvalues

    d <- length(positive_indices)
    print(paste("number of positive eigenvalues:", d))

    if (d == length(eigen_values)) {
      print("Eigenvalues are all positive.")
    }

    # Project the data onto the eigenvectors
    PC.scores <- X.orig %*% v_star

    X_star <- PC.scores

    print(paste("Dim:", dim(X_star)))

    return(list(X = X_star, d = d, lambda = best.lambda))
  }, error = function(e) {
    # Handle errors and retry with reduced npc
    message("Error encountered: ", e$message)
    # if (npc > 1) {
    #  new_npc <- max(3,round(npc/ 4 *3))
    #   message("Retrying with npc = ", new_npc)
    #   return(estimate_x_tilde(X.orig, A, Y, new_npc))
    # } else {
    #   stop("All retries failed. Unable to estimate X_tilde.")
    # }
  })
}


#' Supervised and Covariate-Adjusted Principal Coordinates Analysis (scaPCoA)
#'
#' Main function for performing scaPCoA or scaPCoA-R (with residualization).
#' Finds principal coordinates that maximize association with the outcome Y
#' while minimizing association with confounding covariates A.
#'
#' @param obj A list containing:
#'   \describe{
#'     \item{A}{A data frame or matrix of covariates.}
#'     \item{Y}{A numeric vector of outcome values.}
#'     \item{X_tilde}{A numeric matrix of features (e.g., microbiome counts).}
#'   }
#' @param lambda_ls A numeric vector of candidate lambda values for tuning.
#' @param pcoa Logical. If \code{TRUE}, use Bray-Curtis PCoA decomposition.
#'   If \code{FALSE}, use centered log-ratio (CLR) transformation.
#' @param nPC Integer. Number of principal components to compute. Default is 3.
#' @param residualization Logical. If \code{TRUE}, use the scaPCoA-R variant
#'   with residualized matrices. Default is \code{FALSE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{PC.scores}{A matrix of principal coordinate scores.}
#'     \item{v_hat}{A matrix of estimated eigenvectors.}
#'     \item{best.lambda}{The optimal lambda selected by tuning.}
#'     \item{eigen_values}{A numeric vector of eigenvalues.}
#'     \item{scores.partial.cor}{A matrix of partial correlations between
#'       PC scores and Y/A.}
#'   }
#' @importFrom compositions clr
#' @export
scapcoa <- function(obj, lambda_ls, pcoa, nPC=3, residualization = FALSE){
    results = list()

    A = obj$A
    Y = obj$Y
    X = obj$X_tilde

    Y <- as.matrix(Y)
    A <- as.matrix(A)
    X <- as.matrix(X)

    n = nrow(X)

    if(length(unique(Y)) > 2){
      y_cont = TRUE
    }else{
      y_cont = FALSE
    }

    if(pcoa == FALSE){
      print("s+acPCA (Residualization)")
      X_hat <- clr(X) #centerd log transformation
      X_hat <- X_hat[, apply(X_hat, 2, var) != 0]
      dim(X_hat)
    }else{
      print("s+acPCoA (Residualization)")
      print("Use bray curtis for X decomposition. Obtained X_hat.")
      eigen_output <-  eigen_bray(X, pos=TRUE)
      eigen_vectors <- eigen_output$eigenvectors
      eigen_values = eigen_output$eigenvalues

      X_hat <- eigen_vectors %*% diag(sqrt(eigen_values))

      dim(X_hat)

      X_hat_train = X_hat

    }

    X.orig = X_hat

    if(residualization == TRUE){
      print("Run scaPCoA-R")
      resid.dat = get.resid(X.orig,Y, A)
      X.max = resid.dat$X.star
      X.penalize = resid.dat$X.tilde
    }else{
      print("Run scaPCoA")
      X.max = X_hat
      X.penalize =X_hat
    }

    lambda.tune = pedecure.tune(X.orig = X.orig,
                                  X.max = X.max,
                                  X.penalize = X.penalize,
                                  lambdas = lambda_ls,
                                  A = A,
                                  Y = Y,
                                  nPC = nPC, centerX = T, scaleX = T)

      best.lambda = lambda.tune$lambda_tune
      print(paste('Best Lambda from tuning:',best.lambda))

      pedecure.out = pedecure(X = X.max,
                              X.penalize = X.penalize,
                              A = A,
                              Y = Y,
                              lambda = best.lambda,
                              nPC = nPC, centerX = T, scaleX = T)

      v_hat = pedecure.out$vectors
      PC.scores = X.orig%*%v_hat
      eigen_values <- pedecure.out$values

      # Partial Correlation
      # Look at correlations between the first few PC scores and the nuisance variables (A1, A2, Y)
      cor.scores = partial.cor(PC.scores, A, Y)

      scores.partial.cor = cor.scores$partial$estimates

      l = dim(A)[2]

      # Loop for y and PCs
      for (j in 1:3) {  # Assuming 3 principal components (PC1, PC2, PC3)
        results[[paste0("pcor_y_pc", j)]] <- abs(scores.partial.cor[l + 1, j])
      }

      # Loop for all A (a1 to a10) and PCs
      for (i in 1:l) {  # Loop through A1, A2, ..., Al
        for (j in 1:3) {  # Loop through PC1, PC2, PC3
          results[[paste0("pcor_a", i, "_pc", j)]] <- abs(scores.partial.cor[i, j])
        }
      }

      return(list(PC.scores = PC.scores, v_hat = v_hat, best.lambda = best.lambda, eigen_values = eigen_values, scores.partial.cor= scores.partial.cor))
  }



