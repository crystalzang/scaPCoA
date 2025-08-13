
#zero-inflated beta regression
library(gamlss)
library(gamlss.mx)
library(RSpectra ) #eigs_sym
library(MASS)
library(polycor)
library(mltools) #mcc
library(ltm) #point biserial correlation
library(rstatix) # correlation (ANOVA-based Eta Squared η²)


color_palette_cz <- c("#27548A", "#DDA853", "#3F7D58", 
                      "#C84D6A","#A6C36F", "#3B8D5D")   # teal)



###########################################
######  From Pedecure  ####################
###########################################

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


plsMax = function(v,args){
  X = args$X
  Y = args$Y
  
  out = tcrossprod(crossprod(X,Y))%*%matrix(v,ncol=1)
  # Replace NA values with 0
  out[is.na(out)] <- 0
  
  return(out)
}


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


check_correlation <- function(cor) {
  if (cor < -1 | cor > 1) {
    stop("Error: Correlation is out of bounds (-1 to 1). Check your data.")
  }
}


library(rstatix)

# scores = PC.scores
# A = A_train
# Y = Y_train
partial.cor = function(scores, A, Y){
  if(is.character(Y)){
    print("Change Y from character to numeric")
  }
  # input: scores:matrix of PC scores (e.g., Xv1, Xv2, Xv3)
  # estimate partial (+ marginal) correlations between each PC score and Y|A and A|Y
  # output: list of partial and marginal correlations between each PC score and A and Y
  scores = as.matrix(scores)
  A = as.matrix(A)
  colnames(scores)  = paste0("PC",1:ncol(scores))
  n = nrow(scores)
  df.partial = n-2-ncol(A) # n - 2 - number of variables conditioning on
  
  pcor.Xv.j.Y.given.A = vector(mode = "numeric",length = ncol(scores))
  pcor.Xv.j.A.given.Y = rep(list(vector(mode = "numeric", length = ncol(scores))), ncol(A))
  names(pcor.Xv.j.A.given.Y) = colnames(A)
  
  t.statistics = matrix(nrow=(ncol(A) + 1),ncol=ncol(scores))
  rownames(t.statistics) = c(paste0(colnames(A)),"Y")
  colnames(t.statistics) = paste0("PC",1:ncol(scores))
  
  for (j in 1:ncol(scores)) {
    
    omega.j <- matrix(0, ncol(A) + 2, ncol(A) + 2)  # Initialize omega.j
    colnames(omega.j) <- rownames(omega.j) <- c(colnames(A), "Y", paste0("Xv", j))
    
    # A and Y
    for (a in 1:ncol(A)) {
      # Check if A[, a] has zero standard deviation
      if (is.numeric(A[, a]) && sd(A[, a]) == 0) {
        stop(paste("Error: Column", a, "in A has no standard deviation (all values are the same)."))
      }
      
      is_Y_binary <- length(unique(Y)) == 2
      is_A_binary <- length(unique(A[, a])) == 2
      is_A_categorical <- is.factor(A[, a]) || is.character(A[, a]) || (!is_A_binary && length(unique(A[, a])) <= 10)
      
      
      if (is_A_binary && is_Y_binary) {
        #print("Binary Ai and Binary Y: MCC")
        Y_vector = as.vector(Y)
        Y_vector <- as.numeric(as.character(Y_vector))
        
        cor_i = mcc(A[,a], Y_vector)
      }else if(is_A_binary && !is_Y_binary){
        #print("A is Binary, Y is Continuous: Point-Biserial Correlation")
        cor_i = biserial.cor(Y, A[, a],use = "all.obs", level = 2)

      }else if (!is_A_binary && is_Y_binary) {
        if (is_A_categorical) {
          # Categorical A, Binary Y: Eta Squared
          #print("Categorical A, Binary Y: Eta Squared")
          df_tmp <- data.frame(score = as.vector(Y), group = as.factor(A[, a]))
    
          test <- anova_test(data = df_tmp, dv = score, between = group)
          cor_i <-  test$ges[1]
        } else {
          #print("Continuous A, Binary Y: Point-Biserial")
          cor_i <- biserial.cor(A[, a], Y, use = "all.obs", level = 2)
        }
        
      } else if (is_A_categorical && !is_Y_binary) {
        # Categorical A, Continuous Y: Eta Squared
        df_tmp <- data.frame(score = as.vector(Y), group = as.factor(A[, a]))
        test <- anova_test(data = df_tmp, dv = score, between = group)
        
        cor_i <- test$ges[1]
        
      }else {
        # Both Ai and Y are continuous: Use Pearson Correlation
        #print("Continuous Ai and Continuous Y: Pearson Correlation")
        cor_i = cor(A[, a], Y, method = "pearson")
      }
      check_correlation(cor_i)
      omega.j[a, ncol(A) + 1] <- omega.j[ncol(A) + 1, a] <- cor_i
    }
    
    
    # Compute correlation between scores and Y
    if (length(table(Y)) == 2) {
      #print("Y is Binary: Point-Biserial Correlation with scores")
      cor_i = biserial.cor( scores[, j],Y, use = "all.obs", level = 2)
      check_correlation(cor_i)
      omega.j[ncol(A) + 1, ncol(A) + 2] <- omega.j[ncol(A) + 2, ncol(A) + 1] <- cor_i
    } else {
      #print("Y is Continuous: Pearson Correlation with scores")
      cor_i = cor(Y, scores[, j], method = "pearson")
      check_correlation(cor_i)
      omega.j[ncol(A) + 1, ncol(A) + 2] <- omega.j[ncol(A) + 2, ncol(A) + 1] <- cor_i
    }
    
    # Compute correlation between scores and A
    for (a in 1:ncol(A)) {
      if (length(table(A[, a])) == 2) {
        #print("Ai is Binary: Point-Biserial Correlation with scores")
        cor_i =biserial.cor(scores[, j], A[, a], use = "all.obs", level = 2)
        check_correlation(cor_i)
        omega.j[a, ncol(A) + 2] <- omega.j[ncol(A) + 2, a] <- cor_i
      } else if(length(table(A[, a])) > 2 && length(table(A[, a])) < 10 ){
        # Ai is Categorical with >2 levels: use eta squared (0 to 1)
        print("A is categorical with > 2 levels. Use ANOVA eta score.")
        df_tmp <- data.frame(score = scores[, j], group = A[, a])
        
        result <- tryCatch({
          test <- anova_test(data = df_tmp, dv = score, between = group)
          eta_sq <- test$ges[1]  # Generalized eta squared
          check_correlation(eta_sq)
          eta_sq
        }, error = function(e) {
          warning(paste("ANOVA failed for A column", a))
          NA
        })
        
        omega.j[a, ncol(A) + 2] <- omega.j[ncol(A) + 2, a] <- result
        
      }else {
        #print("Ai is Continuous: Pearson Correlation with scores")
        cor_i = cor(A[, a], scores[, j], method = "pearson")
        check_correlation(cor_i)
        
        omega.j[a, ncol(A) + 2] <- omega.j[ncol(A) + 2, a] <- cor_i
      }
    }
    
    # Fill in the correlations among the A variables
    for (a in 1:ncol(A)) {
      for (b in a:ncol(A)) {  # We only fill in the upper triangle (or lower triangle since it's symmetric)
        if (a != b) {  # Skip diagonal as it's already filled with variance
          if (length(table(A[, a])) == 2 && length(table(A[, b])) == 2) {
            #print("Both Ai and Ab are Binary: Matthews Correlation Coefficient (MCC)")
            omega.j[a, b] <- omega.j[b, a] <- mcc(A[, a], A[, b])
          } else if (length(table(A[, a])) == 2 && length(table(A[, b])) > 2) {
            #print("Ai is Binary, Ab is Continuous: Point-Biserial Correlation")
            cor_i = biserial.cor(A[, b], A[, a], use = "all.obs", level = 2)
            check_correlation(cor_i)
            omega.j[a, b] <- omega.j[b, a] <- cor_i
          } else if (length(table(A[, a])) > 2 && length(table(A[, b])) == 2) {
            #print("Ai is Continuous, Ab is Binary: Point-Biserial Correlation")
            cor_i = biserial.cor(A[, a], A[, b], use = "all.obs", level = 2)
            check_correlation(cor_i)
            omega.j[a, b] <- omega.j[b, a] <- cor_i
          } else {
            #print("Both Ai and Ab are Continuous: Pearson Correlation")
            cor_i = cor(A[, a], A[, b], method = "pearson")
            check_correlation(cor_i)
            omega.j[a, b] <- omega.j[b, a] <- cor_i
          }
        }
      }
    }
    
    # Diagonal elements (variances)
    # omega.j[ncol(A) + 1, ncol(A) + 1] <- var(Y)
    # omega.j[ncol(A) + 2, ncol(A) + 2] <- var(scores[, j])
    # for (a in 1:ncol(A)) {
    #   omega.j[a, a] <- var(A[, a])
    # }
    diag(omega.j) = 1
    
    diag(omega.j) <- diag(omega.j) + 1e-10
    
    # Calculate the inverse of omega.j. (If not singular, a pseudoinverse is calculated)
    det_omega <- det(omega.j)
    if (abs(det_omega) < 1e-10) {
      print(paste("Determinant =", det_omega, "Matrix is nearly singular, performing Moore-Penrose Pseudoinverse"))
      P.j <- ginv(omega.j)
    } else {
      P.j <- solve(omega.j)
    }
    
    # if omega.j is already correlation matrix, skip this step:
    #pcor.estimate.j <- -cov2cor(P.j)  #convert var-covariance matrix to correlation matrix
    
    #pcor.Xv.j.Y.given.A[j] <- pcor.estimate.j["Y", (ncol(A) + 2)]
    
    # pcor.Xv.j.Y.given.A[j] <- P.j["Y", (ncol(A) + 2)]
    # 
    # 
    # for (a in 1:ncol(A)) {
    #  # pcor.Xv.j.A.given.Y[[colnames(A)[a]]][j] <- pcor.estimate.j[colnames(A)[a], (ncol(A) + 2)]
    #    pcor.Xv.j.A.given.Y[[colnames(A)[a]]][j] <- P.j[colnames(A)[a], (ncol(A) + 2)]
    #   
    # }
    
    # Compute partial correlation between scores and Y given A
    pcor.Xv.j.Y.given.A[j] <- -P.j["Y", ncol(A) + 2] / sqrt(P.j["Y", "Y"] * P.j[ncol(A) + 2, ncol(A) + 2])
    
    # Compute partial correlations between scores and each A given Y
    for (a in 1:ncol(A)) {
      pcor.Xv.j.A.given.Y[[colnames(A)[a]]][j] <- -P.j[colnames(A)[a], ncol(A) + 2] / sqrt(P.j[colnames(A)[a], colnames(A)[a]] * P.j[ncol(A) + 2, ncol(A) + 2])
    }
  }
  
  partial.correlations = rbind(do.call("rbind",pcor.Xv.j.A.given.Y),
                               Y = pcor.Xv.j.Y.given.A)
  
  colnames(partial.correlations) = paste0("PC",1:ncol(scores))
  
  out.partial = list(estimates = partial.correlations,
                     statistics = t.statistics,
                     df = df.partial)

  return(list(partial = out.partial))
  
}


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

is_symmetric <- function(mat) {
  all.equal(mat, t(mat)) == TRUE
}


sparsity = function(X){
  length(which(X==0))/prod(dim(X))
}


#sparsity(X)

expit <- function(x) 1 / (1 + exp(-x))


#Nyström Extension
nystrom <- function(A_tt, eigenvectors, eigenvalues) {
  # A_tt: n_train x n_test Gower-centered similarity matrix
  # eigenvectors: Q_M, n_train x M matrix from eigen_bray
  # eigenvalues: Lambda_M, vector of M positive eigenvalues
  
  # Defensive: Make sure eigenvalues is a column vector if needed
  Lambda_sqrt_inv <- diag(1 / sqrt(eigenvalues))  # M x M matrix
  
  # Nyström projection formula: A_tt^T Q_M Λ^{-1/2}
  X_hat_test <- t(A_tt) %*% eigenvectors %*% Lambda_sqrt_inv  # n_test x M
  
  return(X_hat_test)
}

trosset_out_of_sample <- function(X_train_hat, A_tt, return_beta = FALSE) {
  n_test <- ncol(A_tt)
  M <- ncol(X_train_hat)
  n_train <- nrow(X_train_hat)
  X_test_hat <- matrix(0, nrow = n_test, ncol = M)
  
  for (j in 1:n_test) {
    b_j <- A_tt[, j]
    
    if (!return_beta) {
      X_test_hat[j, ] <- MASS::ginv(X_train_hat) %*% b_j
    } else {
      beta_j <- mean(b_j)  # approximate squared norm
      
      objective <- function(y) {
        fit <- X_train_hat %*% y
        err1 <- sum((b_j - fit)^2)
        err2 <- (beta_j - sum(y^2))^2
        return(err1 + err2)
      }
      
      y0 <- MASS::ginv(X_train_hat) %*% b_j
      opt <- optim(par = y0, fn = objective, method = "BFGS")
      X_test_hat[j, ] <- opt$par
    }
  }
  
  return(X_test_hat)
}


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

# X.orig = X_hat_d0
# npc = d

estimate_x_star <- function(X.orig, Y, A, npc, lambda_ls) {
  tryCatch({
    # Scale the input matrix
    X.orig <- scale(X.orig)
    
    # Tune lambda
    lambda.tune <- pedecure.tune( #PeDecURe::pedecure.tune
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
supPCoA <- function(obj, lambda_ls, pcoa, nPC=3){
    results = list()
    lambda_optimal_ls <- list()
    
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
      #X_hat <- X_hat[, apply(X_hat, 2, var) != 0]
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
    
    resid.dat = get.resid(X.orig,Y, A)
    
    X.max = resid.dat$X.star
    X.penalize = resid.dat$X.tilde
      
    lambda.tune = pedecure.tune(X.orig = X.orig,
                                  X.max = X.max,
                                  X.penalize = X.penalize,
                                  lambdas = lambda_ls,
                                  A = A,
                                  Y = Y,
                                  nPC = nPC, centerX = T, scaleX = T)
      
      
      best.lambda = lambda.tune$lambda_tune
      print(paste('Best Lambda from tuning:',best.lambda))
      
      lambda_optimal_ls[["final"]] = best.lambda
      
      # run pedecure:
      pedecure.out = pedecure(X = X.max,
                              X.penalize = X.penalize,
                              A = A,
                              Y = Y,
                              lambda = best.lambda,
                              nPC = nPC, centerX = T, scaleX = T)
      
      v_hat = pedecure.out$vectors
      PC.scores = X.orig%*%v_hat
      eigen_values <- pedecure.out$values
      
      U_df = calc_feature_importance(PC.scores, eigen_values, X)
      eigen_values_out = eigen_values
      

    
      #### 4. Partial Correlation  #### 
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
      
   
      
      ##### PC Plot    ##### 
      if(pcoa ==TRUE){
        title = "s+acPCoA (Residualization)"
      }else{
        title = "s+acPCA (Residualization)"
      }
   
      pc_ls <- pc_train_function(PC.scores, Y, pcoa=pcoa, alpha=0.8, title = title)
      
      pc_plot <- pc_ls$pc_plot
      pc_df <- pc_ls$pc_df 
      
      return(list(pc_plot = pc_plot, pc_df = pc_df, results = results, PC.scores = PC.scores,   v_hat = v_hat, lambda_optimal_ls=lambda_optimal_ls, U_df = U_df,eigen_values = eigen_values_out ))
    
  }
  



generate_data <- function(seed, n,p,l, cor_val, y_type){
  set.seed(seed)
  mu <- rep(0, l)
  
  #### Generate A ####
  # Generate the block-wise correlation matrix
  # Construct correlation matrix
  sigma <- matrix(cor_val, nrow = l, ncol = l) + diag(1 - cor_val, l)
  
  # Force positive semi-definiteness
  sigma <- as.matrix(nearPD(sigma, corr = TRUE)$mat)
  
  
  A <- mvrnorm(n = n, mu = mu, Sigma = sigma)   # Generate Ai,1 and Ai,2
  A[, 2] <- rbinom(n, 1, 0.1) # this works
  
  set.seed(seed)
  n = as.numeric(n)
  N.d = zinLDA::rdu(n = n, min = 1000, max = 10000)
  sim1 = zinLDA::simulateZINLDA(D = n, #number of sample
                                V = p,  #number of feature
                                N = N.d, #vector of length D containing the total number of sequencing readings per sample.
                                K = 5,  #num_cluster
                                Alpha = 0.1, #scalar symmetric hyperparameter of the Dirichlet prior on theta.
                                Pi = 0.7,# Probability of inflation (excess zeros)
                                a = 0.6, # Increase to allow structured feature correlations
                                b = 10#lower more sparse (b=2:prop = 0.89, b=0.2:p=0.965, b=7:prop=0.76)
  )
  
  data <-  as.data.frame(sim1$sampleTaxaMatrix)
  
  sum(colSums(data)==0)/p
  
  #prop_feature_zero = apply(data, 2, function(col) sum(col==0)/length(col))
  #hist(prop_feature_zero)
  
  #prop_sample_zero = apply(data, 1, function(row) sum(row==0)/length(row))
  #hist(prop_sample_zero)
  
  data_c <- data[, colSums(data != 0) > 0]
  data_c <- data_c[, apply(data_c, 2, var) != 0]
  data_c <- data_c[rowSums(data_c != 0) > 0, ]
  #dim(data_c)
  print(paste("Exclude features with all zeros. Exclude features with variance of zero. 
              Number of features = ",dim(data_c)[2]))
  
  
  # normalized
  data_rel = apply(data_c, 2, function(x) x / sum(x))
  X <- as.matrix(data_rel)
  
  hist(cor(X, A[,1]))
  # Add some correlation between X and A.
  num_feature = ncol(X)
  num_sample = nrow(X)
  # 1. Construct features correlated with A1
  n_selected <- ceiling(0.6 * num_feature) #Randomly select 30% of the columns in X
  set.seed(seed)
  selected_cols <- sample(1:num_feature, n_selected)
  # Construct a signal matrix from the selected X columns
  error1 <- rnorm(num_sample, mean = 0, sd = 0.4)
  error1[error1 < 0] <- 0
  for (j in selected_cols) {
    X[, j] <- X[, j] + 0.2 * (A[, 1]-min(A[,1])) + error1
  }
  
  # 2. Construct features correlated with A2 and A3
  set.seed(seed+2)
  selected_cols <- sample(1:num_feature, n_selected)
  # Construct a signal matrix from the selected X columns
  error2 <- rnorm(num_sample, mean = 0, sd = 0.4)
  error2[error2 < 0] <- 0
  for (j in selected_cols) {
    X[, j] <- X[, j] + 0.3*A[,2] + 0.1 * (A[, 3]-min(A[,3])) + error2
  }
  
  
  X = apply(X, 2, function(x) x / sum(x))
  
  cor(X, A[,1])
  hist(cor(X, A[,1]))
  
  # Compute the Pearson correlation matrix
  cor_matrix <- cor(X, method = "pearson")
  heatmap(cor_matrix)
  hist(cor_matrix)
  
  A <- as.data.frame(A)
  colnames(A)  <- paste0("A", 1:l)
  # correlation with covariates
  correlation_feature_A <- cor(X, as.matrix(A))
  hist(correlation_feature_A[,1])
  
  
  # # Step 1: Generate W
  ## Option 1. Select the 10 features with the highest means
  feature_means <- colMeans(X)
  top_features <- names(sort(feature_means, decreasing = TRUE)[1:15])
  
  #### Generate W ####
  W <- X[, top_features] # option 1. (This works)
  print("Top 15 features")

  
  W_scaled = sweep(W, 1, rowSums(W)+1e-20, "/")
  #print("Transform W_scaled as relative abundance.")
  W_scaled <- log(W_scaled + 1e-8)
  print("Transform W_scaled as log(relative abundance).")
  
  #W_scaled <- clr(W) #centerd log transformation, this leads to a decrease in partial correlation with Y in testing
  W_scaled <- W_scaled[, apply(W_scaled, 2, var) != 0]
  
  W_scaled <- as.matrix(W_scaled)
  C_A <- as.matrix(rep(1, ncol(W_scaled)))
  
  
  print(paste("Formula for non-linear relationship between A and W:", formula_aw))
  
  if(formula_aw =="0"){
    print("No modification of A.")
  }
  
  #print("Scale numeric covairates A.")
  for(j in 1:ncol(A)){
    if(length(table(A[,j])) > 2){
      A[,j] <- scale(A[,j], center = TRUE, scale = TRUE)
    }
  }
  #hist(A[,1])
  
  
  # Step 2: Generate Z2_i ~ N(0,1) (n x 1 vector)
  #Z2 <- rnorm(n, mean = 0, sd = 1)
  
  # Step 3: Compute subject-level probabilities π
  expit <- function(x) 1 / (1 + exp(-x))
  
  #### Generate Y ####
  print(paste("Formula to generate Y:", formula_y))
  
  
  linear_predictor <- W_scaled %*% C_A + 0.3*A$A1 - 0.01*A$A2 + 0.1*A$A1 * A$A3
  lp_median <- -median(linear_predictor)
  linear_predictor <- linear_predictor + lp_median


  # Step 4: Generate pi ~ Bernoulli(π), and then generate Y from pi.
  if(y_type == "bin"){
    print("Binary Y")
    pi <- 1 / (1 + exp(-linear_predictor)) # Binary Y
    Y <- rbinom(n, size = 1, prob = pi)  # Draw realizations of binary outcomes
    print(table(Y))
  }else if(y_type == "cont"){
    print("Continuous Y")
    Y = scale(linear_predictor) # continuous Y
  }

  top_features <- sub("Feature", "V", top_features)
  
  
  output = list(
    A = A[,1:l],
    Y = Y,
    W = W,
    X_tilde = X,
    top_features = top_features
  )
  return(output)
}


pc_train_function <- function(PC.scores_train, Y_train, pcoa, alpha = 0.8, title = "PC Plot"){
  train_pc_long <- PC.scores_train%>%
    as.data.frame()%>%
    dplyr::mutate(Y = Y_train)%>%
    dplyr::mutate(train_test = "train")
  
  if (pcoa == TRUE){
    colnames(train_pc_long) <- c("PC1", "PC2", "PC3", "Y", "train_test")
  }
  
  pc_all <- train_pc_long
  
  if (y_type == "bin") {
    pc_all$Y <- factor(pc_all$Y, levels = c(0, 1), labels = c("control", "case"))
  }
  
  pc_plot <- ggplot(pc_all, aes(x = PC1, y = PC2, color = Y)) + #, shape = train_test
    geom_point(size = 3.5, alpha = alpha) +  # Increase point size, adjust transparency
    # stat_ellipse(aes(group = train_test), linetype = "dashed", alpha = 0.4) +  # Add ellipses for better grouping
    #scale_shape_manual(values = c(16, 17)) +  # Define specific shapes
    labs(
      x = "Principal Coordinate 1",
      y = "Principal Coordinate 2",
      color = "Outcome (Y)",
      title = title
    ) +
    theme_minimal(base_size = 14) +  # Slightly larger base font size
    theme(
      legend.position = "right",
      legend.key = element_rect(fill = "white"),
      strip.text = element_text(size = 18, face = "bold"), # facet wrap train/test title size
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5)  # <-- Larger & centered
    )
  if(y_type == "bin"){
    pc_plot <- pc_plot +
      scale_color_manual(
        values = c("control" = color_palette_cz[5], "case" = color_palette_cz[4])
      )
    
  }else{
    pc_plot = pc_plot+scale_color_gradient2(
      low = color_palette_cz[4],       #pink
      mid = "#DCE6DA",         # white at 0
      high = color_palette_cz[6],      #green
      midpoint = 0           # center of gradient
    )
  }
  return(list(pc_df = pc_all, pc_plot = pc_plot))
}

