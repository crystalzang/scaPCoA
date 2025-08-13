
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

supPCoA <- function(obj, lambda_ls, pcoa){
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
      X_hat <- clr(X_train) #centerd log transformation
      X_hat <- X_hat[, apply(X_hat, 2, var) != 0]
      dim(X_hat)
    }else{
      print("s+acPCoA (Residualization)")
      print("Use bray curtis for X decomposition. Obtained X_hat.")
      eigen_output <-  eigen_bray(X_train, pos=TRUE)
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
                                  nPC = 3, centerX = T, scaleX = T)
      
      
      best.lambda = lambda.tune$lambda_tune
      print(paste('Best Lambda from tuning:',best.lambda))
      
      lambda_optimal_ls[["final"]] = best.lambda
      
      # run pedecure:
      pedecure.out = pedecure(X = X.max,
                              X.penalize = X.penalize,
                              A = A,
                              Y = Y,
                              lambda = best.lambda,
                              nPC = 3, centerX = T, scaleX = T)
      
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
   
      pc_ls <- pc_train_function(PC.scores_train, Y_train, pcoa=pcoa, alpha=0.8, title = title)
      
      pc_plot <- pc_ls$pc_plot
      pc_df <- pc_ls$pc_df 
      
      return(list(pc_plot = pc_plot, pc_df = pc_df, results = results, PC.scores = PC.scores,   v_hat = v_hat, lambda_optimal_ls=lambda_optimal_ls, U_df = U_df,eigen_values = eigen_values_out ))
    
  }
 
