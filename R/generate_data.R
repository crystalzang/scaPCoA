#' Generate Simulated Microbiome Data with Confounders
#'
#' Generates simulated microbiome count data using zero-inflated LDA,
#' along with correlated covariates and binary or continuous outcomes
#' for testing scaPCoA methods.
#'
#' @param seed Integer. Random seed for reproducibility.
#' @param n Integer. Number of samples. Default is 100.
#' @param p Integer. Number of features (taxa). Default is 200.
#' @param l Integer. Number of covariates. Default is 3.
#' @param cor_val Numeric. Correlation between continuous covariates A1 and A2.
#'   Default is 0.4.
#' @param y_type Character. Type of outcome variable. Either \code{"bin"} for
#'   binary or \code{"cont"} for continuous. Default is \code{"bin"}.
#' @param w_theta Numeric. Currently unused. Default is 1.
#' @param w_perc Numeric. Proportion of columns in X to randomly select.
#'   Default is 0.5.
#' @param coef_x Numeric. Coefficient for signal features W in the outcome
#'   model. Default is 0.5.
#' @param coef_a Numeric vector of length 3. Coefficients for A1, A3, and
#'   A1*A2 interaction in the outcome model. Default is \code{c(1, -1, 1)}.
#' @param zinLDA_min Integer. Minimum sequencing depth per sample.
#'   Default is 1000.
#' @param zinLDA_max Integer. Maximum sequencing depth per sample.
#'   Default is 30000.
#'
#' @return A list with:
#'   \describe{
#'     \item{A}{A data frame of covariates (A1, A2, A3).}
#'     \item{Y}{A vector of outcome values (binary or continuous).}
#'     \item{W}{A matrix of signal features selected based on high correlations.}
#'     \item{X_tilde}{The full feature matrix after filtering.}
#'     \item{features_signals}{Column names of the signal features in W.}
#'   }
#'
#' @details
#' The data generation process:
#' \enumerate{
#'   \item Covariates A1 and A2 are drawn from a bivariate normal with
#'     correlation \code{cor_val}. A3 is binary (Bernoulli with p = 0.3).
#'   \item Microbiome counts are simulated using \code{\link[zinLDA]{simulateZINLDA}}
#'     with zero inflation.
#'   \item 15 signal features (W) are selected based on high pairwise
#'     Pearson correlations among relative abundances.
#'   \item The outcome Y is generated from a linear predictor involving W and
#'     the covariates, using a logistic link for binary outcomes.
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats cor rbinom rnorm median
#' @importFrom dplyr mutate arrange group_by ungroup filter pull desc slice_head row_number
#' @importFrom tidyr gather
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- generate_data(seed = 1, n = 100, p = 200, zinLDA_min = 1000, zinLDA_max = 10000)
#' str(dat)
#' }
generate_data <- function(seed = 1, n =100, p = 200, l=3, cor_val=0.4, y_type = "bin", w_theta = 1, w_perc = 0.5, coef_x = 0.5, coef_a=c(1, -1, 1), zinLDA_min = 1000, zinLDA_max = 30000){
  set.seed(seed)
  print(paste("Generate data with seed =", seed))

  #### Generate A ####
  mu <- rep(0, 2)
  sigma <- matrix(cor_val, nrow = 2, ncol = 2) + diag(1 - cor_val, 2)
  A <- mvrnorm(n = n, mu = mu, Sigma = sigma)
  A3 <- rbinom(n, 1, 0.3)

  A <- cbind(A, A3)
  A <- as.data.frame(A)
  colnames(A)  <- paste0("A", 1:l)

  A[,1] <- scale(A[,1], center = TRUE, scale = TRUE)
  A[,2] <- scale(A[,2], center = TRUE, scale = TRUE)


  #### Generate X (feature) ####
  set.seed(seed)
  n = as.numeric(n)
  N.d = zinLDA::rdu(n = n, min = zinLDA_min, max = zinLDA_max) #should try 1,000 - 30,000
  print( paste0('N.d:', zinLDA_min,"-",zinLDA_max))
  print("Simulate ZINLDA")
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

  data_c <- data[, colSums(data != 0) > 0]
  data_c <- data_c[, apply(data_c, 2, var) != 0]
  data_c <- data_c[rowSums(data_c != 0) > 0, ]
  print(paste("Exclude features with all zeros. Exclude features with variance of zero.
              Number of features = ",dim(data_c)[2]))

  X <- as.matrix(data_c)

  x_means <- rowMeans(X)   #sample mean
  num_feature = ncol(X)
  num_sample = nrow(X)

  # Construct features correlated with A1
  print(paste("Randomly select", w_perc*100,"% of the columns in X"))
  n_selected <- ceiling(w_perc * num_feature)
  set.seed(seed)
  selected_cols <- sample(1:num_feature, n_selected)


  X_c <- X[, colSums(X != 0) > 0]
  X_c <- X_c[, apply(X_c, 2, var) != 0]
  print(paste("Exclude features with all zeros. Exclude features with variance of zero.
              Number of features = ",dim(X_c)[2]))

  X <- X_c

  #1. count
  X_count <- round(X)
  #2. relative abundance
  X_rel <- X / rowSums(X)

  cor_matrix <- cor(X_rel, method = "pearson")

  cutoff = 0.5
  # Upper triangle only (avoid duplicates & diagonal)
  cor_upper <- cor_matrix
  cor_upper[lower.tri(cor_upper, diag = TRUE)] <- NA

  cor_pairs <- which(!is.na(cor_upper) & abs(cor_upper) > cutoff, arr.ind = TRUE)
  edge_df <- data.frame(
    feature1 = rownames(cor_matrix)[cor_pairs[,1]],
    feature2 = colnames(cor_matrix)[cor_pairs[,2]],
    correlation = cor_matrix[cor_pairs]
  )

  feature_scores <- edge_df%>%
    mutate(correlation = abs(correlation))%>%
    gather(node, feature, -correlation)%>%
    arrange(desc(correlation)) %>%
    group_by(feature) %>%
    slice_head(n = 1) %>%
    ungroup()%>%
    mutate(selected = row_number() <=15)

  print("W: 15 features with high correlations.")
  high_cor_features <- feature_scores%>%
    filter(selected == TRUE)%>%
    pull(feature)

  ## Option 3. Select features that are in a community with other features and with high correlations
  W <- X[, high_cor_features, drop = FALSE]
  W <- W[, apply(W, 2, var) != 0]
  W <- as.matrix(W)
  C_A <- as.matrix(rep(coef_x, ncol(W)))

  #### Generate Y ####
  linear_predictor <-  W %*% C_A + coef_a[1]*A$A1 + coef_a[2]*A$A3 + coef_a[3]*A$A1 * A$A2

  lp_median <- -median(linear_predictor)
  linear_predictor <- linear_predictor + lp_median

  if(y_type == "bin"){
    pi <- 1 / (1 + exp(-linear_predictor)) # Binary Y
    Y <- rbinom(n, size = 1, prob = pi)  # Draw realizations of binary outcomes
  }else if(y_type == "cont"){
    error_y <- rnorm(num_sample, mean = 0, sd = 0.01) #gaussian error
    Y = scale(linear_predictor) + error_y
  }

  output = list(
    A = A,
    Y = Y,
    W = W,
    X_tilde = X,
    features_signals = colnames(W)
  )
  return(output)
}

