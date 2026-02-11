#' Check That a Correlation Value Is Within Valid Bounds
#'
#' @param cor A numeric correlation value.
#'
#' @return Invisible NULL. Throws an error if \code{cor} is outside the range.
#'
#' @keywords internal
check_correlation <- function(cor) {
  if (cor < -1 | cor > 1) {
    stop("Error: Correlation is out of bounds (-1 to 1). Check your data.")
  }
}


#' Compute Partial Correlations Between PC Scores, Outcome, and Covariates
#'
#' Estimates partial correlations between each principal coordinate score and
#' the outcome Y conditional on covariates A, and between each score and each
#' covariate conditional on Y. Automatically selects the appropriate correlation
#' measure based on variable types (continuous, binary, categorical).
#'
#' @param scores A numeric matrix of principal coordinate scores (rows = samples,
#'   columns = principal coordinates).
#' @param A A numeric matrix of covariates. Column names are used to label
#'   results.
#' @param Y A numeric vector of outcome values (binary or continuous).
#'
#' @return A list containing:
#'   \describe{
#'     \item{partial}{A list with:
#'       \describe{
#'         \item{estimates}{A matrix of partial correlations. Rows correspond to
#'           each covariate and Y. Columns correspond to each PC.}
#'         \item{statistics}{A matrix of t-statistics (currently not populated).}
#'         \item{df}{Degrees of freedom for partial correlation tests
#'           (\code{n - 2 - ncol(A)}).}
#'       }
#'     }
#'   }
#'
#' @details
#' The function constructs a mixed correlation matrix for each PC score using
#' the appropriate measure for each variable pair:
#' \describe{
#'   \item{Continuous vs Continuous}{Pearson correlation}
#'   \item{Binary vs Continuous}{Point-biserial correlation via
#'     \code{\link[ltm]{biserial.cor}}}
#'   \item{Binary vs Binary}{Matthews Correlation Coefficient via
#'     \code{\link[mltools]{mcc}}}
#'   \item{Categorical vs Continuous/Binary}{Eta squared via
#'     \code{\link[rstatix]{anova_test}}}
#' }
#'
#' Partial correlations are derived by inverting the correlation matrix. If the
#' matrix is nearly singular, the Moore-Penrose pseudoinverse is used via
#' \code{\link[MASS]{ginv}}.
#'
#' @importFrom MASS ginv
#' @importFrom mltools mcc
#' @importFrom ltm biserial.cor
#' @importFrom rstatix anova_test
#' @importFrom stats cor sd
#' @export
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
