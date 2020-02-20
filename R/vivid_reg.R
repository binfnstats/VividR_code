

#' Regression method used in the VIVID method
#'
#' @param weight A p length vector
#' @param x input matrix, rows represent observations and columns represent
#' features. column names and row names are recommended.
#' @param y response variable, restricted to two classes.
#' @param crossfold Number of crossfold samples.
#' @param lambda Used in cv.glmnet.
#'
#' @return
#' @export
#'
#' @examples 1

vivid_reg = function(weight, x, y, crossfold = 10, lambda = "lambda.min") {
  # Fit a ridge regression with observation weights

  nFolds = crossfold
  foldId = base::sample(base::rep(x = seq(nFolds),
                                  length.out = nrow(x)))
  ridgeCV = tryCatch(glmnet::cv.glmnet(x = x,
                                      y = y,
                                      standardize = TRUE,
                                      alpha = 0,
                                      family = "binomial",
                                      weights = weight),
      error = function(weight, x, y){
        base::print(
          "Error in predmat[which, seq(nlami)] <- preds : replacement has length zero. Fixed lambda used."
        )
        ridgeCVtest = glmnet::cv.glmnet(
          x = x,
          y = y,
          standardize = TRUE,
          alpha = 0,
          family = "binomial",
          weights = weight,
          lambda = base::exp(base::seq(
            from = log(0.001),
            to = log(50),
            length.out = 100)))
      return(ridgeCVtest)
      }
  )


  # Compute coefficients
  ridgeCoef = glmnet::coef.glmnet(object = ridgeCV$glmnet.fit,
                                  s = ridgeCV[[lambda]])

  estCoef = ridgeCoef[-1,]

  return(estCoef)
}
