
#' Information Criterion
#'
#' @param model Regression model.
#' @param x input matrix, rows represent observations and columns represent
#' features. column names and row names are recommended.
#' @param y response variable, restricted to two classes.
#' @param family GLM family.
#' @param lambda Used in cv.glmnet.
#' @param metric AIC,  AICc, BIC or EBIC.
#' @param gamma value for EBIC.
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#'
#' @return The value of the metric when calculated from the model.
#'
#' @examples
#' 
#' library('ropls')
#' data("sacurine") #Load sacurine dataset from the 'ropls' package
#' 
#' dat = sacurine$dataMatrix
#' outcomes = sacurine$sampleMetadata$gender
#' 
#' inf_criterion(model = rep(TRUE, NCOL(dat)), x = dat, y = outcomes)

inf_criterion = function(model,
                         x,
                         y,
                         family = "binomial",
                         lambda = 'lambda.1se',
                         metric = "BIC",
                         gamma = gamma) {
  n = base::length(y)
  P = base::NCOL(x)

  xNew = x[, model]

  CVGlmFit = glmnet::cv.glmnet(
    x = xNew,
    y = y,
    standardize = TRUE,
    alpha = 0,
    family = family
  )

  GlmFit = glmnet::glmnet(
    x = xNew,
    y = y,
    standardize = TRUE,
    alpha = 0,
    family = family,
    lambda = CVGlmFit[[lambda]]
  )

  XScaled = scale(x)
  ld = CVGlmFit[[lambda]] * diag(ncol(XScaled))
  H = XScaled %*% solve(t(XScaled) %*% XScaled + ld) %*% t(XScaled)
  df = fBasics::tr(H)
  
  deviance = (1L - GlmFit$"dev.ratio") * GlmFit$"nulldev"

  base::switch(
    metric,
    "AIC" = deviance + 2 * df,
    "AICC" = deviance + 2 * df + 2 * (df ^ 2 + df) / (n - df - 1),
    "BIC" = deviance + log(n) * df,
    "EBIC" = deviance + log(n) * df + 2 * gamma * lchoose(n = P, k = df),
    NA
  )

}
