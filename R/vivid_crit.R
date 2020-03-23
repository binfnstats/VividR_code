
#' Change the information criteria used in VIVID
#'
#' @param vividObj An object passed from the vivid() funciton.
#' @param x input matrix, rows represent observations and columns represent
#' features. column names and row names are recommended.
#' @param y response variable, restricted to two classes.
#' @param metric AIC,  AICc, BIC or EBIC.
#' @param lambda Used in cv.glmnet
#'
#' @return A list of results:
#' \itemize{
#' \item{coefficients}: A numeric matrix of estimated coefficient values over re-samples.
#' \item{varMat}: A numeric matrix of the estimated variance matrix explain in the VIVID method.
#' \item{varClust}: The hclust() output of the clustering performed in the VIVID method.
#' \item{selection}: A boolean matrix where each row indicates a model identified as a potential final model.
#' \item{sizes}: A integer vector of model sizes in the seleciton path.
#' \item{compareMethod}: A character value indicating the Information Criterion used for the results.
#' \item{compareValues}: A numeric vector of Infomation Criterion values for each model in the selection path.
#' \item{optModel}: A boolean vector indicating whether or not each feature is include in the final model.
#' \item{optFeatures}: A character vector containing all features in the final model.
#' \item{vividSplit}: A boolean value indicating whether or not this output is from the vivid_split() funciton.
#' }
#' @export
#'
#' @examples 
#' library('ropls')
#' data("sacurine") #Load sacurine dataset from the 'ropls' package
#'
#' dat = sacurine$dataMatrix
#' outcomes = sacurine$sampleMetadata$gender
#' 
#' vividResults = vivid(x = dat,
#'       y = outcomes)
#'       
#' newVivid = vivid_crit(vividObj = vividResults, x = dat, y = outcomes, metric = 'AIC')
#' 
#' newVivid$optFeatures

vivid_crit = function(vividObj, x, y, metric, lambda = 'lambda.1se') {
  selectionPath = vividObj$selection

  compareValues = base::apply(
    X = selectionPath,
    MARGIN = 1,
    FUN = inf_criterion,
    x = x,
    y = y,
    lambda = lambda,
    metric = metric
  )
  optModel = selectionPath[which.min(compareValues), ]

  optFeatures = base::names(x = optModel)[base::unlist(x = optModel)]

  vividObj$compareMethod = metric
  vividObj$compareValues = compareValues
  vividObj$optModel = optModel
  vividObj$optFeatures = optFeatures

  return(vividObj)
}
