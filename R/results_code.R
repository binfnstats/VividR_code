library(devtools)
devtools::install_github("binfnstats/VividR")
library('VIVID')

library("sigFeature")
library('caret')

#######################

check_data = function(x, y, B, k, repB, setSeed = 1234567, nCores = parallel::detectCores() - 1, methods, topN = 50, compareMethod = "BIC", data){
  
  p = ncol(x)
  df = base::data.frame(x, y)
  nMethods = length(methods)
  kFolds = list()
  
  output = list()
  
  # K-folds
  set.seed(1234)
  
  for(i in 1:repB){
    kFold = createFolds(y, 
                        k = k, 
                        list = TRUE, 
                        returnTrain = FALSE)
    kFolds[[i]] = kFold
  }
  
  saveRDS(p,paste0(data,"_folds_seed_",setSeed,".rds",collapse=""))
  
  
  for(l in 1:nMethods){
    accVals = matrix(data = 0,
                     nrow = k,
                     ncol = repB)
    aucVals = matrix(data = 0,
                     nrow = k,
                     ncol = repB)
    timeVals = matrix(data = 0,
                      nrow = k,
                      ncol = repB)
    features = list()
    
    for(i in 1:repB){
      kFold = kFolds[[i]]
      
      features2 = list()
      
      for(j in 1:k){
        trainTest = unlist(kFolds[[i]][-j])
        xTrain = x[trainTest,]
        yTrain = y[trainTest]
        xTest = x[-trainTest,]
        yTest = y[-trainTest]
        preProcValues = caret::preProcess(xTrain, method = c("center", "scale"))
        xTrain = stats::predict(preProcValues, xTrain)
        xTest = stats::predict(preProcValues, xTest)
        
        if(methods[l] == "vivid"){
          tic()
          vivid = VIVID::vivid(x = xTrain,
                               y = yTrain,
                               bootstraps = B,
                               cores = nCores,
                               seed = setSeed,
                               compareMethod = compareMethod)
          time = toc()
          saveRDS(vivid,paste0(data,"_", methods[l],"_fold_",i,"_",j,".rds",collapse=""))
          newX = xTrain[,unlist(vivid$optModel) == 1]
          vividGLM = glmnet::cv.glmnet(newX, 
                                       yTrain, 
                                       alpha = 0, 
                                       family = "binomial")
          newXTest = xTest[,unlist(vivid$optModel) == 1]
          vividPred = predict(vividGLM,
                              s = "lambda.1se",
                              newx = newXTest)
          fitted = c(exp(vividPred)/(1+exp(vividPred)))
          features2[[j]] = vivid$optFeatures
        }
        
        if(methods[l] == "rf"){
          tic()
          rf = randomForest::randomForest(x = xTrain,
                                          y = yTrain,
                                          ntree = 500)
          time = toc()
          saveRDS(rf,paste0(data,"_", methods[l],"_fold_",i,"_",j,".rds",collapse=""))
          fitted = predict(rf, 
                           newdata = xTest, 
                           type = "prob")[,2]
        }
        
        if(methods[l] == "boruta"){
          tic()
          boruta = Boruta::Boruta(x = xTrain,
                                  y = yTrain,
                                  doTrace = 0,
                                  maxRuns = 500)
          time = toc()
          
          saveRDS(boruta,paste0(data,"_", methods[l],"_fold_",i,"_",j,".rds",collapse=""))
          newX = xTrain[,boruta$finalDecision == "Confirmed"]
          borutaGLM = glmnet::cv.glmnet(newX, 
                                        yTrain, 
                                        alpha = 0, 
                                        family = "binomial")
          newXTest = xTest[,boruta$finalDecision == "Confirmed"]
          borutaPred = predict(borutaGLM,
                               s = "lambda.1se",
                               newx = newXTest)
          fitted = c(exp(borutaPred)/(1+exp(borutaPred)))
          features2[[j]] = colnames(x)[which(boruta$finalDecision == "Confirmed")]
        }
        
        if(methods[l] == "RFE"){
          tic()
          RFE = sigFeature::sigFeature(X = xTrain,
                                       Y = yTrain)
          time = toc()
          saveRDS(RFE,paste0(data,"_", methods[l],"_fold_",i,"_",j,".rds",collapse=""))
          newX = xTrain[,RFE > (p - topN)]
          RFEGLM = glmnet::cv.glmnet(newX, 
                                     yTrain, 
                                     alpha = 0, 
                                     family = "binomial")
          newXTest = xTest[,RFE > (p - topN)]
          RFEPred = predict(RFEGLM,
                            s = "lambda.1se",
                            newx = newXTest)
          fitted = c(exp(RFEPred)/(1+exp(RFEPred)))
          features2[[j]] = base::colnames(x)[RFE > (p - topN)]
          
        }
        
        modelPred = 1*(fitted > 0.5)
        binary = 1*(yTest == levels(y)[2])
        accVals[j, i] = mean(1*(modelPred == binary))
        roc = pROC::roc(binary, fitted)
        aucVals[j, i] = roc$auc
        timeVals[j, i] = time$toc - time$tic
        
      }
      features[[i]] = features2
    }
    
    output[[l]] = list(
      acc = accVals,
      auc = aucVals,
      time = timeVals,
      features = features
    )
  }
  names(output) = methods
  
  return(output)
  
}

######################

library(RankProd)
library(BioMark)
library(parallel)
library(tidyverse)
library(latticeExtra)
library(biosigner)
library(tictoc)
library(pROC)

data(sacurine)

x <- sacurine$dataMatrix
y <- sacurine$sampleMetadata$gender

sacurine_check = check_data(x = x, y = y, B = 100, k = 10 , repB = 10, methods = c("vivid","boruta","RFE"), topN = 10, data = "sacurine")
saveRDS(sacurine_check,"sacurine_check.rds")

#############################

data("diaplasma")

x <- diaplasma$dataMatrix
y <- diaplasma$sampleMetadata$type

###########################

eval_plot = function(object, repB){
  models = names(object)
  nModels = length(models)
  auc = matrix(data = NA_real_, 
               nrow = repB, 
               ncol = nModels)
  base::colnames(auc) = models
  acc = matrix(data = NA_real_, 
               nrow = repB, 
               ncol = nModels)
  base::colnames(acc) = models
  
  for (i in 1:nModels){
    auc[,i] = base::apply(X = object[[i]]$auc,
                          MARGIN = 2,
                          FUN = mean)
    acc[,i] = base::apply(X = object[[i]]$auc,
                          MARGIN = 2,
                          FUN = mean)
  }
  
  plotDataAUC = reshape2::melt(auc, factorsAsStrings = TRUE)
  plotDataACC = reshape2::melt(acc, factorsAsStrings = TRUE)
  base::colnames(plotDataAUC) = c("Rep", "Method", "Value")
  base::colnames(plotDataACC) = c("Rep", "Method", "Value")  
  plotDataAUC$'Measure' = "AUC"
  plotDataACC$'Measure' = "ACC"
  plotData = rbind(plotDataAUC, plotDataACC)
  output = ggplot2::ggplot(data = plotData,
                           mapping = aes(x = Measure, 
                                         y = Value,
                                         fill = Method)) +
    geom_boxplot()
  
  return(output)
}


###############################


