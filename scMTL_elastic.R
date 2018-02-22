args <- commandArgs(trailingOnly=T)
set.seed(1)

features.path <- args[1] ## e.g, /MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/K562/01_K562_newBPs_plus_NFR_minus_divergent.x, but the file should have a header
reponse.path <- args[2] ## e.g, paste("/MMCI/MS/ExpRegulation/work/data/singleCell/K562/transcript/K562_scell_",response.name,"_transcript_expression.txt",sep="")
response.name <- args[3] ## Either "plus" or "minus" to specify which file to be used
feature.name.path <- args[4] ## path to the TF's name
cell.name <- args[5] ## Either K562 or HepG2, at the moment
logTransY <- args[6] ## Logical, whether the response variable (Y) should be log transformed or not
logTransX <- args[7] ## Logical, whether the response variable (Y) should be log transformed or not

if(length(response.name) == 0)
  stop("The argument is missing. Please specify the response type, by calling the script with either \"plus\" or \"minus\"")
print(paste("response:",response.name))

load(paste("/MMCI/MS/ExpRegulation/work/data/singleCell/",cell.name,"/",cell.name,"_single_cell_expr_swapped_transcript.RData",sep=""))
load(paste("/MMCI/MS/ExpRegulation/work/data/singleCell/",cell.name,"/",cell.name,"_hclust_transcript_4.RData",sep=""))
features <- read.table(features.path,header=T,stringsAsFactors=F)

response <- read.table(reponse.path ,header=T,row.names = 1)
response$transcript_id <- NULL
feature.names <- readLines(feature.name.path)
library(pheatmap)
library(glmnet)
library(gplots)
library(parallel)
library(doParallel)
library(ggplot2)
registerDoParallel(8)

############## Functions ##############
train.indiv.models <- function(x,y,alpha_seq){
  K <- ncol(y)
  all.pred.res <- list()
  all.glmnet.res <- list()
  for(i in seq(K)){
    cv.glmnet_res <- mclapply(alpha_seq,function(alpha){cv.glmnet(x=as.matrix(x),y=as.numeric(y[,i]),family="gaussian",alpha=alpha,parallel=T)},mc.cores=20)
    best.error <- Inf
    for(alpha in seq(length(cv.glmnet_res))){
      cvm.best <- min(cv.glmnet_res[[alpha]]$cvm)
      if(best.error > cvm.best){
        best.model <- cv.glmnet_res[[alpha]]
        best_alpha <- alpha_seq[alpha]
        best.error <- cvm.best
      }
    }
    glmnet_res <- glmnet(x=as.matrix(x),y=as.numeric(y[,i]),family="gaussian",lambda=best.model$lambda.min,alpha=best_alpha)
    pred.res <- predict(glmnet_res,as.matrix(x),s=best.model$best.model$lambda.min,alpha=best_alpha)
    all.glmnet.res[[i]] <- glmnet_res
    all.pred.res[[i]] <- pred.res
  }
  return(list(all.pred.res=all.pred.res,all.glmnet.res=all.glmnet.res))
}

#########################################

partition <- function(x,y,percent=0.9){
  train.cnt <- floor(nrow(x)*percent)
  train.idx <- sample(nrow(x))[seq(train.cnt)]
  train.x <- x[train.idx,]
  test.x <- x[-train.idx,]
  if(!is.null(y)){
    train.y <- y[train.idx,]
    test.y <- y[-train.idx,]
  }else{
    train.y <- y[train.idx]
    test.y <- y[-train.idx]
  }
  return(list(train=list(x=train.x,y=train.y),test=list(x=test.x,y=test.y)))
}

#########################################

get.outerCV.partitions <-function(x,y,n.folds=10){
  fold.size <- ceiling(nrow(x)/n.folds)
  partitions <- list();
  for(i in seq(n.folds)){
    start <- (i-1)*fold.size+1;
    end <- i*fold.size;
    if(end > nrow(x))
      end <- nrow(x);
    test.idx <- seq(start,end)
    train.idx <- seq(nrow(x))[-test.idx]
    if(is.null(dim(y))){
      partitions[[i]] <- list(train=list(x=x[train.idx,],y=y[train.idx]),test=list(x=x[test.idx,],y=y[test.idx]))
    }
    else{
      partitions[[i]] <- list(train=list(x=x[train.idx,],y=y[train.idx,]),test=list(x=x[test.idx,],y=y[test.idx,]))
    }
  }
  return(partitions)
}

#########################################

train_elasticNet <- function(fold,alpha_seq,fold_cnt){
  if(nrow(fold$train$x) > 50){
    nfold <- 10
  }else{
    #nfold <- nrow(fold$train$x) - 1
    nfold <- 10
  }
  ################################
  ################################
  ### MTL ###

  best.error <- Inf ## To keep track of the best model out of several alphas
  cv.models <- mclapply(alpha_seq,function(alpha){
    cv.glmnet(x=as.matrix(fold$train$x),y=as.matrix(fold$train$y),family="mgaussian",parallel=T,alpha=alpha,nfolds=nfold)
  },mc.cores=20) ## Performs MTL on the data above in parallel
  
  #### Choose the best model from "models"
  for(m in seq(length(alpha_seq))){
    cvm.best <- min(cv.models[[m]]$cvm)
    if(best.error > cvm.best){
      best.model.idx <- m
      best.error <- cvm.best
      best.alpha <- alpha_seq[m]
    }
  }

  glmnet_res <- glmnet(x=as.matrix(fold$train$x),y=as.matrix(fold$train$y),family="mgaussian",lambda=cv.models[[best.model.idx]]$lambda.min,alpha=best.alpha) ## Perform the MTL on the training set using the optimal parameters lambda and alpha

  pred.res <- predict(glmnet_res,as.matrix(fold$train$x)) ## Get the prediction on the training set
  pred.res.test <- predict(glmnet_res,as.matrix(fold$test$x)) ## Get the prediction on the test set


  ## Compute covariances
  train.cov <- cov(as.matrix(fold$train$y), pred.res[,,1])
  test.cov <- cov(as.matrix(fold$test$y), pred.res.test[,,1])

  ## Compute RMSEs
  train.RMSE <- sqrt(sum((as.matrix(fold$train$y) - pred.res[,,1])^2) / length(pred.res))
  test.RMSE <- sqrt(sum((as.matrix(fold$test$y) - pred.res.test[,,1])^2) / length(pred.res.test))
  print(paste("train RMSE=",train.RMSE))
  print(paste("test RMSE=",test.RMSE))

  ### Plot all predictions vs all responses for a cluster
  ## Train
  pred.res.vectorized <- as.vector(pred.res[,,1])
  y.vectorized <- as.vector(as.matrix(fold$train$y))
  cor.train <- round(cor(pred.res.vectorized,y.vectorized),3)
  train.diff <- pred.res[,,1] - as.matrix(fold$train$y)

  ## Test
  pred.res.test.vectorized <- as.vector(pred.res.test[,,1])
  y.test.vectorized <- as.vector(as.matrix(fold$test$y))
  cor.test <- round(cor(pred.res.test.vectorized,y.test.vectorized),3)
  test.diff <- pred.res.test[,,1] - as.matrix(fold$test$y)

  ### Compute correlation per task
  ## Train
  train.cors <- sapply(seq(ncol(pred.res[,,1])),function(i)cor(pred.res[,i,1], fold$train$y[,i]))

  ## Test
  test.cors <- sapply(seq(ncol(pred.res.test[,,1])),function(i)cor(pred.res.test[,i,1], fold$test$y[,i]))

  ## Create data.frame to show train and test correlations as barplots (each bar for a task
  df <- data.frame(cor=c(train.cors,test.cors),dataset=c(rep("train",length(train.cors)),rep("test",length(test.cors))),cell=c(seq(length(train.cors)),seq(length(test.cors))))

  ################################
  ################################
  ### Individual models ###
  
  indiv.preds <- train.indiv.models(fold$train$x,fold$train$y,alpha_seq) ## Get the predictions for each cell
  indiv.betas.mat <- NULL; ## To aggegates the beta coefficients spitted from the individual models

  print(length(indiv.preds$all.glmnet.res))

  ## assign the coefficients to indiv.betas.mat
  for(cell in seq(length(indiv.preds$all.glmnet.res)))
    indiv.betas.mat <- cbind(indiv.betas.mat,as.numeric(indiv.preds$all.glmnet.res[[cell]]$beta))

  models[[i]] <- glmnet_res ## To keep the MTL model obtained for each cluster

  betas <- NULL; ## MTL betas
  ## Combine the beta coefficients associated to each task
  for(j in seq(length(glmnet_res$beta)))
    betas <- cbind(betas,as.numeric(glmnet_res$beta[[j]]))

  pdf(paste("/MMCI/MS/ExpRegulation/work/data/singleCell/",cell.name,"/scMTL_hclust_",k,"_cluster_",i,"_response_",response.name,"_fold_",fold_cnt,".pdf",sep=""))
  if(sum(betas) == 0){
    print(paste("the MTL for cluster",i,"was not successfully trained. All the coefficients are zero!"))
    next;
  }
  heatmap.2(betas,Rowv=F,Colv=F,dendrogram="n",trace="n",col=my_palette1,key=T,cexRow=.5,main=paste("MTL coefficients\n clust_",i,"_lambda*=",round(cv.models[[best.model.idx]]$lambda.min,3)," alpha*=",best.alpha,sep=""),labRow=colnames(fold$train$x))
  ###
  #heatmap.2(train.cov, dendrogram="n",trace="n",col=my_palette1,key=T,main=paste("Covariance btw. predictions & measured\n training data; clust_",i,"",sep=""))
  #heatmap.2(test.cov, dendrogram="n",trace="n",col=my_palette1,key=T,main=paste("Covariance btw. predictions & measured\n test data; clust_",i,"",sep=""))
  ###
  heatmap.2(train.diff, dendrogram="n",trace="n",col=my_palette1,key=T,main=paste("train;pred - measured\nclust_",i,"\nRMSE=",train.RMSE,sep=""))
  heatmap.2(test.diff, dendrogram="n",trace="n",col=my_palette1,key=T,main=paste("test; pred - measured\nclust_",i,"\nRMSE=",test.RMSE,sep=""))
  ###
  plot(pred.res.vectorized,y.vectorized,pch=20,col="blue",main=paste("train, cor=",cor.train),xlab="MTL prediction",ylab="measured scRNA")
  abline(seq(range(pred.res.vectorized)),seq(range(y.vectorized)),col="red")
  ###
  plot(pred.res.test.vectorized,y.test.vectorized,pch=20,col="blue",main=paste("test, cor=",cor.test),xlab="MTL prediction",ylab="measured scRNA")
  abline(seq(range(pred.res.test.vectorized)),seq(range(y.test.vectorized)),col="red")
  ###
  print(ggplot(data=df,aes(x=cell,y=cor,fill=factor(dataset))) + geom_bar(stat="identity",position="dodge") + ggtitle(cell.name) + theme(plot.title=element_text(hjust=.5)))
  if(F){
    for(j in seq(feature.cnt)){
      beta_hm <- betas[seq((j-1)*bin.cnt+1,j*bin.cnt),]
      beta_hm_indiv <- indiv.betas.mat[seq((j-1)*bin.cnt+1,j*bin.cnt),]
      if(sum(sum(beta_hm))){
        heatmap.2(beta_hm,Rowv=F,Colv=F,dendrogram="n",trace="n",col=my_palette1,key=T,cexRow=.5,main=paste("MTL; clust_",i,"_lambda*=",round(cv.models[[best.model.idx]]$lambda.min,3)," alpha*=",best.alpha,sep=""),labRow=paste(feature.names[j],bins,sep="_"))
      }
      if(sum(sum(beta_hm_indiv)))
        heatmap.2(beta_hm_indiv,Rowv=F,Colv=F,dendrogram="n",trace="n",col=my_palette1,key=T,cexRow=.5,main=paste("individual model",sep=""),labRow=paste(feature.names[j],bins,sep="_"))
    }
  }
  write.table(glmnet_res$a0,paste("/MMCI/MS/ExpRegulation/work/data/singleCell/",cell.name,"/scMTL_a0_",i,"_response_",response.name,"_fold_",fold_cnt,".txt",sep=""))
  dev.off()

  return(list(indiv.preds = indiv.preds, models = models, betas = betas, MTL_pred.res = pred.res, MTL_pred.res.test = pred.res.test))
}
############ End of Functions ############
##########################################
##########################################

outerCV_folds <- 5
feature.cnt <- length(feature.names)
bin.cnt <- ncol(features)/feature.cnt
#feature.names <- c("H3K4me1","H3K4me3","H3K27me3","H3K36me3","H3K9me3","H3K27ac")
bins <- seq(-floor(bin.cnt/2),floor(bin.cnt/2))
models <- list() ## To keep the MTL models run for each cluster
all.betas <- list() ## A list, of size equal to the number of clusters, to keep the beta matrices computed for each cluster
my_palette1 <- colorRampPalette(c("blue","white","red"))(256)
indiv.preds <- list() ## To keep the glmnet models run for each cell separately. Used to compare the MTL with individual models
tf.names <- colnames(features)


if(logTransY == "TRUE"){
  response <- log2(1 + response)
}

if(logTransX == "TRUE"){
  features <- log2(1 + features)
}
alpha_seq <- seq(0,1,.05)



for(i in seq(k)){

  print(i) ## Print the cluster number
  response_i <- response[which(clusts == i),] ## Select the response (scell expressions) corresponding to the cluster i
  features_i <- features[which(clusts == i),] ## Select the input features corresponding to the cluster i

  outer_CV_partitions <- get.outerCV.partitions(features_i,response_i,n.folds=outerCV_folds) ## Generate the outer CV folds

  scMTL.outerCV <- lapply(seq(length(outer_CV_partitions)),function(fold)train_elasticNet(outer_CV_partitions[[fold]],alpha_seq,fold))


  all.ocv.betas <- NULL
  folds_number <- NULL;
  test.responses <- NULL
  pred.res.vectorized <- NULL
  pdf(paste("/MMCI/MS/ExpRegulation/work/data/singleCell/",cell.name,"beta_distribution_outerCV_clust_",i,".pdf",sep=""))
  par(mfrow=c(2,3))
  for(f in seq(outerCV_folds)){
    nonZeroCoef.idx <- which(rowSums(abs(scMTL.outerCV[[f]]$betas)) != 0)
    boxplot(t(scMTL.outerCV[[f]]$betas[nonZeroCoef.idx,]),names=tf.names[nonZeroCoef.idx],las=2,ylab="beta distribution across cells",main=paste("outer CV fold",f))

    all.ocv.betas <- cbind(all.ocv.betas, scMTL.outerCV[[f]]$betas)

    test.responses <- c(test.responses,as.vector(as.matrix(outer_CV_partitions[[f]]$test$y)))
    pred.res.vectorized <- c(pred.res.vectorized,as.vector(scMTL.outerCV[[f]]$MTL_pred.res.test[,,1]))
    folds_number <- c(folds_number,rep(f,prod(dim(outer_CV_partitions[[f]]$test$y))));
  }

  all.betas[[i]] <- all.ocv.betas
  nonZeroCoef.idx <- which(rowSums(abs(all.betas[[i]])) != 0)
  boxplot(t(all.betas[[i]][nonZeroCoef.idx,]),names=tf.names[nonZeroCoef.idx],las=2,ylab="beta distribution across cells",main=paste("accumulation of all outer CV folds"))

  df <- data.frame(predicted=pred.res.vectorized,measured=test.responses,fold=folds_number)


  cor_test_allFolds <- cor(test.responses,pred.res.vectorized)
  print(ggplot(df,aes(x=predicted,y=measured,fill=factor(fold))) + geom_point(aes(colour= factor(fold),shape= factor(fold))) + theme_bw() + xlim(range(c(df$predicted,df$measured))) + ylim(range(c(df$predicted,df$measured))) + ggtitle(paste("corr=",round(cor_test_allFolds,3))))
  dev.off()

  save(scMTL.outerCV, indiv.preds,models,all.betas,file=paste("/MMCI/MS/ExpRegulation/work/data/singleCell/",cell.name,"/scMTL_hclust_",k,"_cluster_",i,"_response_",response.name,".RData",sep=""))
}

