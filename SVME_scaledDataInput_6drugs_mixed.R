library(e1071)
library(caTools)
library(caret)
library(data.table)
library(miscTools)
library(foreach)
library(doParallel)
library(iterators)
library(doMC)
library(cowplot)

source("functions_3C.R")
source("SVMEnsemble_12.R")
source("SVMEnsemble_11.R")
source("SVMEnsemble_10.R")
source("SVMEnsemble_9.R")
source("SVMEnsemble_8.R")
source("SVMEnsemble_7.R")
source("SVMEnsemble_6.R")
source("SVMEnsemble_5.R")
source("SVMEnsemble_4.R")
source("SVMEnsemble_3.R")
source("SVMEnsemble_2.R")


args = commandArgs(trailingOnly=TRUE)

#system.time({

GE <- fread("GE_Standardized.csv")
GE<-as.data.frame(GE)
response = fread("response.csv")
mu <- fread("mu.csv")
CNV <- fread("CNV.csv")
Targets = read.csv("Targets.csv")

names_CNV = colnames(CNV)
names_GE = colnames(GE)
names_mu = colnames(mu)
names_response = colnames(response)

names_vector = c()

#Cell line Annotation
####Get shortest names in four data forms###
for (i in 2:611) {
  if (names_CNV[i] == names_GE[i] && names_CNV[i] == names_mu[i] && names_CNV[i] == names_response[i]) {
    names_vector = c(names_vector, names_CNV[i])
  } else {
    str = c(names_CNV[i], names_GE[i], names_mu[i], names_response[i])
    str_len = c(nchar(names_CNV[i]), nchar(names_GE[i]), nchar(names_mu[i]), nchar(names_response[i]))
    min = min(str_len)
    index = which(str_len == min)
    if (length(index) > 1) {
      index = index[1]
    }
    names_vector = c(names_vector, str[index])
  } 
}

###change all names to same###
names_vector= c("name", names_vector)
names(CNV) = names_vector
names(GE) = names_vector
names(response) = names_vector
names(mu) = names_vector

CNV = as.data.frame(CNV)

GE = as.data.frame(GE)
#remove ACTB, other cannot center and scale during PCA
delACTB<-which(GE$name=='ACTB')
GE<-GE[-delACTB,]

mu = as.data.frame(mu)
response = as.data.frame(response)

Target_names = Targets$Target_Symbol_321
intersect_name = intersect(Target_names, GE$name)
intersect_name_index = which(GE$name %in% intersect_name)


GEvalues<-GE[,-1]
GEvalues<-(GEvalues-min(GEvalues))/(max(GEvalues)-min(GEvalues))
GE[,-1]<-GEvalues
pcaInput<-t(GEvalues)
pcaInput<-as.data.frame(pcaInput)
rownames(pcaInput)<-colnames(GE)[-1]
colnames(pcaInput)<-GE$name


CNV_test = c()
GE_test = GE[intersect_name_index,]

#mu_test = mu[1:3,]
mu_test = c()
CNV_names = get_CNV_rownames(CNV_test)
CNV_names = as.vector(CNV_names)
GE_names = get_GE_rownames(GE_test)
GE_names = as.vector(GE_names)
mu_names = get_mu_rownames(mu_test)


response$name =  gsub("-", "", response$name)
##delete white spaces in names of response##
response$name =  gsub(" ", "", response$name)
drug_names = response[, 1]

for (i in 1:481) {
  drug_names[i] = paste("DRUG", drug_names[i], sep = "")
}
drug_names[436] = 'DRUG1S3RRSL3'
response$name = drug_names


index<-c()
for (row in 1:nrow(pcaInput)){
  if (is.na(pcaInput[row, ncol(pcaInput)])==TRUE){
    index<-c(index,row)
  }
}
if(length(index)>0) pcaInput<-pcaInput[-index,]

myPCA1<-prcomp(pcaInput, scale. = T, center = T)
pc1<-myPCA1$x[,1]


reals<-c()
predicteds<-c()
titles<-c()
mse_ensemble<-c()
mse_simple<-c()
drug_index<-as.integer(args[1])


print(paste('drug_index: ',drug_index,sep = ''))
drug<-response[drug_index,-1]
delD<-which(drug=='NaN')
drug<-drug[-delD]  #removing cell-line where AUC value is NaN
auc<-unlist(drug)
y<-pc1[-delD]   #removing those cell-line PC1 values whose AUC value is NaN
yName<-rownames(pcaInput)[-delD]


ITERATION = 1
DEBUG<-FALSE



ALLDATA = get_data_for_model2(GE_names, mu_names, CNV_names, yName, CNV, GE, mu, response[drug_index, ])

DIST<-matrix(nrow=length(yName),ncol=2)
DIST<-as.data.frame(DIST)
DIST[,1]<-auc
DIST[,2]<-y
rownames(DIST)<-yName
colnames(DIST)<-c('auc','pc1')

TEST = ALLDATA

set.seed(12345)


#percentage of data to be used for testing
percTesting<-0.25

#SELECTION=0: means using the maximum avg Spearman correlation with all cell-lines across all blocks
#SELECTION=1: means using the maximum avg Spearman correlation with neighbors across all blocks
#SELECTION=2: means using avg Spearman correlation with neighbors as weights and returning weighted average
#SELECTION=3: means using only number of neughbors
#SELECTION=4: means using number of neighbors as weights and returning weighted average
SELECTION<-1

#type of cross-validation CV:0 means ordinary CV, CV:1 means repeated CV
CV<-0 
fold<-5 #no. of cross-validation folds
#if CV=1, number of repetitions
N<-5




#where to store results
PATH<-args[2]


formula = paste(drug_names[drug_index], '~.', sep = '')

starttime<-Sys.time()


system.time({
  svmTrI<-sample(1:nrow(ALLDATA),(1-percTesting)*nrow(ALLDATA))
  svmTr<-ALLDATA[svmTrI,]
  svmTe<-ALLDATA[-svmTrI,]
  svm <- train(as.formula(formula), data = svmTr, method = "svmRadial", tuneLength = 10)
  train_svm_mse<-min(svm$results$RMSE)
  svm_pred<-predict(svm,svmTe)
  svm_pred<-as.numeric(svm_pred)
  test_svm_mse<-MSE(svmTe[,ncol(svmTe)],svm_pred)
  tbw<-data.frame(c(1),train_svm_mse,test_svm_mse,mean(train_svm_mse,test_svm_mse))
  colnames(tbw)<-c('Iteration/Fold','Train_MSE','Test_MSE','Avg_Train_Test_MSE')
  write.csv(tbw,paste(PATH,'SVMEnsemble_1_',drug_names[drug_index],'.csv',sep = ''), row.names = F, quote = F)
  })
simple<-predict(svm,TEST)
tbp<-as.data.frame(simple)
colnames(tbp)<-'svm'
tbp$real<-TEST[,ncol(TEST)]
tbp$error<-abs(tbp$svm-tbp$real)
tbp$pc1<-DIST$pc1
rownames(tbp)<-rownames(TEST)


S<-sort.int(tbp$real, index.return = T)
tbp1<-tbp[S$ix,]

X<-tbp1[,2:3]


bestMse<-mean(train_svm_mse,test_svm_mse)
bestC<-1

testMSE<-c(bestMse)
clstr<-c(1)

bestParams<--1
bestPartition<-NULL

for(clusters in 2:12){
  
  clstr<-c(clstr,clusters)
  kme<-kmeans(X,clusters)
  tbp1$clid<-kme$cluster
  blocks<-vector("list", clusters)
  for(k in 1:clusters){
    indx<-which(kme$cluster==k)
    blocks[[k]]<-rownames(X)[indx]
  }
  print(paste('k=',clusters))
  
  if(clusters==12) {
    system.time({
      params12<-SVMEnsemble_12(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params12[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params12
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==11) {
    system.time({
      params11<-SVMEnsemble_11(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params11[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params11
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==10) {
    system.time({
      params10<-SVMEnsemble_10(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params10[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params10
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==9) {
    system.time({
      params9<-SVMEnsemble_9(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params9[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params9
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==8) {
    system.time({
      params8<-SVMEnsemble_8(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params8[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params8
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==7) {
    system.time({
      params7<-SVMEnsemble_7(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params7[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params7
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==6) {
    system.time({
      params6<-SVMEnsemble_6(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params6[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params6
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==5) {
    system.time({
      params5<-SVMEnsemble_5(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params5[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params5
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==4) {
    system.time({
      params4<-SVMEnsemble_4(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params4[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params4
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==3) {
    system.time({
      params3<-SVMEnsemble_3(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params3[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params3
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  } else if(clusters==2) {
    system.time({
      params2<-SVMEnsemble_2(blocks,percTesting,SELECTION,fold,N,CV,DEBUG,ALLDATA,DIST,y,yName,PATH)
      mse<-params2[[1]][1]
      if(mse<bestMse){
        bestC<-clusters
        bestParams<-params2
        bestMse<-mse
        bestPartition<-kme$cluster
      }})
  }
  testMSE<-c(testMSE,mse)
  
}#clusters

stoptime<-Sys.time()
print(paste('Model_Setup_Time: ',difftime(stoptime, starttime, units = "secs")[[1]],sep=''))

minIndx<-which(testMSE==min(testMSE))
drug_names[drug_index]
testMSE[minIndx]
clstr[minIndx]
bestMse
bestC

fn1<-paste(PATH,drug_names[drug_index],'_bestParams_',bestC,'clusters.rds',sep = '')
saveRDS(bestParams, file = fn1)
fn2<-paste(PATH,drug_names[drug_index],'_bestPartitions_',bestC,'clusters.rds',sep = '')
saveRDS(bestPartition, file = fn2)
#} #drug_index 

tbp1$clid<-bestPartition

p1<-ggplot(tbp1, aes(real,error, color=factor(clid))) + geom_point() + scale_color_manual(values=c("blue", "darkolivegreen4","darkorchid3", "orange","red","black","cyan","pink","green","maroon","gray","peachpuff")) + theme_classic() + labs(x = 'Real AUC' , y = 'Predicted AUC Error')
p2<-ggplot(tbp1, aes(real,pc1, color=factor(clid))) + geom_point() + scale_color_manual(values=c("blue", "darkolivegreen4","darkorchid3", "orange","red","black","cyan","pink","green","maroon","gray","peachpuff")) + theme_classic() + labs(x = 'Real AUC' , y = 'PC1')


blocks<-vector("list", bestC)
for(k in 1:bestC){
  indx<-which(bestPartition==k)
  blocks[[k]]<-rownames(X)[indx]
}

if(bestC==12) {
  system.time({SVMEnsemble_12_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_12_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==11) {
  system.time({SVMEnsemble_11_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_11_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==10) {
  system.time({SVMEnsemble_10_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_10_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==9) {
  system.time({SVMEnsemble_9_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_9_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==8) {
  system.time({SVMEnsemble_8_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_8_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==7){
  system.time({SVMEnsemble_7_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_7_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==6) {
  system.time({SVMEnsemble_6_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_6_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==5) {
  system.time({SVMEnsemble_5_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_5_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==4) {
  system.time({SVMEnsemble_4_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_4_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==3) {
  system.time({SVMEnsemble_3_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_3_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
} else if(bestC==2) {
  system.time({SVMEnsemble_2_Replicate(bestParams,blocks,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
  system.time({pred<-SVMEnsemble_2_Predict(bestParams,blocks,TEST,SELECTION,DEBUG,ALLDATA,DIST,y,yName)})
}


tbp$pred<-pred
tbp$predError<-abs(tbp$pred-tbp$real)

combined_plot<-data.frame(rep(tbp$real,2))
names(combined_plot)<-'real'
combined_plot$pred<-c(tbp$svm,tbp$pred)
combined_plot$error<-c(tbp$error,tbp$predError)
combined_plot$model<-c(rep('original_svm',nrow(tbp)),rep('ensemble_svm',nrow(tbp)))

p3<-ggplot(combined_plot, aes(real,error, color=factor(model), shape = factor(model))) + geom_point() + scale_color_manual(values=c('original_svm' = "red", 'ensemble_svm' = 'green')) + theme_classic() + labs(x = 'Real AUC' , y = 'Predicted AUC Error')
p4<-ggplot(data = combined_plot , aes(x = real, y = pred, color = factor(model), shape = factor(model))) + geom_point() + geom_abline(slope=1, intercept=0) +
  theme_classic() + labs(x = 'Real AUC' , y = 'Predicted AUC') +
  scale_color_manual(values=c('original_svm' = "red", 'ensemble_svm' = 'green')) +
  xlim(0, max(combined_plot$real,combined_plot$pred)) + ylim(0, max(combined_plot$real,combined_plot$pred))

png(filename = paste(PATH,drug_names[drug_index],"_all4Figs_",bestC,"clusters.png",sep=""), width = 10, height = 10, units = "in", res = 200)
plot_grid(p1,p2,p3,p4, ncol = 2, labels = "AUTO", label_size = 16)
dev.off()


##########################adding code for Kmeans performance##########################

K<-c()
yAvgTestMSE<-c()
yAvgTTMSE<-c()

d<-drug_index

for(k in 1:12){
  K<-c(K,k)
  fn<-paste(PATH,'SVMEnsemble_',k,'_',drug_names[d],'.csv',sep='')
  df<-fread(fn)
  df<-as.data.frame(df)
  yAvgTestMSE<-c(yAvgTestMSE,mean(df$Test_MSE))
  yAvgTTMSE<-c(yAvgTTMSE,mean(df$Avg_Train_Test_MSE))
}

tbp2<-data.frame(K,yAvgTestMSE,yAvgTTMSE)
#p<-ggplot(tbp2, aes(K,yAvgTestMSE)) + geom_point(color="black") + geom_line(color="black") + theme_classic() + labs(x = 'K' , y = 'Mean Test MSE')
p5<-ggplot(tbp2, aes(K,yAvgTTMSE)) + geom_point(color="black") + geom_line(color="black") + theme_classic() + labs(x = 'K' , y = 'Mean (Train+Test) MSE')





png(filename = paste(PATH,drug_names[drug_index],"_all5Figs_",bestC,"clusters.png",sep=""), width = 10, height = 10, units = "in", res = 200)
plot_grid(p5,p1,p2,p3,p4, ncol = 2, labels = "AUTO", label_size = 16)
dev.off()

if(bestC==1){
  p6<-ggplot(tbp, aes(real,svm)) + geom_point(color="red") + geom_abline(slope=1, intercept=0) + theme_classic() + labs(x = 'Real AUC' , y = 'Predicted AUC') + xlim(0, max(tbp$real,tbp$svm)) + ylim(0, max(tbp$real,tbp$svm))
  png(filename = paste(PATH,drug_names[drug_index],"_all2Figs_",bestC,"clusters.png",sep=""), width = 10, height = 10, units = "in", res = 200)
  plot_grid(p5,p6, ncol = 2, labels = "AUTO", label_size = 16)
  dev.off()
}
#a+(A-min(A))*(b-a)/(max(A)-min(A))
