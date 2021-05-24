###functions to help to get patient features###
get_CNV_rownames = function(patient_CNV_data) {
  CNV_rownames = patient_CNV_data[, 1]
  return(CNV_rownames)
}

get_mu_rownames = function(patient_mu_data) {
  mu_rownames = patient_mu_data[, 1]
  return(mu_rownames)
}

get_GE_rownames = function(patient_GE_data) {
  GE_rownames = patient_GE_data[, 1]
  return(GE_rownames)
}

###get data frame that contains the features of patients and combine the data together###
get_data_for_model = function(GE_rownames, mu_rownames, CNV_rownames, CNV, GE, mu, response) {
  CNV_new = CNV[which(CNV$name %in% CNV_rownames), ]
  GE_new = GE[which(GE$name %in% GE_rownames), ]
  mu_new = mu[which(mu$name %in% mu_rownames), ]
  combine_new = rbind(CNV_new, GE_new, mu_new, response)
  
  combine_new = as.data.frame(combine_new)
  n = combine_new$name
  combine_new = t(combine_new[, -1])
  
  colnames(combine_new) = n
  combine_new = as.data.frame(combine_new)
  return(combine_new)
}


get_data_for_model_onlyGE = function(GE_rownames, mu_rownames, CNV_rownames, CNV, GE, mu, response) {
  
  GE_new = GE[which(GE$name %in% GE_rownames), ]
  
  combine_new = rbind(GE_new, response)
  
  combine_new = as.data.frame(combine_new)
  n = combine_new$name
  combine_new = t(combine_new[, -1])
  
  colnames(combine_new) = n
  combine_new = as.data.frame(combine_new)
  return(combine_new)
}


get_data_for_model2 = function(GE_rownames, mu_rownames, CNV_rownames, cellLines, CNV, GE, mu, response) {
  CNV_new = CNV[which(CNV$name %in% CNV_rownames), ]
  GE_new = GE[which(GE$name %in% GE_rownames), ]
  mu_new = mu[which(mu$name %in% mu_rownames), ]
  
  combine_new = rbind(CNV_new, GE_new, mu_new, response)
  
  combine_new = as.data.frame(combine_new)
  n = combine_new$name
  index<-which(colnames(combine_new) %in% cellLines)
  combine_new = t(combine_new[,index])
  
  colnames(combine_new) = n
  combine_new = as.data.frame(combine_new)
  return(combine_new)
}

get_data_for_model2_onlyGE = function(GE_rownames, cellLines, GE, response) {
  
  GE_new = GE[which(GE$name %in% GE_rownames), ]
  
  
  combine_new = rbind(GE_new, response)
  
  combine_new = as.data.frame(combine_new)
  n = combine_new$name
  index<-which(colnames(combine_new) %in% cellLines)
  combine_new = t(combine_new[,index])
  
  colnames(combine_new) = n
  combine_new = as.data.frame(combine_new)
  return(combine_new)
}

###function to compute MSE###
MSE <- function(obs,pred)
{
  error<-obs-pred
  mean(error^2) 
}

RSquare<-function(obs,pred){
  num<-0
  mn<-mean(obs)
  num<-sum((pred-mn)^2)
  denominator<-var(obs)
  num/denominator
}

##N-fold division
getFoldIndices<-function(L,N){
  size<-round(L/N)
  indices<-c(1:L)  #indices content will be the index values itself
  
  foldIndices<-vector("list", N)
  
  i<-1
  while(i<=N-1){
    tmp<-sample(indices,size = size)
    foldIndices[[i]]<-tmp
    i<-i+1
    a<-which(indices %in% tmp)
    indices<-indices[-a]
  }
  tmp<-sample(indices,size = length(indices))
  foldIndices[[i]]<-tmp
  
  foldIndices
}




##my personal histogram calculator
getHistogram<-function(auc, interval){
  min_auc<-round(min(auc))-1
  max_auc<-round(max(auc))+1
  bin<-c(min_auc)
  bb<-bin[1]+interval
  i<-1
  while(bb<max_auc){
    bin<-c(bin,bb)
    i<-i+1
    #print(bin)
    bb<-bin[i]+interval
    #print(bb)
  }
  bin<-c(bin,max_auc)
  
  binfreq<-c()
  for(i in 2:length(bin)){
    tmp<-auc[which(auc<bin[i])]
    binfreq<-c(binfreq,length(which(tmp>=bin[i-1])))
    #binfreq<-c(binfreq,length(intersect(which(auc>=bin[i-1]), which(auc<bin[i]))))
  }
  list(bin,binfreq)
}



trainOrdinarySVM<-function(TEST){
  
  #drug_index<-1
  TUNELENGTH<-10 
  model11 = get_data_for_model(GE_names, mu_names, CNV_names, CNV, GE, mu, response[drug_index, ])
  #remove NAs
  index<-c()
  for (row in 1:nrow(model11)){
    if (is.na(model11[row, ncol(model11)])==TRUE){
      index<-c(index,row)
    }
  }
  if(length(index)>0) model11<-model11[-index,]
  
  formula = paste(drug_names[drug_index], '~.', sep = '')
  
  sigma11=-1
  C11=-1
  ix11 = sample(1:nrow(model11), nrow(model11) * perc)
  train11 = model11[ix11, ]
  test11 = model11[-ix11, ]
  
  svm11 <- train(as.formula(formula), data = train11, method = "svmRadial", tuneLength = 10)
  sigma11<-svm11$bestTune$sigma
  C11<-svm11$bestTune$C
  
  L<-nrow(TEST) 
  
  #vector containing the final prediction values
  pred<-c()
  
  for(test_index in 1:L){
    
    #test_index<-1
    prd = predict(svm11, TEST[test_index,])
    prd = as.numeric(prd)
    pred<-c(pred,prd)  
    
  } #test_index
  
  pred<-as.numeric(pred)
  pred
  #mse_svmB = MSE(TEST[1:L,ncol(test_normal)],pred)
  #print(c(xline,yline,radius,1/mse_svmB))
  #list(c(xline,yline,radius),c(1/mse_svmB),c(sigma11,C11),ix11)
  
} #trainOrdinarySVM()


