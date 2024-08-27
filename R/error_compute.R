true_rm = function(qtl_list,snp_position_list,matrix_position_list,radius=NULL) {
  if(is.null(radius)){radius = (matrix_position_list[[1]][2] - matrix_position_list[[1]][1])/2}
  lapply(1:length(matrix_position_list),function(chr){
    qtl_list_chr = lapply(qtl_list,function(qtl){
      if(qtl$chromosome ==chr) list(position = snp_position_list[[chr]][qtl$position],effect = qtl$y_effect)
      else {list(position = 0 ,effect=rep(0,length(qtl$y_effect)))}
    })
    lapply(matrix_position_list[[chr]],function(center){
      window = c(center - radius,center+radius)
      effect_in_window = lapply(qtl_list_chr,function(qtl){
        if((qtl$position <= window[2] )& (qtl$position >= window[1])) {qtl$effect}
        else {0}
      })
      result = var(Reduce("+",effect_in_window))
      if(is.na(result)) {return(0)}
      else {return(result)}
    })
  })
}

second_true_rm = function(qtl_list,snp_position_list,matrix_position_list,radius=NULL) {
  if(is.null(radius)){radius = (matrix_position_list[[1]][2] - matrix_position_list[[1]][1])/2}
  lapply(1:length(qtl_list),function(origin){
    lapply(1:length(matrix_position_list),function(chr){
      qtl_list_chr = lapply(qtl_list[[origin]],function(qtl){
        if(qtl$chromosome ==chr) list(position = snp_position_list[[chr]][qtl$position],effect = qtl$y_effect)
        else {list(position = 0 ,effect = rep(0,length(qtl$y_effect)))}
      })
      lapply(matrix_position_list[[chr]],function(center){
        window = center + c(-radius,radius)
        effect_in_window = lapply(qtl_list_chr,function(qtl){
          if(qtl$position <= window[2] & qtl$position >= window[1]) {qtl$effect}
          else {rep(0,length(qtl$effect))}
        })
        var(Reduce("+",effect_in_window))
      })
    })
  })
}

pred_rm = function(variance_prediction,matrix_position_list,radius=NULL){
  if(is.null(radius)){radius = (matrix_position_list[[1]][2] - matrix_position_list[[1]][1])/2}
  lapply(1:length(variance_prediction),function(chr){
    lapply(1:length(variance_prediction[[chr]]),function(i){
      center = matrix_position_list[[chr]][i]
      window = center + c(-radius,radius)
      window = matrix_position_list[[chr]]>=window[1] & matrix_position_list[[chr]]<=window[2]
      result = sum(variance_prediction[[chr]][window])
      if(is.na(result)) {return(0)}
      else {return(result)}
    })
  })
}

get_bglr_pred = function(y,X,matrix_position_list,snp_position_list,chromosome_vector,chr_lists){
  pb=0.2 # probIn 0.1/0.01/0.001 ici 1marker/1000 dans la distrib non nulle
  ct=2 # counts 1/2/5
  ni=10000 #nIter 60K
  bi=2000 #burnIn
  th=2 #thin 12 nbr iterations à garder 1/12
  #calculs initialisation
  DF<-5
  S<-var(y)/2*(DF-2)

  idx=paste("pb",pb,"_ct",ct,"_ni",ni,"_bi",bi,"_th",th, sep="")
  #Selection
  ETA<-list(
    PA = list(X=X,model='BayesC', lambda=25, type='gamma', saveEffects=TRUE,
              rate=1e-4,shape=0.55,probIn = pb, counts = ct))
  fit_bglr <-BGLR::BGLR(y=y,ETA=ETA,nIter=ni,burnIn=bi,thin=th,
                  df0=DF,S0=S, verbose = FALSE,saveAt="~/BGLR/")
  beta_hat = fit_bglr$ETA$PA$b
  return(bglr_rm(beta_hat,snp_position_list,matrix_position_list,chromosome_vector,chr_list))
}



bglr_rm = function(bglr_pred,snp_position_list,matrix_position_list,chromosome_vector,chr_list,radius=NULL) {
  if(is.null(radius)){radius = (matrix_position_list[[1]][2] - matrix_position_list[[1]][1])/2}
  lapply(1:length(snp_position_list),function(chr){
    bglr_chr = bglr_pred[chromosome_vector == chr]
    sapply(matrix_position_list[[chr]],function(center){
      window = center + c(-1,1)*radius
      window = snp_position_list[[chr]]>=window[1] & snp_position_list[[chr]]<=window[2]
      effect = scale(chr_list[[chr]][,window])%*%bglr_chr[window]
      result = var(effect)
      if(is.na(result)){result=0}
      return(result)
    })
  })
}



sample_error_compute = function(variance_prediction_list,qtl_list,snp_position_list,matrix_position_list,radius=NULL,k_max=10) {
  true_fun = unlist(true_rm(qtl_list,snp_position_list,matrix_position_list,radius))
  positives = true_fun>0
  lapply(variance_prediction_list,function(pred){
    pred_fun = unlist(pred_rm(pred,matrix_position_list,radius))
    pred_cat = ceiling(rank(pred_fun) * k_max /length(pred_fun))
    error_function = pred_fun - true_fun
    overestimation = abs(sum(error_function[error_function>0]))/sum(pred_fun)
    underestimation = abs(sum(error_function[error_function<0]))/sum(true_fun)
    misestimation = sum(abs(pred_fun/sum(pred_fun) - true_fun/sum(true_fun)))/2
    correlation = cor(true_fun,pred_fun)
    mcc_result = numeric(k_max)
    tpr = numeric(k_max)
    fpr = numeric(k_max)
    for(k in 1:k_max){
      pred_cat_k = pred_cat>=k
      tpr[k] = sum(positives & pred_cat_k) / sum(positives)
      fpr[k] = sum(!positives & pred_cat_k) / sum(!positives)
      mcc_result[k] = mltools::mcc(positives,pred_cat_k)
    }
    list(misestimation=misestimation,underestimation=underestimation,overestimation=overestimation,correlation=correlation,mcc=mcc_result,tpr=tpr,fpr=fpr)
  })
}

second_sample_error_compute = function(variance_prediction_list,qtl_list,snp_position_list,matrix_position_list,radius=NULL,k_max=10) {
  true_rm_list = second_true_rm(qtl_list,snp_position_list,matrix_position_list,radius)
  true_fun = unlist(Reduce("+",lapply(true_rm_list,unlist)))
  positives = true_fun>0
  lapply(variance_prediction_list,function(pred){
    pred_fun = unlist(pred_rm(pred,matrix_position_list,radius))
    pred_cat = ceiling(rank(pred_fun) * k_max /length(pred_fun))
    error_function = pred_fun - true_fun
    overestimation = abs(sum(error_function[error_function>0]))/sum(pred_fun)
    underestimation = abs(sum(error_function[error_function<0]))/sum(true_fun)
    misestimation = sum(abs(pred_fun/sum(pred_fun) - true_fun/sum(true_fun)))/2
    correlation = cor(true_fun,pred_fun)

    origin_underestimation = sapply(true_rm_list,function(true_rm_origin){
      true_rm_origin_fun = unlist(true_rm_origin)
      error_function_origin = pred_fun - true_rm_origin_fun
      abs(sum(error_function_origin[error_function_origin<0]))/sum(true_rm_origin_fun)
    })
    mcc_result = numeric(k_max)
    tpr = numeric(k_max)
    fpr = numeric(k_max)
    for(k in 1:k_max){
      pred_cat_k = pred_cat>=k
      tpr[k] = sum(positives & pred_cat_k) / sum(positives)
      fpr[k] = sum(!positives & pred_cat_k) / sum(!positives)
      mcc_result[k] = mltools::mcc(positives,pred_cat_k)
    }
    names(origin_underestimation) = c("Y","L","D","SNP")
    list(misestimation=misestimation,underestimation=underestimation,overestimation=overestimation,correlation=correlation,mcc=mcc_result,origin_underestimation=origin_underestimation,tpr=tpr,fpr=fpr)
  })
}


plot_sim = function(simulation,s,methods=c("Horseshoe","Fused","Fusion")) {
  sample = simulation$samples[[s]]
  qtl_position_list = sapply(sample$qtl_list,function(qtl){
    simulation$position_list[[qtl$chromosome]][qtl$position]
  })
  qtl_size_list = sapply(sample$qtl_list,function(qtl)var(qtl$y_effect))
  par(mfrow=c(length(methods),1))
  for(m in 1:length(methods)){
    pred = sample$result[[m]]$lambda2_median
    plot(unlist(simulation$matrix_position_list),unlist(pred),type="l",col="black",main=paste0(methods[m]),ylab="Variance estimée",xlab="Position génétique (cM)")
    abline(v=simulation$chr_start,col="gray",lty=3)
    points(qtl_position_list,qtl_size_list,pch=19)
  }
}

plot_error_sim = function(simulation,errors,s,radius,methods=c("Horseshoe","Fused","Fusion","BGLR")) {
  preds = errors[[s]]$preds
  matrix_positions = unlist(simulation$matrix_position_list)
  qtl_list = simulation$samples[[s]]$qtl_list
  qtl_data_frame = data.frame()
  true_rm_fun = unlist(true_rm(qtl_list,simulation$position_list,simulation$matrix_position_list,radius=radius))
  for(qtl in qtl_list){
    qtl_data_frame = rbind(qtl_data_frame,data.frame(position=simulation$position_list[[qtl$chromosome]][qtl$position]))
  }
  par(mfrow=c(max(1,min(2,length(methods)%/%2)),max(1,length(methods)%/%2)))
  for(m in 1:length(methods)) {
    pred_rm_fun = unlist(pred_rm(preds[[m]],simulation$matrix_position_list,radius=radius))
    plot(matrix_positions,pred_rm_fun,type="l",main = paste0(methods[m]),ylab="Variance estimée",xlab="Position génétique (cM)",ylim=c(0,max(max(pred_rm_fun),max(true_rm_fun))))
    lines(matrix_positions,true_rm_fun,col="blue",lty=3)
    points(qtl_data_frame$position,rep(0,nrow(qtl_data_frame)),pch=19)
    abline(v=simulation$chr_start,col="gray",lty=3)
  }
}

second_plot_error_sim = function(simulation,errors,s,radius,methods=c("Horseshoe","Fused","Fusion","BGLR")) {
  preds = errors[[s]]$preds
  matrix_positions = unlist(simulation$matrix_position_list)
  qtl_list = simulation$samples[[s]]$qtl_list
  qtl_data_frame = data.frame()
  true_rm_fun_list = second_true_rm(qtl_list,simulation$position_list,simulation$matrix_position_list,radius=radius)
  true_rm_fun = Reduce("+",lapply(true_rm_fun_list,unlist))
  for(o in qtl_list) {
    for(q in o) {
      qtl_data_frame = rbind(qtl_data_frame,data.frame(position=simulation$position_list[[q$chromosome]][q$position],size=q$size,genitor=q$genitor))
    }
  }
  qtl_data_frame$genitor =as.factor(qtl_data_frame$genitor)
  par(mfrow=c(max(1,min(2,length(methods)%/%2)),max(1,length(methods)%/%2)))
  for(m in 1:length(methods)) {
    pred_rm_fun = unlist(pred_rm(preds[[m]],simulation$matrix_position_list,radius=radius))
    plot(matrix_positions,pred_rm_fun,type="l",main = paste0(methods[m]),ylab="Variance estimée",xlab="Position génétique (cM)",ylim=c(0,max(max(pred_rm_fun),max(true_rm_fun))))
    lines(matrix_positions,true_rm_fun,col="blue",lty=3)
    points(qtl_data_frame$position,qtl_data_frame$size,col=qtl_data_frame$genitor,pch=19)
    abline(v=simulation$chr_start,col="gray",lty=3)
    legend("topright",legend=levels(qtl_data_frame$genitor),col=1:length(levels(qtl_data_frame$genitor)),pch=1)
  }
}


second_plot_sim = function(simulation,s,methods=c("Horseshoe","Fused","Fusion")) {
  sample = simulation$samples[[s]]

  qtl_data = data.frame()
  for(origin in 1:length(sample$qtl_list)){
    for(qtl in sample$qtl_list[[origin]]){
      result = list(chromosome = qtl$chromosome,position = simulation$position_list[[qtl$chromosome]][qtl$position],size = qtl$size,origin = qtl$genitor)
      qtl_data = rbind(qtl_data,result)
    }
  }
  qtl_data$origin = as.factor(qtl_data$origin)
  par(mfrow=c(length(methods),1))
  for(m in 1:length(methods)){
    pred = sample$result[[m]]$lambda2_median
    plot(unlist(simulation$matrix_position_list),unlist(pred),type="l",col="black",main=paste0(methods[m]),ylab="Variance estimée",xlab="Position génétique (cM)")
    abline(v=simulation$chr_start,col="gray",lty=3)
    points(qtl_data$position,qtl_data$size/3,col=qtl_data$origin,pch=19)
    legend("topright",legend=levels(qtl_data$origin),col=1:length(levels(qtl_data$origin)),pch=1)
  }
}

my_error_boxplot = function(error_comparison,quantiles=c(0.3,0.5,0.7)){
  par(mfrow=c(2,2))
  error_comparison$radius = as.factor(error_comparison$radius)
  error_comparison$method = as.factor(error_comparison$method)
  radii = levels(error_comparison$radius)
  methods = levels(error_comparison$method)
  for(e in c("correlation","misestimation","overestimation","underestimation")){
    plot(rep(0,length(radii)),rep(0,length(radii)),pch=0,ylim=c(0,1),xlim=c(2.5,max(as.numeric(radii))),main=switch(e,
                                                                                                             "correlation"="1 - Corrélation",
                                                                                                             "misestimation"="Mésestimation",
                                                                                                             "overestimation"="Sur-estimation",
                                                                                                             "underestimation"="Sous-estimation"),ylab="Médiane de l'échantillon",xlab="Rayon de la somme mobile(cM)")
    for(m in 1:length(methods)){
      for(q in 1:length(quantiles)){
        temp =rep(0,length(radii))
        for(r in 1:length(radii)){
          temp[r] =quantile(error_comparison[error_comparison$method==methods[m]&error_comparison$radius==radii[r],][[e]],quantiles[q])
        }
        lines(as.numeric(radii),temp,col=m,lty=1 + 3*(q%%2),type="o",pch=19,lwd = 2 - (q%%2))
      }
    }
    legend("topright",legend=methods,col=1:length(methods),pch=19)
  }
}


