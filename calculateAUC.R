#install.packages("MESS")
library(MESS)

rocAUC = function(term_list,score_list,score_list_labels,top_rank){
  
  #term_list is a character vector of terms -- for example the genes from a ChEA or ENCODE experiment
  #score_list is a vector of scores (e.g. pvals, correlation coefficients) with genes as rownames
  #score_list_labels is a character vector of terms corresponding to score_list (e.g. gene names)
  #top_rank=c("least","greatest") -- indicates whether the top rank is min(score_list) or max(score_list), respectively
  
  score_list_group = rep(0, length(score_list))
  overlap_idx = is.element(score_list_labels,term_list)
  n_overlap = sum(overlap_idx)
  score_list_group[overlap_idx] = rep(1,n_overlap)
  if (top_rank == "greatest") score_list_group = score_list_group[order(-score_list)]
  else if(top_rank == "least") score_list = score_list_group[order(score_list)] 
  else{
    print("Error: Enter a valid argument for top_rank (i.e. 'greatest' or 'least')") 
    stop()
  } 
  cum_labels = cumsum(score_list_group)/n_overlap
  
  return(auc((1:length(cum_labels))/length(cum_labels),cum_labels))
  
}

dataset1=creedsFisher
dataset2=coexp500

combos<-combn(union(colnames(dataset1),colnames(dataset2)),2,simplify=F)
combos1<-combn(colnames(dataset1),2,simplify=F)
combos2<-combn(colnames(dataset2),2,simplify=F)
combos<-setdiff(combos,combos1)
combos<-setdiff(combos,combos2)

calculateAUC<-function(d) { 
  
  tf1<-dataset1[[d[1]]]
  tf1<-tf1[!is.na(tf1)]
  
  tf2<-dataset2[[d[2]]]
  tf2<-tf2[!is.na(tf2)]
  
  auc<-rocAUC(term_list = tf1,score_list = tf2,score_list_labels = rownames(dataset2),top_rank = "greatest")
  
  row<-data.frame(d[1],d[2],auc,stringsAsFactors = FALSE)
  colnames(row)=c("tf1","tf2","auc")
  return(row)
}

datasetAUC<-as.data.frame(do.call(rbind,(lapply(combos, function(d) calculateAUC(d)))))
