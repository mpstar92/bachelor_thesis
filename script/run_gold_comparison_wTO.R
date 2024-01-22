suppressWarnings({
  
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(PRROC)
  library(pracma)
  
})



## change the input file
csd_results_test<-data.table::fread("data/control_mayo_TCX_wTO.txt")

###naming the output files

# number of cutoff for the precision recall curve
maxd<-max(csd_results_test$wTO)-0.00000001
cutoffs<-seq(0.0,maxd, length.out = 192)

length(cutoffs)
48*10

exit_file<- "wTO_control"
#directory
dir_1<-"results/wTO/"

fileout<-paste(dir_1,"outfile_parallel_",exit_file,"_cutoffs",length(cutoffs),".txt",sep="")
fileplot<-paste(dir_1,"plot_",exit_file,"_cutoffs",length(cutoffs),".png",sep="")
filetable<-paste(dir_1,"values_",exit_file,"_cutoffs",length(cutoffs),".txt",sep="")




csd_results_test <- csd_results_test %>% select(-(pval:pval.adj))

########
brain_gold<-data.table::fread("data/gold_standard/brain_goldst.dat",header = FALSE)
#####

print(head(csd_results_test))
print(head(brain_gold))

print(cat("length of co expression network",length(csd_results_test$Node.1)))
print(cat("length of Gold Standard Brain",length(brain_gold$V1)))

print(cat("number of cutoffs to be run in parallel ",length(cutoffs) ))
print("")

############



####parallel getting tp etc values
######


getting_values<-function(brain_gold,csd_results_test,cutoff,i){
  brain_gold<-brain_gold
  csd_results_test<-csd_results_test
  cutoff<-cutoff
  
  go_1<-brain_gold$V1[1]
  csd_results_test2<-filter(csd_results_test, Node.1 == go_1 | Node.2 == go_1  )
  tp=0
  fp=0
  fn=0
  tn=0
  
  br<-length(brain_gold$V1)
  
  for(i in 1:br){
    if(go_1 == brain_gold$V1[i]){
      csd_results_test3<-filter(csd_results_test2, Node.1 == brain_gold$V2[i] | Node.2 == brain_gold$V2[i])
      if(is.na(csd_results_test3$Node.1[1]) == FALSE ) {
        #calculate tp etc
        if(abs(csd_results_test3$wTO[1]) >  cutoff & brain_gold$V3[i] == 1) tp=tp+1
        if(abs(csd_results_test3$wTO[1]) >  cutoff & brain_gold$V3[i] == 0) fp=fp+1
        if(abs(csd_results_test3$wTO[1]) <= cutoff & brain_gold$V3[i] == 1) fn=fn+1
        if(abs(csd_results_test3$wTO[1]) <= cutoff & brain_gold$V3[i] == 0) tn=tn+1
      }
    }
    else{
      go_1 <- brain_gold$V1[i]
      csd_results_test2<-filter(csd_results_test, Node.1 == go_1 | Node.2 == go_1  )
      csd_results_test3<-filter(csd_results_test2, Node.1 == brain_gold$V2[i] | Node.2 == brain_gold$V2[i])
      
      if(is.na(csd_results_test3$Node.1[1]) == FALSE ) {
        #calculate tp etc
        if(abs(csd_results_test3$wTO[1]) >  cutoff & brain_gold$V3[i] == 1) tp=tp+1
        if(abs(csd_results_test3$wTO[1]) >  cutoff & brain_gold$V3[i] == 0) fp=fp+1
        if(abs(csd_results_test3$wTO[1]) <= cutoff & brain_gold$V3[i] == 1) fn=fn+1
        if(abs(csd_results_test3$wTO[1]) <= cutoff & brain_gold$V3[i] == 0) tn=tn+1
      }
    }
    
  }
  u<-c(tp,fp,fn,tn)
}
##

#parallel
cl <- parallel::makeCluster(detectCores(),outfile=fileout)

print(cl)
# Activate cluster for foreach library
doParallel::registerDoParallel(cl)


out<-foreach(i=1:length(cutoffs),.packages = "dplyr" )%dopar%{
  
  #cutoff_1<-cutoffs[i]
  
  getting_values(brain_gold,csd_results_test,cutoffs[i],i)
}

stopCluster(cl)
#stop parallel

tp_s<-c()
fp_s<-c()
fn_s<-c()
tn_s<-c()
for(i in 1: length(out)){
  tp_s<-append(tp_s,out[[i]][1])
  fp_s<-append(fp_s,out[[i]][2])
  fn_s<-append(fn_s,out[[i]][3])
  tn_s<-append(tn_s,out[[i]][4])
}

tp_s
fp_s
fn_s
tn_s

values<-tibble(tp_s,fp_s,fn_s,tn_s)
values

precision = tp_s/ (tp_s +fp_s )


area_under_curve <- trapz(recall,precision)
print(area_under_curve)

####outfile the tp etc values and the precision recall curve
fwrite(values,filetable, sep ="\t")
