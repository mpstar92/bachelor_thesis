suppressWarnings({
  
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(pracma)
  
})

## change the input file
csd_results_test<-data.table::fread("results/csd_results_subsampling_spearmen_100_25.txt")


maxd<-max(csd_results_test$rho1)-0.00000001
cutoffs<-seq(0.0,maxd, length.out = 192)

length(cutoffs)

## for output file

exit_file<- "sub_spear"
#directory
dir_1<-"results/sub_spear/"

#outfile, precision recall plot and table of tp,tn,fp,fn

fileout<-paste(dir_1,"outfile_parallel_complete_",exit_file,"_cutoffs",length(cutoffs),".txt",sep="")
fileplot<-paste(dir_1,"plot_3",exit_file,"_cutoffs",length(cutoffs),".png",sep="")
filetable<-paste(dir_1,"values_complete",exit_file,"_cutoffs",length(cutoffs),".txt",sep="")



# cut off the things from the data that is not used
csd_results_test <- csd_results_test %>% select(-(rho2:dVal))

########
brain_gold<-data.table::fread("data/gold_standard/brain_goldst.dat",header = FALSE)
brain_gold<-brain_gold[order(V1)]
#####
#for output
print(exit_file)

print(cat("number of cutoffs to be run in parallel ",length(cutoffs) ))
print("")

############



####parallel getting tp etc values
#### function: compares goldstandard to co-expression network and
#### returns the confusion matrix at a cutoff




getting_values<-function(brain_gold,csd_results_test,cutoff,i){
  brain_gold<-brain_gold
  csd_results_test<-csd_results_test
  cutoff<-cutoff
  
  go_1<-brain_gold$V1[1]
  csd_results_test2<-filter(csd_results_test, Gene1 == go_1 | Gene2 == go_1  )
  
  tp=0
  fp=0
  fn=0
  tn=0
  
  br<-length(brain_gold$V1)
  
  for(i in 1:br){
    if(go_1 == brain_gold$V1[i]){
      csd_results_test3<-filter(csd_results_test2, Gene1 == brain_gold$V2[i] | Gene2 == brain_gold$V2[i])
      if(is.na(csd_results_test3$Gene1[1]) == FALSE ) {
        #calculate tp etc
        if(abs(csd_results_test3$rho1[1]) >  cutoff & brain_gold$V3[i] == 1) tp=tp+1
        if(abs(csd_results_test3$rho1[1]) >  cutoff & brain_gold$V3[i] == 0) fp=fp+1
        if(abs(csd_results_test3$rho1[1]) <= cutoff & brain_gold$V3[i] == 1) fn=fn+1
        if(abs(csd_results_test3$rho1[1]) <= cutoff & brain_gold$V3[i] == 0) tn=tn+1
      }
    }
    else{
      go_1 <- brain_gold$V1[i]
      csd_results_test2<-filter(csd_results_test, Gene1 == go_1 | Gene2 == go_1  )
      csd_results_test3<-filter(csd_results_test2, Gene1 == brain_gold$V2[i] | Gene2 == brain_gold$V2[i])
      if(is.na(csd_results_test3$Gene1[1]) == FALSE ) {
        #calculate tp etc
        if(abs(csd_results_test3$rho1[1]) >  cutoff & brain_gold$V3[i] == 1) tp=tp+1
        if(abs(csd_results_test3$rho1[1]) >  cutoff & brain_gold$V3[i] == 0) fp=fp+1
        if(abs(csd_results_test3$rho1[1]) <= cutoff & brain_gold$V3[i] == 1) fn=fn+1
        if(abs(csd_results_test3$rho1[1]) <= cutoff & brain_gold$V3[i] == 0) tn=tn+1
      }
    }
    
  }
  u<-c(tp,fp,fn,tn)
}
##

#parallel use maximum cores
cl <- parallel::makeCluster(detectCores(),outfile=fileout)

# Activate cluster for foreach library
doParallel::registerDoParallel(cl)

# call function in parallel
out<-foreach(i=1:length(cutoffs),.packages = "dplyr" )%dopar%{
  
  #cutoff_1<-cutoffs[i]
  
  getting_values(brain_gold,csd_results_test,cutoffs[i],i)
}

stopCluster(cl)
#stop parallel


## add the confusion matrices to one table
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



values<-tibble(tp_s,fp_s,fn_s,tn_s)

####outfile the tp etc values 
fwrite(values,filetable, sep ="\t")

