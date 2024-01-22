# preciscion recall new
suppressWarnings({
  
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(pracma)
  
})
#change input
csd_network<-data.table::fread("data/control_mayo_TCX_wTO.txt")


exit_file<- "wTO_new_t1"
dir_1<-"results/naive/"

#maxd<-max(csd_network$rho1)-0.00000001
cutoffs<-seq(0.0,1, length.out = 192)

filetable<-paste(dir_1,"values_complete_",exit_file,"_cutoffs",length(cutoffs),".txt",sep="")

#########

csd_network <- csd_network %>% select(-(pval:pval.adj))
colnames(csd_network) <- c("Gene1","Gene2","rho1")


naive_gold<-data.table::fread("data/gold_standard/gobp_naive_standard.dat",header = FALSE)


##

merged_df1 <- merge(csd_network, naive_gold, by.x = c("Gene1", "Gene2"), by.y = c("V1", "V2"))
merged_df2 <- merge(csd_network, naive_gold, by.x = c("Gene2", "Gene1"), by.y = c("V1", "V2"))
gene_brain<- rbind(merged_df1,merged_df2)
gene_brain <- gene_brain %>% select(-(rho1))

csd_network <- merge(csd_network, gene_brain, by.x = c("Gene1", "Gene2"), by.y = c("Gene1", "Gene2"), all.x = TRUE) #gives to the network if it is in gold standard


print(head(csd_network))
print(paste(" Number of cutoffs",length(cutoffs)))

tp=0
fp=0
fn=0
tn=0

for (i in 1:length(cutoffs)) {
  
  tp[i] <- sum(abs(csd_network$rho1) > cutoffs[i] & csd_network$V3 == 1, na.rm = TRUE)
  fp[i] <- sum(abs(csd_network$rho1) > cutoffs[i] & csd_network$V3 == 0, na.rm = TRUE)
  fn[i] <- sum(abs(csd_network$rho1) <= cutoffs[i] & csd_network$V3 == 1, na.rm = TRUE)
  tn[i] <- sum(abs(csd_network$rho1) <= cutoffs[i] & csd_network$V3 == 0, na.rm = TRUE)
}



values<-data.frame(tp = tp, fp = fp, fn = fn, tn = tn)
values

fwrite(values,filetable, sep ="\t")

