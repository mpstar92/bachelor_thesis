# friedmann test/anova

suppressWarnings({
  
  library(data.table)
  library(magrittr)
  library(dplyr)
  
  
  
})

brain_gold<-data.table::fread("data/gold_standard/brain_goldst.dat",header = FALSE)

#first network
csd_s_s_25<-data.table::fread("results/csd_results_subsampling_spearmen_100_25.txt")
csd_s_s_25 <- csd_s_s_25 %>% select(-(rho2:dVal))


##


###
merged_df1 <- merge(csd_s_s_25, brain_gold, by.x = c("Gene1", "Gene2"), by.y = c("V1", "V2"))
merged_df2 <- merge(csd_s_s_25, brain_gold, by.x = c("Gene2", "Gene1"), by.y = c("V1", "V2"))
gene_brain<- rbind(merged_df1,merged_df2)
gene_brain <- gene_brain %>% select(-(rho1))
####

csd_s_s_25$version <- 1   #version of subsampling as int
csd_s_s_25$ID <- 1:nrow(csd_s_s_25) #ID for every genepair
csd_s_s_25$rank_rho <- rank(abs(csd_s_s_25$rho1)) #ranks for all correlations for better comparison with wto

csd_s_s_25 <- merge(csd_s_s_25, gene_brain, by.x = c("Gene1", "Gene2"), by.y = c("Gene1", "Gene2"), all.x = TRUE) #gives to the network if it is in gold standard
#######

table_prep<-function(csd_network,gene_brain, version){
  
  csd_network$version <- version   #version of subsampling as int
  csd_network$ID <- 1:nrow(csd_network) #ID for every genepair
  csd_network$rank_rho <- rank(abs(csd_network$rho1)) #ranks for all correlations for better comparison with wto
  csd_network <- merge(csd_network, gene_brain, by.x = c("Gene1", "Gene2"), by.y = c("Gene1", "Gene2"), all.x = TRUE) #gives to the network if it is in gold standard
  #csd_network <- csd_network %>% select(-(Gene1:Gene2))
}

#all other networks

csd_s_p_25<-data.table::fread("results/csd_results_subsampling_pearson_100_25.txt")
csd_s_p_25 <- csd_s_p_25 %>% select(-(rho2:dVal))
csd_s_p_25<-table_prep(csd_s_p_25,gene_brain,version = 2)
all_networks <- rbind(csd_s_s_25,csd_s_p_25)
rm(csd_s_s_25,csd_s_p_25)

csd_b_p<-data.table::fread("results/csd_results_bootstrap_pearson_1000_final.txt")
csd_b_p <- csd_b_p %>% select(-(rho2:dVal))
csd_b_p<-table_prep(csd_b_p,gene_brain,version = 3)
all_networks <- rbind(all_networks,csd_b_p)
rm(csd_b_p)

csd_b_s<-data.table::fread("results/csd_results_bootstrap_spearmen_1000_final.txt")
csd_b_s <- csd_b_s %>% select(-(rho2:dVal))
csd_b_s<-table_prep(csd_b_s,gene_brain,version = 4)
all_networks <- rbind(all_networks,csd_b_s)
rm(csd_b_s)

csd_wTO<-data.table::fread("data/control_mayo_TCX_wTO.txt")
csd_wTO <- csd_wTO %>% select(-(pval:pval.adj))
colnames(csd_wTO) <- c("Gene1","Gene2","rho1")
csd_wTO<-table_prep(csd_wTO,gene_brain,version = 5)
all_networks <- rbind(all_networks,csd_wTO)
rm(csd_wTO)
###############


print(head(all_networks[order(-all_networks$rho1), ],30))
all_networks_1 <- all_networks[all_networks$V3 == 1, ]
file1<- paste0("all_networks_go_1.txt")
fwrite(all_networks_1,file1, sep ="\t")

all_networks_0 <- all_networks[all_networks$V3 == 0, ]

file1<- paste0("all_networks_go_0.txt")
fwrite(all_networks_0,file1, sep ="\t")

print(head(all_networks_1[order(-all_networks_1$rho1), ],30))

#############

#all_networks %>%
# group_by(version) %>%
#  summarize(Mdn = median(rank_rho1),
#           Q1 = quantile(rank_rho1, probs = .25),
#          Q3 = quantile(rank_rho1, probs = .75)) %>%
#  as.data.frame()

#all_networks %>%
# group_by(version) %>%
#summarize(Mdn = median(abs(rho1)),
#         Q1 = quantile(abs(rho1), probs = .25),
#        Q3 = quantile(abs(rho1), probs = .75)) %>%
#as.data.frame()
#############
print("Version 1 = Subsampling Spearman 25")
print("Version 2 = Subsampling Pearson 25 ")
print("Version 3 = Bootstrap Spearman")
print("Version 4 = Bootstrap Pearson")
print("Version 5 = wTO")

#####
print("test on actual correlation scores")

friedman.test(y=all_networks_1$rho1, groups=all_networks_1$version, blocks=all_networks_1$ID)

### pairwise test
pairwise.wilcox.test(all_networks_1$rho1, all_networks_1$version, p.adj = "bonf")


print("ranks test")
######
friedman.test(y=all_networks_1$rank_rho, groups=all_networks_1$version, blocks=all_networks_1$ID)


### pairwise test ranks
pairwise.wilcox.test(all_networks_1$rank_rho, all_networks_1$version, p.adj = "bonf")

