

## check for linearity

# Example data
set.seed(42)  # Setting seed for reproducibility


gene1<-"ENSG00000233008"
gene2<-"ENSG00000238164"



x<- Matrix_mayo_TCX_ctrl[,gene1]
y<- Matrix_mayo_TCX_ctrl[,gene2]
# Create a scatterplot
plot(x,y, main =paste("Scatterplot of ", gene1, " and ", gene2),xlab = gene1, ylab = gene2)
abline(lm(y ~ x), col = "red")



# check for normalization

gene<-"ENSG00000238009"

sample_data <-Matrix_mayo_TCX_ctrl[,gene]


# Create a normal probability plot (Q-Q plot)
qqnorm(sample_data, main = paste("Q-Q Plot of Gene", gene))
qqline(sample_data, col = "red")



####################

# Create a histogram
hist(sample_data, main = paste("Histogram of Gene ", gene), xlab = "Expression Value", col = "lightblue", border = "black")



#####################

#plots for subsamlping sets





library(ggplot2)

mse_results<-read.delim("results/mse/mse_results.txt")

mse_results_sub_10<-mse_results[mse_results$n_sub == "10", ]
 


mse_results_10<-read.delim("results/mse/mse_results_sub10.txt")

mse_results_10 <- rbind(mse_results_10 , mse_results_sub_10)


plot(mse_results_10$n_it,mse_results_10$mse_1,ylim = c(0.05, 0.2))

############
mse_results_20<-read.delim("results/mse/mse_results_sub20.txt")
mse_results_sub_20<-mse_results[mse_results$n_sub == "20", ]

mse_results_20 <- rbind(mse_results_20 , mse_results_sub_20)

plot(mse_results_20$n_it,mse_results_20$mse_1,ylim = c(0.05, 0.2))

############################

mse_results_30<-read.delim("results/mse/mse_results_sub30.txt")
mse_results_sub_30<-mse_results[mse_results$n_sub == "30", ]
 
mse_results_30 <- rbind(mse_results_30 , mse_results_sub_30)

plot(mse_results_30$n_it,mse_results_30$mse_1,ylim = c(0.05, 0.2))


###############


mse_results_all <- rbind(mse_results_30 , mse_results_20 , mse_results_10)
####


mse_results_it_10<-mse_results[mse_results$n_it == "10", ]
mse_results_it_20<-mse_results[mse_results$n_it == "20", ]
mse_results_it_30<-mse_results[mse_results$n_it == "30", ]
mse_results_it_40<-mse_results[mse_results$n_it == "40", ]
mse_results_it_50<-mse_results[mse_results$n_it == "50", ]
mse_results_it_100<-mse_results[mse_results$n_it == "100", ]




plot(mse_results_it_20$n_sub,mse_results_it_20$first,ylim = c(0, 0.3))

ggplot(mse_results_all, aes(x = n_it, y = mse_1, color = factor(n_sub))) +
  geom_point() +
  scale_y_continuous(limits = c(0.07, 0.13)) +
  labs(x = "n_it", y = "first", color = "n_sub") +
  ggtitle("Scatter Plot of 'first' vs 'n_it' with colored points by 'n_sub'")



ggplot(mse_results_all, aes(x = n_it, y = mse_1, color = factor(n_sub))) +
  geom_point() +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(x = "n_it", y = "first", color = "n_sub") +
  ggtitle("Scatter Plot of 'first' vs 'n_it' with colored points by 'n_sub'")


ggplot(mse_results_all, aes(x = n_sub, y = mse_1, color = factor(n_it))) +
  geom_point() +
  scale_y_continuous(limits = c(0.07, 0.13)) +
  labs(x = "n_sub", y = "first", color = "n_sub") +
  ggtitle("Scatter Plot of 'first' vs 'n_sub' with colored points by 'n_it'")

ggplot(mse_results_all, aes(x = n_sub, y = mse_1, color = factor(n_it))) +
  geom_point() +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(x = "n_sub", y = "first", color = "n_sub") +
  ggtitle("Scatter Plot of 'first' vs 'n_sub' with colored points by 'n_it'")


#########################
library(pracma)
library(data.table)











values_wTO<-fread("results/wTO/values_wTO_control_cutoffs192.txt")

#### wTO
tp_w<-values_wTO$tp_s
fp_w<-values_wTO$fp_s
fn_w<-values_wTO$fn_s
tn_w<-values_wTO$tn_s

values_wTO[, precision := ifelse((tp_w + fp_w) == 0, 1, tp_w / (tp_w + fn_w))]
values_wTO[, recall := ifelse((tp_w + fn_w) == 0, 1, tp_w / (tp_w + fn_w))] 

values_wTO[, true_positive_rate :=  tp_w / (tp_w + fn_w)] 
values_wTO[, false_positive_rate :=  fp_w / (fp_w + tn_w)] 


plot(values_wTO$false_positive_rate,values_wTO$true_positive_rate ,type = "l",col = "blue")

area_under_curve_roc <- abs(trapz(values_wTO$false_positive_rate,values_wTO$true_positive_rate))
area_under_curve_roc


naive_s_p<-fread("results/naive/values_naive_sub_pear_cutoffs79.txt")
tp_n_s_p<-naive_s_p$tp_s
fp_n_s_p<-naive_s_p$fp_s
fn_n_s_p<-naive_s_p$fn_s
tn_n_s_p<-naive_s_p$tn_s
naive_s_p[, precision := ifelse((tp_n_s_p + fp_n_s_p) == 0, 1, tp_n_s_p / (tp_n_s_p + fp_n_s_p))]
naive_s_p[, recall := ifelse((tp_n_s_p + fn_n_s_p) == 0, 1, tp_n_s_p / (tp_n_s_p + fn_n_s_p))] 
naive_s_p[, true_positive_rate :=  tp_n_s_p / (tp_n_s_p + fn_n_s_p)] 
naive_s_p[, false_positive_rate :=  fp_n_s_p / (fp_n_s_p + tn_n_s_p)] 


plot(naive_s_p$recall,naive_s_p$precision,type = "l",col = "blue")
lines(x,y)

area_under_curve_n_sp <- abs(trapz(naive_s_p$recall,naive_s_p$precision))
print(cat("AUC naive sub pear  is ",abs(area_under_curve_n_sp)))


plot(naive_s_p$false_positive_rate,naive_s_p$true_positive_rate ,type = "l",col = "blue")

area_under_curve_roc_naive <- abs(trapz(naive_s_p$false_positive_rate,naive_s_p$true_positive_rate))
area_under_curve_roc_naive



##################

train<-c(5,6,7,5,4,8,7,5)
test<-c(4,5,6,7,8,9,10,9)

mean(train)
mean(test)
sum((test-mean(train))^2)/length(train)

mse<-(sum((test-mean(train))^2))/length(train)
mse

sum((test-mean(train))^2)/length(train)


test-mean(train)
(test-mean(train))^2
sum((test-mean(train))^2)

(mean(test)-mean(train))^2


mse2<-(sum((test-train)^2))/length(train)
mse2

mse3<-((sum(test)-(length(train)*mean(train)))^2)/length(train)
mse3
mse4<-((mean(test)-mean(train))^2)
mse4

x<-(5-3.5)^2+(6-3.5)^2
x
y<-((5-3)^2+(6-3)^2)-2*(0.5^2)
y
4-2.25
1.75^2



subset_matrix_ctrl <- Matrix_mayo_TCX_ctrl[1:42, 1:1000]

n_res <- ncol(subset_matrix_ctrl)
gene_names <- colnames(subset_matrix_ctrl)


res_matrix2 <- array(0, c(n_res, n_res),
                    dimnames = list(gene_names, gene_names)
)
dim(res_matrix)


dim(Matrix_mayo_TCX_ctrl)
dim(subset_matrix_ctrl)


n <- 1  # Set the desired length of the vector

# Creating the matrix with vectors of length n with all zeroes
res_matrix <- array(
  data = rep(list(rep(0, n)), n_res^2),
  dim = c(n_res, n_res),
  dimnames = list(gene_names, gene_names)
)

res_matrix[[1,1]][2]

#####################

subset_matrix_ctrl <- Matrix_mayo_TCX_ctrl[1:42, 1:10000]
gene_names <- colnames(subset_matrix_ctrl)
n_res <- length(gene_names)

# Set the third dimension size
n_dim3 <- 3 #n_it

# Create a 3D array
res_matrix_3d <- array(
  data = rep(0, n_res^2 * n_dim3),  # Initialize with vectors of length n_res^2 filled with zeros
  dim = c(n_res, n_res, n_dim3),
  dimnames = list(gene_names, gene_names, paste0("Dim3_", 1:n_dim3))
)



res_matrix <- array(0, c(n_res, n_res),
                    dimnames = list(gene_names, gene_names)
)
for(i in 1:n_res){
  for(j in 1:n_res){
    res_matrix[i,j]<-sum((res_matrix_3d[i,j,]-rho[i,j])^2)
  }
}
res_matrix<-res_matrix/n_dim3



rm(res_matrix_3d)




subset_matrix_ctrl <- Matrix_mayo_TCX_ctrl[1:3, 1:1000]
subset_matrix_ad <- Matrix_mayo_TCX_ad[1:2, 1:1000]


dim(subset_matrix_ctrl)
dim(Matrix_mayo_TCX_ctrl)



library(pracma)
library(data.table)



gold_brain_sp_100_20_tp<-read.delim("results/sub_pear/values_final_sub_pearson_100_20_cutoffs192_tp.txt")
gold_brain_sp_100_20_fp<-read.delim("results/sub_pear/values_final_sub_pearson_100_20_cutoffs192_fp.txt")
gold_brain_sp_100_20_fn<-read.delim("results/sub_pear/values_final_sub_pearson_100_20_cutoffs192_fn.txt")
gold_brain_sp_100_20_tn<-read.delim("results/sub_pear/values_final_sub_pearson_100_20_cutoffs192_tn.txt")


gold_brain_sp_100_20<-data.table(gold_brain_sp_100_20_tp,gold_brain_sp_100_20_fp,gold_brain_sp_100_20_fn,gold_brain_sp_100_20_tn)

fwrite(gold_brain_sp_100_20,"results/sub_pear/values_final_sub_pearson_100_20_cutoffs192.txt", sep ="\t")


gold_brain_ss_100_20_tp<-read.delim("results/sub_spear/values_final_sub_spearmen_100_20_cutoffs192_tp.txt")
gold_brain_ss_100_20_fp<-read.delim("results/sub_spear/values_final_sub_spearmen_100_20_cutoffs192_fp.txt")
gold_brain_ss_100_20_fn<-read.delim("results/sub_spear/values_final_sub_spearmen_100_20_cutoffs192_fn.txt")
gold_brain_ss_100_20_tn<-read.delim("results/sub_spear/values_final_sub_spearmen_100_20_cutoffs192_tn.txt")


gold_brain_ss_100_20<-data.table(gold_brain_ss_100_20_tp,gold_brain_ss_100_20_fp,gold_brain_ss_100_20_fn,gold_brain_ss_100_20_tn)

fwrite(gold_brain_ss_100_20,"results/sub_spear/values_final_sub_spearmen_100_20_cutoffs192.txt", sep ="\t")


values[, precision := ifelse((tp_n + fp_n) == 0, 1, tp_n / (tp_n + fp_n))]
values[, recall := ifelse((tp_n + fn_n) == 0, 1, tp_n / (tp_n + fn_n))]

values[, tp_rate =  tp_n / (tp_n + fn_n)]
values[, fp_rate =  fp_n / (fp_n + tn_n)] 



tp_n<-gold_brain_ss_100_20$tp_s
fn_n<-gold_brain_ss_100_20$fn_s

gold_brain_ss_100_20[, tp_rate :=  tp_n / (tp_n + fn_n)]

baseline<-1/2
x<-seq(0,1,by=0.01  )
y<-rep(baseline*x,times =length(x))



plot(x,y)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "X-axis", ylab = "Y-axis")

# Add a line from (0,0) to (1,1)
lines(c(0, 1), c(0, 1), col = "blue")




cutoffs <-seq(0.01,0.99, by=0.00513)


max_c<-0.6
maxd<-max(csd_results3$rho1)-0.00000001

seq1<-seq(0.0,maxd, length.out = 192)
seq2<-seq(0, max_c, length.out = 192)
length(seq1)
length(seq2)

print(paste("ddd",max_c))


1783 +173708
length(brain_gold$V1)
