​
library(tidyverse)
​
setwd("~/Dropbox (Sanjana Lab)/_Neuron reporters HR/__TF screen manuscript/Revision 18Jul2022/New analysis/Network analysis/NewAnalysis_9_5_22/GeneBodyExpanded")
## install.packages('R.utils')
​
## df <- data.table::fread("GeneBody_expanded.tsv.gz")
## write.csv(df, 'GeneBody_expanded.csv') 
​
## df <- data.table::fread("IntergenicAndGeneBody_expanded.tsv.gz")
## write.csv(df, 'IntergenicAndGeneBody_expanded.tsv.csv') 
​
##### PART I, build networks ###########################################################################
#####################################################################################################
######################## Prepare source file for analysis (from Chao 08/25/22)
​
​
df <- read.csv('GeneBody_expanded.csv') ## 8 more columns compared to "IntergenicAndGeneBody_expanded.csv"
​
df <- df[c(2:19, 28:43)]
​
df$rATAC_16H <- df$atac.cnt.H16/df$atac.cnt.ES
df$rATAC_D1 <- df$atac.cnt.D1/df$atac.cnt.ES
df$rATAC_D4 <- df$atac.cnt.D4/df$atac.cnt.ES
​
write.csv(df, 'df_GeneBody_expanded.csv') 
​
​
########################################################################################
##### PART I, Step 1, built network at 15H #############################################
########################################################################################
​
# load file for analysis
​
df <- read.csv('df_GeneBody_expanded.csv')
​
df <- df[c(2:38)]
​
############################################ Look for NEUROG2/1's target at 15h 
​
NEUROG21 <- filter(df, regulator_gene %in% c("NEUROG1","NEUROG2"))
​
Network_15H <- subset(NEUROG21, abs(NEUROG21$rna.tar.l2fc.H16) > 1 ## target gene expression > 2 fold increase
                       & NEUROG21$rna.tar.fdr.H16 < 0.05      ###### try different cutoff here
                       & NEUROG21$atac.cnt.H16 > 50   ###### try different cutoff here
                       & NEUROG21$rATAC_16H >1
                       & NEUROG21$rna.tar.cnt.H16 > 10)
​
distinct(Network_15H, target_gene) %>% count() #46
​
write_csv(Network_15H, 'Network_15H.csv')  ######## save network at 15H
​
​
########################################################################################
##### PART I, Step 2, built network at D1 #############################################
########################################################################################
​
# load file for analysis
df <- read.csv('df_GeneBody_expanded.csv')
df <- df[c(2:38)]
Network_15H <- read.csv('Network_15H.csv')
​
############################################ Find regulators from 15H network to build D1 network
######## NEUROG2/1 plus the up-regulated targets in 15H network serve as regulators in D1 network
​
df_sub <- subset(Network_15H, Network_15H$rna.tar.l2fc.H16 >1 ) 
list_reg <- df_sub[c(2)]
list_reg <- distinct(list_reg, target_gene)
​
names(list_reg)[1] <- "regulator_gene"
​
##add NEUROG1 and NEUROG2 into the list
reg_NGNs <- data.frame (regulator_gene = c("NEUROG1", "NEUROG2"))
​
list_reg_D1 <- rbind(list_reg, reg_NGNs)
distinct(list_reg_D1, regulator_gene) %>% count()
​
## 41 TFs plus NEUROG2 and NEUROG1 as regulators for D1 network
​
############################################ build D1 network
​
dfDay1 <- subset(df, abs(df$rna.tar.l2fc.D1) > 1 
               & df$rna.tar.fdr.D1 < 0.05
               & df$atac.cnt.D1 > 50
               & df$rATAC_D1 >1
               & df$rna.tar.cnt.D1 > 10)
​
xy = merge(list_reg_D1, dfDay1, by ='regulator_gene', all.x = TRUE)
xy <- na.omit(xy) ## remove rows with NA
​
write.csv(xy, 'Network_D1.csv')  ######## save network at D1
​
​
########################################################################################
##### PART I, Step 3, built network at D4 #############################################
########################################################################################
​
# load file for analysis
df <- read.csv('df_GeneBody_expanded.csv')
df <- df[c(2:38)]
Network_D1 <- read.csv('Network_D1.csv')
Network_D1 <- Network_D1[c(2:38)]
​
############################################ Get regulators from D1 network to build D4 network
#### NEUROG2/1 plus the up-regulated targets in D1 network serve as regulators in D4 network####
​
df_sub <- subset(Network_D1, Network_D1$rna.tar.l2fc.D1 >1 ) 
list_reg <- df_sub[c(2)]
list_reg <- distinct(list_reg, target_gene)
names(list_reg)[1] <- "regulator_gene"  
​
##add NEUROG1 and NEUROG2 into the list
reg_NGNs <- data.frame (regulator_gene = c("NEUROG1", "NEUROG2"))
​
list_reg_D4 <- rbind(list_reg, reg_NGNs)
list_reg_D4 <- distinct(list_reg_D4, regulator_gene)
​
## 75 TFs plus NEUROG2 and NEUROG1 as regulators for D4 network
​
############################################ build D4 network
​
dfDay4 <- subset(df, abs(df$rna.tar.l2fc.D4) > 1 
                 & df$rna.tar.fdr.D4 < 0.05
                 & df$atac.cnt.D4 > 50
                 & df$rATAC_D4 >1
                 & df$rna.tar.cnt.D4 > 10)
​
xy = merge(list_reg_D4, dfDay4, by ='regulator_gene', all.x = TRUE)
xy <- na.omit(xy) ## remove rows with NA
​
write.csv(xy, 'Network_D4.csv')  ######## save network at D4
​
​
######## PART II . network analysis ##########################################################################
##############################################################################################################
###########################################
​
## load files for analysis:
df15H <- read.csv("Network_15H.csv")
df15H <- df15H[c(1,2,25,28,33,34)]
dfD1 <- read.csv("Network_D1.csv")
dfD1 <- dfD1[c(2,3,27,30,34,35)]
dfD4 <- read.csv("Network_D4.csv")
dfD4 <- dfD4[c(2,3,28,31,34,35)]
​
​
## to get the number of regulators and targets
df15H_reg <- distinct(df15H, regulator_gene) ## 2
df15H_tar <- distinct(df15H, target_gene)    ## 46
​
dfD1_reg <- distinct(dfD1, regulator_gene) ## 17
dfD1_tar <- distinct(dfD1, target_gene)    ## 112
​
dfD4_reg <- distinct(dfD4, regulator_gene) ## 12
dfD4_tar <- distinct(dfD4, target_gene)   ##  327
​
## to get the number of activation or repression
​
df15H_act <- subset(df15H, df15H$rna.tar.l2fc.H16 >0)  ## 48
df15H_rep <- subset(df15H, df15H$rna.tar.l2fc.H16 <0)  ## 7
​
dfD1_act <- subset(dfD1, dfD1$rna.tar.l2fc.D1 >0) ## 279
dfD1_rep <- subset(dfD1, dfD1$rna.tar.l2fc.D1 <0)    ## 134
​
dfD4_act <- subset(dfD4, dfD4$rna.tar.l2fc.D4 >0) ## 439
dfD4_rep <- subset(dfD4, dfD4$rna.tar.l2fc.D4 <0)   ##  497
