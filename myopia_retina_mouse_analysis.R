library(readxl)
library(MetaboAnalystR)
library(dplyr)
library(data.table)
library(EnhancedVolcano)

# 
# # Set 1
# abund_s1 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S1_abundance_ratios"))                                                                            
# abund_log_s1 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S1_abundance_log2"))                                                                            
# abund_group_s1 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S1_abundance_grouped"))                                                                            
# 
# # Set 2
# abund_s2 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S2_abundance_ratios"))                                                                            
# abund_log_s2 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S2_abundance_log2"))                                                                            
# abund_group_s2 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S2_abundance_grouped"))                                                                            

############# Metaboanalyst ###############

# initializing object
mSet <- InitDataObjects("pktable", "stat", FALSE)
# loading in data
mSet <- Read.TextData(mSet, "abundance_S1_grouped_input.csv", "colu", "disc");
# data check
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet);

# filtering features
mSet <- FilterVariable(mSet, "none", "F", 25)

# median normalization, log10 transformation, pareto scaling
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "MedianNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet <- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

# plot volcano plot
mSet <- Volcano.Anal(mSet, FALSE, 1.5, 0, F, 0.05, TRUE, "raw")
mSet <- PlotVolcano(mSet, "volcano_abundance_S1",1, 0, "png", 72, width=NA)

# save transformed dataset
mSet <- SaveTransformedData(mSet)

############## 
# test <- fread("data_normalized_S1.csv")

# reading in csv of 3 groups

grouped_s1 <- fread("grouped_1.csv",sep=',')
grouped_s2 <- fread("grouped_2.csv",sep=',')
grouped_s3 <- fread("grouped_3.csv",sep=',')

grouped_combined <- left_join(grouped_s3,grouped_s2,by = 'Accession')
grouped_combined_2 <- left_join(grouped_combined,grouped_s1,by = 'Accession')
grouped_combined_2_no_na <- na.omit(grouped_combined_2)

# median_LI_1hr <- median(grouped_combined_2_no_na$S1_LI_1hr,grouped_combined_2_no_na$S2_LI_1hr,grouped_combined_2_no_na$S3_LI_1hr,na.rm = T)

# creating row medians from 3 sets for LI group
grouped_combined_2_no_na$median_LI_0hr <- apply(grouped_combined_2_no_na[,c('S1_LI_0hr','S2_LI_0hr','S3_LI_0hr')], 1, median)
grouped_combined_2_no_na$median_LI_1hr <- apply(grouped_combined_2_no_na[,c('S1_LI_1hr','S2_LI_1hr','S3_LI_1hr')], 1, median)
grouped_combined_2_no_na$median_LI_6hr <- apply(grouped_combined_2_no_na[,c('S1_LI_6hr','S2_LI_6hr','S3_LI_6hr')], 1, median)
grouped_combined_2_no_na$median_LI_9hr <- apply(grouped_combined_2_no_na[,c('S1_LI_9hr','S2_LI_9hr','S3_LI_9hr')], 1, median)
grouped_combined_2_no_na$median_LI_D1 <- apply(grouped_combined_2_no_na[,c('S1_LI_D1','S2_LI_D1','S3_LI_D1')], 1, median)
grouped_combined_2_no_na$median_LI_D14 <- apply(grouped_combined_2_no_na[,c('S1_LI_D14','S2_LI_D14','S3_LI_D14')], 1, median)
grouped_combined_2_no_na$median_LI_D3 <- apply(grouped_combined_2_no_na[,c('S1_LI_D3','S2_LI_D3','S3_LI_D3')], 1, median)
grouped_combined_2_no_na$median_LI_D7 <- apply(grouped_combined_2_no_na[,c('S1_LI_D7','S2_LI_D7','S3_LI_D7')], 1, median)

# creating row medians from 3 sets for NL group
grouped_combined_2_no_na$median_NL_0hr <- apply(grouped_combined_2_no_na[,c('S1_NL_0hr','S2_NL_0hr','S3_NL_0hr')], 1, median)
grouped_combined_2_no_na$median_NL_1hr <- apply(grouped_combined_2_no_na[,c('S1_NL_1hr','S2_NL_1hr','S3_NL_1hr')], 1, median)
grouped_combined_2_no_na$median_NL_6hr <- apply(grouped_combined_2_no_na[,c('S1_NL_6hr','S2_NL_6hr','S3_NL_6hr')], 1, median)
grouped_combined_2_no_na$median_NL_9hr <- apply(grouped_combined_2_no_na[,c('S1_NL_9hr','S2_NL_9hr','S3_NL_9hr')], 1, median)
grouped_combined_2_no_na$median_NL_D1 <- apply(grouped_combined_2_no_na[,c('S1_NL_D1','S2_NL_D1','S3_NL_D1')], 1, median)
grouped_combined_2_no_na$median_NL_D14 <- apply(grouped_combined_2_no_na[,c('S1_NL_D14','S2_NL_D14','S3_NL_D14')], 1, median)
grouped_combined_2_no_na$median_NL_D3 <- apply(grouped_combined_2_no_na[,c('S1_NL_D3','S2_NL_D3','S3_NL_D3')], 1, median)
grouped_combined_2_no_na$median_NL_D7 <- apply(grouped_combined_2_no_na[,c('S1_NL_D7','S2_NL_D7','S3_NL_D7')], 1, median)

# combine LI and NL medians into 1 dataframe
grouped_median = cbind(grouped_combined_2_no_na$Accession,grouped_combined_2_no_na[,53:68])

# convert from character to numeric
grouped_median <- grouped_median %>%
  rename('Accession' = 1)

# calculate fold change
FC_0hr <- c(grouped_median$median_LI_0hr/grouped_median$median_NL_0hr)
FC_1hr <- grouped_median$median_LI_1hr/grouped_median$median_NL_1hr
FC_6hr <- grouped_median$median_LI_6hr/grouped_median$median_NL_6hr
FC_9hr <- grouped_median$median_LI_9hr/grouped_median$median_NL_9hr
FC_D1 <- grouped_median$median_LI_D1/grouped_median$median_NL_D1
FC_D14 <- grouped_median$median_LI_D14/grouped_median$median_NL_D14
FC_D3 <- grouped_median$median_LI_D3/grouped_median$median_NL_D3
FC_D7 <- grouped_median$median_LI_D7/grouped_median$median_NL_D7

# test <- c()
# 
# for (i in 1:nrow(grouped_median)) {
#   
#  test[i] <- t.test(grouped_median$median_LI_1hr[i], grouped_median$median_NL_1hr[i])$p.value
#   
# }

pval_0hr <- t.test(grouped_median$median_LI_0hr, grouped_median$median_NL_0hr)$p.value
pval_1hr <- t.test(grouped_median$median_LI_1hr, grouped_median$median_NL_1hr)$p.value
pval_6hr <- t.test(grouped_median$median_LI_6hr, grouped_median$median_NL_6hr)$p.value
pval_9hr <- t.test(grouped_median$median_LI_9hr, grouped_median$median_NL_9hr)$p.value
pval_D1 <- t.test(grouped_median$median_LI_D1, grouped_median$median_NL_D1)$p.value
pval_D14 <- t.test(grouped_median$median_LI_D14, grouped_median$median_NL_D14)$p.value
pval_D3 <- t.test(grouped_median$median_LI_D3, grouped_median$median_NL_D3)$p.value
pval_D7 <- t.test(grouped_median$median_LI_D7, grouped_median$median_NL_D7)$p.value

# combine FC and pval into each time dataframes
grouped_0hr <- data.frame(FC_0hr,pval_0hr)
grouped_1hr <- data.frame(FC_1hr,pval_1hr)
grouped_6hr <- data.frame(FC_6hr,pval_6hr)
grouped_9hr <- data.frame(FC_9hr,pval_9hr)
grouped_D1 <- data.frame(FC_D1,pval_D1)
grouped_D14 <- data.frame(FC_D14,pval_D14)
grouped_D3 <- data.frame(FC_D3,pval_D3)
grouped_D7 <- data.frame(FC_D7,pval_D7)

#plot each time 
plot(grouped_1hr$FC_1hr, grouped_1hr$pval_1hr)

# # plot volcano plot
# EnhancedVolcano(grouped_1hr,
#                 lab = rownames(grouped_median_pval),
#                 x = 'FC_1hr',
#                 y = 'pval_1hr',
#                 title = 'LI / NL',
#                 pCutoff = 0.05,
#                 FCcutoff = 1.0,
#                 pointSize = 3.0,
#                 labSize = 6.0)

#######################################################################

# grouped_median[,2:17] = log2(as.numeric(grouped_median[,2:17]))

         
# grouped_median[2:17] <- grouped_median[2:17] %>% mutate_if(is.character,as.numeric)
# grouped_median[2:17] <- as.character(grouped_median[2:17])

# grouped_median[2:17] <- apply(grouped_median[,2:17], 2, FUN=log2)

##########################################################################


########### group vs group comparisons (used mean to compare) ############
grouped_median$group1mean <- as.numeric(apply(grouped_median[,2:9], 1, FUN=mean))
grouped_median$group2mean <- as.numeric(apply(grouped_median[,10:17], 1, FUN=mean))
grouped_median$FC <- (group1mean/group2mean)
grouped_median$log2FC <- log(FC,2)

# pvalue <- apply(grouped_median, 1, function(x) { t.test(x[,2:9], x[,10:17])$p.value } )

# calculate p values for each row across control/disease groups
pval <- c()

for (i in 1:nrow(grouped_median)) {
  pval[i] <- t.test(grouped_median[,2:9][i], grouped_median[,10:17][i])$p.value
}

# combine pval to df and calculate -log10(pval)
pval <- data.frame(pval)
grouped_median_pval = cbind(grouped_median,pval)
grouped_median_pval$pval_log <- -log10(grouped_median_pval$pval)

plot(grouped_median_pval$log2FC, grouped_median_pval$pval_log)
# 
# grouped_median_pval$change <- ifelse(grouped_median_pval$log2FC> 1.5 & grouped_median_pval$pval_log<0.05,"UP",
#                               ifelse(grouped_median_pval$log2FC< 0.67 & grouped_median_pval$pval_log>0.05,"DOWN",
#                                      "nochange"))

# plot volcano plot
EnhancedVolcano(grouped_median_pval,
                lab = rownames(grouped_median_pval),
                x = 'log2FC',
                y = 'pval',
                title = 'LI / NL',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0)
