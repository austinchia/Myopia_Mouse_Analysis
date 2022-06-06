library(readxl)
library(MetaboAnalystR)
library(dplyr)
library(data.table)
library(EnhancedVolcano)
library(Mfuzz)
library(berryFunctions)
library(destiny)


#============ Metaboanalyst ================
{
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

}

# reading in csv of 3 groups

grouped_s1 <- fread("grouped_1.csv",sep=',')
grouped_s2 <- fread("grouped_2.csv",sep=',')
grouped_s3 <- fread("grouped_3.csv",sep=',')

grouped_combined <- left_join(grouped_s3,grouped_s2,by = 'Accession')
grouped_combined_2 <- left_join(grouped_combined,grouped_s1,by = 'Accession')
grouped_combined_2_no_na <- na.omit(grouped_combined_2)

write.csv(grouped_combined_2_no_na,"grouped_combined.csv", row.names = FALSE)

#============ Creates 7 Volcano Plots - replicates Metaboanalyst ===================

grouped_combined_GS <- fread("grouped_combined_GS_accounted.csv",sep=',')

# create groupings for LI and NL to insert at row 1

# selects columns for each timing (ie 1hr)
{
  # selects columns for 0hr
  grouped_0hr <- grouped_combined_GS %>%
    select(c('Gene Symbol', 'S1_LI_0hr', 'S2_LI_0hr', 'S3_LI_0hr', 'S1_NL_0hr', 'S2_NL_0hr', 'S3_NL_0hr')) %>%
    
    # gets mean of LI and NL
    group_by(`Gene Symbol`) %>%
    mutate('LI_0hr_mean' = mean(c(`S1_LI_0hr`,`S2_LI_0hr`,`S3_LI_0hr`))) %>%
    mutate('NL_0hr_mean' = mean(c(`S1_NL_0hr`,`S2_NL_0hr`,`S3_NL_0hr`))) %>%
    
    # calculates log2 fold change
    mutate('FC' = LI_0hr_mean/NL_0hr_mean) %>%
    mutate('Log2FC' = log2(FC)) %>%
    
    # calculates p value  
    mutate('p_value' = as.numeric(t.test(c(`S1_LI_0hr`,`S2_LI_0hr`,`S3_LI_0hr`), c(`S1_NL_0hr`,`S2_NL_0hr`,`S3_NL_0hr`))$p.value))
  
  
  # selects columns for 1hr
  grouped_1hr <- grouped_combined_GS %>%
    select(c('Gene Symbol', 'S1_LI_1hr', 'S2_LI_1hr', 'S3_LI_1hr', 'S1_NL_1hr', 'S2_NL_1hr', 'S3_NL_1hr')) %>%
    
    # gets mean of LI and NL
    group_by(`Gene Symbol`) %>%
    mutate('LI_1hr_mean' = mean(c(`S1_LI_1hr`,`S2_LI_1hr`,`S3_LI_1hr`))) %>%
    mutate('NL_1hr_mean' = mean(c(`S1_NL_1hr`,`S2_NL_1hr`,`S3_NL_1hr`))) %>%
    
    # calculates log2 fold change
    mutate('FC' = LI_1hr_mean/NL_1hr_mean) %>%
    mutate('Log2FC' = log2(FC)) %>%
  
    # calculates p value  
    mutate('p_value' = as.numeric(t.test(c(`S1_LI_1hr`,`S2_LI_1hr`,`S3_LI_1hr`), c(`S1_NL_1hr`,`S2_NL_1hr`,`S3_NL_1hr`))$p.value))
  
  # selects columns for 6hr
  grouped_6hr <- grouped_combined_GS %>%
    select(c('Gene Symbol', 'S1_LI_6hr', 'S2_LI_6hr', 'S3_LI_6hr', 'S1_NL_6hr', 'S2_NL_6hr', 'S3_NL_6hr')) %>%
    
    # gets mean of LI and NL
    group_by(`Gene Symbol`) %>%
    mutate('LI_6hr_mean' = mean(c(`S1_LI_6hr`,`S2_LI_6hr`,`S3_LI_6hr`))) %>%
    mutate('NL_6hr_mean' = mean(c(`S1_NL_6hr`,`S2_NL_6hr`,`S3_NL_6hr`))) %>%
    
    # calculates log2 fold change
    mutate('FC' = LI_6hr_mean/NL_6hr_mean) %>%
    mutate('Log2FC' = log2(FC)) %>%
    
    # calculates p value  
    mutate('p_value' = as.numeric(t.test(c(`S1_LI_6hr`,`S2_LI_6hr`,`S3_LI_6hr`), c(`S1_NL_6hr`,`S2_NL_6hr`,`S3_NL_6hr`))$p.value))
  
  # selects columns for 9hr
  grouped_9hr <- grouped_combined_GS %>%
    select(c('Gene Symbol', 'S1_LI_9hr', 'S2_LI_9hr', 'S3_LI_9hr', 'S1_NL_9hr', 'S2_NL_9hr', 'S3_NL_9hr')) %>%

    # gets mean of LI and NL
    group_by(`Gene Symbol`) %>%
    mutate('LI_9hr_mean' = mean(c(`S1_LI_9hr`,`S2_LI_9hr`,`S3_LI_9hr`))) %>%
    mutate('NL_9hr_mean' = mean(c(`S1_NL_9hr`,`S2_NL_9hr`,`S3_NL_9hr`))) %>%
    
    # calculates log2 fold change
    mutate('FC' = LI_9hr_mean/NL_9hr_mean) %>%
    mutate('Log2FC' = log2(FC)) %>%
    
    # calculates p value  
    mutate('p_value' = as.numeric(t.test(c(`S1_LI_9hr`,`S2_LI_9hr`,`S3_LI_9hr`), c(`S1_NL_9hr`,`S2_NL_9hr`,`S3_NL_9hr`))$p.value))
  
  # selects columns for D1
  grouped_D1 <- grouped_combined_GS %>%
    select(c('Gene Symbol', 'S1_LI_D1', 'S2_LI_D1', 'S3_LI_D1', 'S1_NL_D1', 'S2_NL_D1', 'S3_NL_D1')) %>%
    
    # gets mean of LI and NL
    group_by(`Gene Symbol`) %>%
    mutate('LI_D1_mean' = mean(c(`S1_LI_D1`,`S2_LI_D1`,`S3_LI_D1`))) %>%
    mutate('NL_D1_mean' = mean(c(`S1_NL_D1`,`S2_NL_D1`,`S3_NL_D1`))) %>%
    
    # calculates log2 fold change
    mutate('FC' = LI_D1_mean/NL_D1_mean) %>%
    mutate('Log2FC' = log2(FC)) %>%
    
    # calculates p value  
    mutate('p_value' = as.numeric(t.test(c(`S1_LI_D1`,`S2_LI_D1`,`S3_LI_D1`), c(`S1_NL_D1`,`S2_NL_D1`,`S3_NL_D1`))$p.value))
  
  # selects columns for D14
  grouped_D14 <- grouped_combined_GS %>%
    select(c('Gene Symbol', 'S1_LI_D14', 'S2_LI_D14', 'S3_LI_D14', 'S1_NL_D14', 'S2_NL_D14', 'S3_NL_D14')) %>%
    
    # gets mean of LI and NL
    group_by(`Gene Symbol`) %>%
    mutate('LI_D14_mean' = mean(c(`S1_LI_D14`,`S2_LI_D14`,`S3_LI_D14`))) %>%
    mutate('NL_D14_mean' = mean(c(`S1_NL_D14`,`S2_NL_D14`,`S3_NL_D14`))) %>%
    
    # calculates log2 fold change
    mutate('FC' = LI_D14_mean/NL_D14_mean) %>%
    mutate('Log2FC' = log2(FC)) %>%
    
    # calculates p value  
    mutate('p_value' = as.numeric(t.test(c(`S1_LI_D14`,`S2_LI_D14`,`S3_LI_D14`), c(`S1_NL_D14`,`S2_NL_D14`,`S3_NL_D14`))$p.value))
  
  # selects columns for D3
  grouped_D3 <- grouped_combined_GS %>%
    select(c('Gene Symbol', 'S1_LI_D3', 'S2_LI_D3', 'S3_LI_D3', 'S1_NL_D3', 'S2_NL_D3', 'S3_NL_D3')) %>%

    # gets mean of LI and NL
    group_by(`Gene Symbol`) %>%
    mutate('LI_D3_mean' = mean(c(`S1_LI_D3`,`S2_LI_D3`,`S3_LI_D3`))) %>%
    mutate('NL_D3_mean' = mean(c(`S1_NL_D3`,`S2_NL_D3`,`S3_NL_D3`))) %>%
    
    # calculates log2 fold change
    mutate('FC' = LI_D3_mean/NL_D3_mean) %>%
    mutate('Log2FC' = log2(FC)) %>%
    
    # calculates p value  
    mutate('p_value' = as.numeric(t.test(c(`S1_LI_D3`,`S2_LI_D3`,`S3_LI_D3`), c(`S1_NL_D3`,`S2_NL_D3`,`S3_NL_D3`))$p.value))
  
  # selects columns for D7
  grouped_D7 <- grouped_combined_GS %>%
    select(c('Gene Symbol', 'S1_LI_D7', 'S2_LI_D7', 'S3_LI_D7', 'S1_NL_D7', 'S2_NL_D7', 'S3_NL_D7')) %>%
    
    # gets mean of LI and NL
    group_by(`Gene Symbol`) %>%
    mutate('LI_D7_mean' = mean(c(`S1_LI_D7`,`S2_LI_D7`,`S3_LI_D7`))) %>%
    mutate('NL_D7_mean' = mean(c(`S1_NL_D7`,`S2_NL_D7`,`S3_NL_D7`))) %>%
    
    # calculates log2 fold change
    mutate('FC' = LI_D7_mean/NL_D7_mean) %>%
    mutate('Log2FC' = log2(FC)) %>%
    
    # calculates p value  
    mutate('p_value' = as.numeric(t.test(c(`S1_LI_D7`,`S2_LI_D7`,`S3_LI_D7`), c(`S1_NL_D7`,`S2_NL_D7`,`S3_NL_D7`))$p.value))

}


# plots volcano plot
vol_plot <- function(x) {
  EnhancedVolcano(x,
                  lab = rownames(x),
                  x = 'Log2FC',
                  y = 'p_value',
                  title = 'LI / NL',
                  pCutoff = 0.05,
                  FCcutoff = 1.0,
                  pointSize = 3.0,
                  labSize = 3.0)
}

# exports plots volcano plots for 7 timings (ie 1hr)
{ 

png(file="vplot_0hr.png")
vol_plot(grouped_0hr)
dev.off()

png(file="vplot_1hr.png")
vol_plot(grouped_1hr)
dev.off()

png(file="vplot_6hr.png")
vol_plot(grouped_6hr)
dev.off()

png(file="vplot_9hr.png")
vol_plot(grouped_9hr)
dev.off()

png(file="vplot_D1.png")
vol_plot(grouped_D1)
dev.off()

png(file="vplot_D3.png")
vol_plot(grouped_D3)
dev.off()

png(file="vplot_D7.png")
vol_plot(grouped_D7)
dev.off()

png(file="vplot_D14.png")
vol_plot(grouped_D14)
dev.off()
  
}

#================ Mfuzz ==========================

#======== fuzz on all sets
# convert dataframe to expression set
grouped_combined_GS_eSet <- as.ExpressionSet(grouped_combined_GS)

# scaling data
grouped_combined_GS_eSet.s <- standardise(grouped_combined_GS_eSet)

# estimating the fuzzifier
m1 <- mestimate(grouped_combined_GS_eSet.s)

# determining no of clusters
Dmin(grouped_combined_GS_eSet.s, m=m1, crange=seq(2,20,1), repeats=3, visu=TRUE)

cl <- mfuzz(grouped_combined_GS_eSet,c=10,m=1.25)
mfuzz.plot(grouped_combined_GS_eSet,cl=cl,mfrow=c(5,5),time.labels=grouped_combined_GS[2,])

#========== only fuzz on set 3
grouped_combined_GS_S3 <- fread("grouped_combined_GS_accounted_Mfuzz.csv",sep=',')
# covert dataframe to expression set
grouped_combined_GS_S3_eSet <- as.ExpressionSet(grouped_combined_GS_S3)

# scaling data
grouped_combined_GS_S3_eSet.s <- standardise(grouped_combined_GS_S3_eSet)

# estimating the fuzzifier
m1 <- mestimate(grouped_combined_GS_S3_eSet.s)

# determining no of clusters
Dmin(grouped_combined_GS_S3_eSet.s, m=m1, crange=seq(2,6), repeats=3, visu=TRUE)

# plot fuzz plot
time_mfuzz <- c(0,60,360,540,1440,4320,10080,20160)
cl <- mfuzz(grouped_combined_GS_S3_eSet.s,c=4,m=1.25)

mfuzz.plot(grouped_combined_GS_S3_eSet.s,
           cl=cl,
           mfrow=c(3,3),
           time.labels = seq(0,20160,100),
           min.mem=0.5)
