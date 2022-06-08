library(readxl)
library(MetaboAnalystR)
library(dplyr)
library(data.table)
library(EnhancedVolcano)
library(Mfuzz)
library(berryFunctions)
library(destiny)
library(tidyr)
library(stringr)
library(marray)
library(ggplot2)
library(janitor)

# ============== 1. Reads Raw Data ===================

# reads S1 raw data
Retina_WP_S1 <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S1', na = c("", "NA")) %>%
  select(c(`Accession`,`Abundance Ratios`))
  
# reads S2 raw data
Retina_WP_S2 <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S2', na = c("", "NA")) %>%
  select(c(`Accession`,`Abundance Ratios`))

# reads S3 raw data
Retina_WP_S3 <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S3', na = c("", "NA")) %>%
  select(c(`Accession`,`Abundance Ratios`))

# ============== 2. Data Manipulation ===============
# S1 Data manipulation
{
  # splits string into columns
  Retina_WP_S1_split <- str_split_fixed(as.character(Retina_WP_S1$`Abundance Ratios`), ';',15)
  
  # adds columns to original data
  Retina_WP_S1 <- cbind(Retina_WP_S1,Retina_WP_S1_split)
  colnames(Retina_WP_S1) <- c('Accession', 'Abundance Ratios', 'S1_LI_1hr','S1_LI_6hr','S1_LI_9hr','S1_LI_D1','S1_LI_D14','S1_LI_D3','S1_LI_D7','S1_NL_0hr','S1_NL_1hr','S1_NL_6hr','S1_NL_9hr','S1_NL_D1','S1_NL_D14','S1_NL_D3','S1_NL_D7')
  
  # removes abundance column
  Retina_WP_S1 <- select(Retina_WP_S1, -c(`Abundance Ratios`)) %>%
    mutate(S1_LI_0hr = S1_NL_0hr) %>%
    relocate(S1_LI_0hr, .after = `Accession`) %>%
    relocate(S1_LI_D14, .after = `S1_LI_D7`) %>%
    relocate(S1_NL_D14, .after = `S1_NL_D7`)
  
}
# S2 Data manipulation
{
  # splits string into columns
  Retina_WP_S2_split <- str_split_fixed(as.character(Retina_WP_S2$`Abundance Ratios`), ';',15)
  
  # adds columns to original data
  Retina_WP_S2 <- cbind(Retina_WP_S2,Retina_WP_S2_split)
  
  colnames(Retina_WP_S2) <- c('Accession', 'Abundance Ratios', 'S2_LI_1hr','S2_LI_6hr','S2_LI_9hr','S2_LI_D1','S2_LI_D14','S2_LI_D3','S2_LI_D7','S2_NL_0hr','S2_NL_1hr','S2_NL_6hr','S2_NL_9hr','S2_NL_D1','S2_NL_D14','S2_NL_D3','S2_NL_D7')
  
  # removes abundance column
  Retina_WP_S2 <- select(Retina_WP_S2, -c(`Abundance Ratios`)) %>%
    mutate(S2_LI_0hr = S2_NL_0hr) %>%
    relocate(S2_LI_0hr, .after = `Accession`) %>%
    relocate(S2_LI_D14, .after = `S2_LI_D7`) %>%
    relocate(S2_NL_D14, .after = `S2_NL_D7`)
  }
# S3 Data manipulation
{
  # splits string into columns
  Retina_WP_S3_split <- str_split_fixed(as.character(Retina_WP_S3$`Abundance Ratios`), ';',15)
  
  # adds columns to original data
  Retina_WP_S3 <- cbind(Retina_WP_S3,Retina_WP_S3_split)
  
  colnames(Retina_WP_S3) <- c('Accession', 'Abundance Ratios', 'S3_LI_1hr','S3_LI_6hr','S3_LI_9hr','S3_LI_D1','S3_LI_D14','S3_LI_D3','S3_LI_D7','S3_NL_0hr','S3_NL_1hr','S3_NL_6hr','S3_NL_9hr','S3_NL_D1','S3_NL_D14','S3_NL_D3','S3_NL_D7')
  
  # removes abundance column
  Retina_WP_S3 <- select(Retina_WP_S3, -c(`Abundance Ratios`)) %>%
    mutate(S3_LI_0hr = S3_NL_0hr) %>%
    relocate(S3_LI_0hr, .after = `Accession`) %>%
    relocate(S3_LI_D14, .after = `S3_LI_D7`) %>%
    relocate(S3_NL_D14, .after = `S3_NL_D7`)
}

# combines all 3 sets into 1 dataset 
# (left joins to S3 because it has most no of proteins)
ratio_combined <- left_join(Retina_WP_S3, Retina_WP_S2, by = 'Accession') %>%
  left_join(Retina_WP_S1, by = 'Accession') %>%
  na.omit()

ratio_combined <- ratio_combined %>% 
  mutate_all(na_if,"")

# splits accession number (ie Q9JHU4-1)
ratio_combined$Accession <- sapply(strsplit(ratio_combined$Accession,"-"), `[`, 1)

# exports accession numbers to upload to Uniprot
fwrite(data.frame(ratio_combined$Accession), "test.csv", sep = ",")

gene_symbol <- fread("test_map.csv",sep=',')

# splits gene symbol by break
gene_symbol_map <- data.frame(str_split_fixed(gene_symbol$`From	To`, '\t',2))
colnames(gene_symbol_map) <- c("Accession", "Gene Symbol") 

ratio_combined_test <- left_join(ratio_combined, gene_symbol_map, by="Accession") %>%
  relocate(`Gene Symbol`, .after = `Accession`)


# =========== Volcano Plots - replaces Metaboanalyst ============= ===================

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



# =========== Mfuzz Plots (Uses Grouped Abundance) ==========================
# ============ 1. Selects Columns From Main Grouped Matrix =========
#== selects S1 for fuzz
fuzz_S1_LI <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S1_LI_0hr`,	`S1_LI_1hr`,	`S1_LI_6hr`,	`S1_LI_9hr`,	`S1_LI_D1`,	`S1_LI_D14`,	`S1_LI_D3`,	`S1_LI_D7`)
fuzz_S1_NL <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S1_NL_0hr`,	`S1_NL_1hr`,	`S1_NL_6hr`,	`S1_NL_9hr`,	`S1_NL_D1`,	`S1_NL_D14`,	`S1_NL_D3`,	`S1_NL_D7`)

#== selects S2 for fuzz
fuzz_S2_LI <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S2_LI_0hr`,	`S2_LI_1hr`,	`S2_LI_6hr`,	`S2_LI_9hr`,	`S2_LI_D1`,	`S2_LI_D14`,	`S2_LI_D3`,	`S2_LI_D7`)
fuzz_S2_NL <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S2_NL_0hr`,	`S2_NL_1hr`,	`S2_NL_6hr`,	`S2_NL_9hr`,	`S2_NL_D1`,	`S2_NL_D14`,	`S2_NL_D3`,	`S2_NL_D7`)

#== selects S3 for fuzz
fuzz_S3_LI <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S3_LI_0hr`,	`S3_LI_1hr`,	`S3_LI_6hr`,	`S3_LI_9hr`,	`S3_LI_D1`,	`S3_LI_D14`,	`S3_LI_D3`,	`S3_LI_D7`)
fuzz_S3_NL <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S3_NL_0hr`,	`S3_NL_1hr`,	`S3_NL_6hr`,	`S3_NL_9hr`,	`S3_NL_D1`,	`S3_NL_D14`,	`S3_NL_D3`,	`S3_NL_D7`)

# ============ 2. Creates Timepoints and Binds to Original Dataframe ====
# == S3 NL
# timepoint <- data.frame(t(c("NA",0,1,6,9,24,72,168,336)))
# colnames(timepoint) <- colnames(fuzz_S3_NL)
# temp_table <- rbind(timepoint, fuzz_S3_NL)
# row.names(temp_table)[1]<-"time"
# tmp <- tempfile()
# write.table(temp_table,file=tmp, sep='\t', quote = F,col.names=NA)

# == function to create timepoints, convert to eset

create_timepoints <- function(x) {
  # creates timepoints
  timepoint <- data.frame(t(c("NA",0,1,6,9,24,72,168,336)))
  colnames(timepoint) <- colnames(x)
  
  # creates temp table
  temp_table <- rbind(timepoint, x) 
  
  # sets rownames
  row.names(temp_table)[1]<-"time" 
  
  # stores as tmp format to read into table2eset
  tmp <- tempfile() 
  write.table(temp_table,file=tmp, sep='\t', quote = F,col.names=NA)
  x <- table2eset(file=tmp)
  
}

# runs create_timepoints
S1_LI_eSet <- create_timepoints(fuzz_S1_LI)
S1_NL_eSet <- create_timepoints(fuzz_S1_NL)

S2_LI_eSet <- create_timepoints(fuzz_S2_LI)
S2_NL_eSet <- create_timepoints(fuzz_S2_NL)

S3_LI_eSet <- create_timepoints(fuzz_S3_LI)
S3_NL_eSet <- create_timepoints(fuzz_S3_NL)

# ============ 3. Estimates Fuzzifier (ie m1) ================
m1_S1_LI <- mestimate(S1_LI_eSet)
m1_S1_NL <- mestimate(S1_NL_eSet)
m1_S2_LI <- mestimate(S2_LI_eSet)
m1_S2_NL <- mestimate(S2_NL_eSet)
m1_S3_LI <- mestimate(S3_LI_eSet)
m1_S3_NL <- mestimate(S3_NL_eSet)

# ============ 4. Plots Mfuzz Plots ======================

plot_mfuzz <- function(x) {
  cl <- mfuzz(x,c=12,m=m1)
  mfuzz.plot2(x,
              cl=cl,
              mfrow=c(3,3),
              time.labels = c(0,1,6,9,24,72,168,336),
              col.main = ,
              min.mem=0.5,
  )
}

# plots mfuzz for Set 1
plot_mfuzz(S1_LI_eSet)
plot_mfuzz(S1_NL_eSet)

# plots mfuzz for Set 2
plot_mfuzz(S2_LI_eSet)
plot_mfuzz(S2_NL_eSet)

# plots mfuzz for Set 3
plot_mfuzz(S3_LI_eSet)
plot_mfuzz(S3_NL_eSet)

# creates correlation matrix between cluster centroids
# (no more than 0.85)
correlation_matrix <- data.frame(cor(t(cl[[1]])))

#extracts membership values 
acore <- acore(grouped_combined_GS_S3_eSet.s,cl,min.acore=0)

# pull out the scores for the cluster assignments 
# (where the assignment is based on the top scoring cluster)
acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))


# =========== Metaboanalyst Volcano Plot ================
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

