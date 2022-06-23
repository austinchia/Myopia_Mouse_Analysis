# About
# This script reads in raw data > preprocesses data > to form a data matrix
# Data matrix is used for further stats analysis
# ============== Loads Packages =======
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
library(IMIFA)
library(tidyverse)
library(ggVennDiagram)
library(ggvenn)

# ====== A) Prepares Data Matrix for Grouped Abundance =====

# reads S1 raw data
Retina_WP_S1_grouped <- read_excel('Myopia_retina_Phospho_results.xlsx', sheet = 'Phospho_Ret_S1', na = c("", "NA")) %>%
  select(c(`Annotated Sequence`,
           `Modifications`,
           `Master Protein Accessions`,
           `Abundances (Grouped)`)) %>%
  na.omit()


# reads S2 raw data
Retina_WP_S2_grouped <- read_excel('Myopia_retina_Phospho_results.xlsx', sheet = 'Phospho_Ret_S2', na = c("", "NA")) %>%
  select(c(`Annotated Sequence`,
           `Modifications`,
           `Master Protein Accessions`,
           `Abundances (Grouped)`)) %>%
  na.omit()

# reads S3 raw data
Retina_WP_S3_grouped <- read_excel('Myopia_retina_Phospho_results.xlsx', sheet = 'Phospho_Ret_S3', na = c("", "NA")) %>%
  select(c(`Annotated Sequence`,
           `Modifications`,
           `Master Protein Accessions`,
           `Abundances (Grouped)`)) %>%
  na.omit()

# ============== 2. Manipulates Data  ===============
# S1 Data manipulation
{
  # splits string into columns
  Retina_WP_S1_grouped_split <- str_split_fixed(as.character(Retina_WP_S1_grouped$`Abundances (Grouped)`), ';',16)
  
  # adds columns to original data
  Retina_WP_S1_grouped <- cbind(Retina_WP_S1_grouped, Retina_WP_S1_grouped_split)
  colnames(Retina_WP_S1_grouped) <- c("Annotated_Sequence",
                                      "Modifications",
                                      "Master_Protein_Accessions",
                                      "Abundances (Grouped)",
                                      "GC_S1",
                                      'S1_LI_1hr','S1_LI_6hr','S1_LI_9hr','S1_LI_D1','S1_LI_D14','S1_LI_D3','S1_LI_D7','S1_NL_0hr','S1_NL_1hr','S1_NL_6hr','S1_NL_9hr','S1_NL_D1','S1_NL_D14','S1_NL_D3','S1_NL_D7')
  
  # removes abundance column
  Retina_WP_S1_grouped <- select(Retina_WP_S1_grouped, -c(`Abundances (Grouped)`)) %>%
    mutate(S1_LI_0hr = S1_NL_0hr) %>%
    relocate(S1_LI_0hr, .after = `Master_Protein_Accessions`) %>%
    relocate(S1_LI_D14, .after = `S1_LI_D7`) %>%
    relocate(S1_NL_D14, .after = `S1_NL_D7`)
  
}

# S2 Data manipulation
{
  # splits string into columns
  Retina_WP_S2_grouped_split <- str_split_fixed(as.character(Retina_WP_S2_grouped$`Abundances (Grouped)`), ';',16)
  
  # adds columns to original data
  Retina_WP_S2_grouped <- cbind(Retina_WP_S2_grouped, Retina_WP_S2_grouped_split)
  colnames(Retina_WP_S2_grouped) <- c("Annotated_Sequence",
                                      "Modifications",
                                      "Master_Protein_Accessions",
                                      "Abundances (Grouped)",
                                      "GC_S2",
                                      'S2_LI_1hr','S2_LI_6hr','S2_LI_9hr','S2_LI_D1','S2_LI_D14','S2_LI_D3','S2_LI_D7','S2_NL_0hr','S2_NL_1hr','S2_NL_6hr','S2_NL_9hr','S2_NL_D1','S2_NL_D14','S2_NL_D3','S2_NL_D7')
  
  # removes abundance column
  Retina_WP_S2_grouped <- select(Retina_WP_S2_grouped, -c(`Abundances (Grouped)`)) %>%
    mutate(S2_LI_0hr = S2_NL_0hr) %>%
    relocate(S2_LI_0hr, .after = `Master_Protein_Accessions`) %>%
    relocate(S2_LI_D14, .after = `S2_LI_D7`) %>%
    relocate(S2_NL_D14, .after = `S2_NL_D7`)
  
}

# S3 Data manipulation
{
  # splits string into columns
  Retina_WP_S3_grouped_split <- str_split_fixed(as.character(Retina_WP_S3_grouped$`Abundances (Grouped)`), ';',16)
  
  # adds columns to original data
  Retina_WP_S3_grouped <- cbind(Retina_WP_S3_grouped, Retina_WP_S3_grouped_split)
  colnames(Retina_WP_S3_grouped) <- c("Annotated_Sequence",
                                      "Modifications",
                                      "Master_Protein_Accessions",
                                      "Abundances (Grouped)", 
                                      "GC_S3",
                                      'S3_LI_1hr','S3_LI_6hr','S3_LI_9hr','S3_LI_D1','S3_LI_D14','S3_LI_D3','S3_LI_D7','S3_NL_0hr','S3_NL_1hr','S3_NL_6hr','S3_NL_9hr','S3_NL_D1','S3_NL_D14','S3_NL_D3','S3_NL_D7')
  
  # removes abundance column
  Retina_WP_S3_grouped <- select(Retina_WP_S3_grouped, -c(`Abundances (Grouped)`)) %>%
    mutate(S3_LI_0hr = S3_NL_0hr) %>%
    relocate(S3_LI_0hr, .after = `Master_Protein_Accessions`) %>%
    relocate(S3_LI_D14, .after = `S3_LI_D7`) %>%
    relocate(S3_NL_D14, .after = `S3_NL_D7`)
  
}

# ============== 3. Combines 3 Sets And Cleans Data ======

# combines all 3 sets into 1 dataset 
# (left joins to S1 because it has most no of proteins, followed by S1)
ratio_combined <- Retina_WP_S2_grouped %>%
  
  # merges S1
  left_join(Retina_WP_S1_grouped, by = 'Annotated_Sequence') %>%
  
  #merges S2
  left_join(Retina_WP_S3_grouped, by = 'Annotated_Sequence') %>%
  
  # removes duplicate peptides
  distinct(`Annotated_Sequence`, .keep_all = TRUE) %>%
  na.omit()

# convert blanks to NA
ratio_combined <- ratio_combined %>% 
  mutate_all(na_if,"")

# splits accession by ";" delimiter (ie "Q9JHU4-1; Q9JHU4" --> "Q9JHU4-1")
ratio_combined$Master_Protein_Accessions.x <- sapply(strsplit(ratio_combined$Master_Protein_Accessions.x,";"), `[`, 1)

# splits accession by "-" delimiter (ie "Q9JHU4-1; Q9JHU4" --> "Q9JHU4-1")
ratio_combined$Master_Protein_Accessions.x <- sapply(strsplit(ratio_combined$Master_Protein_Accessions.x,"-"), `[`, 1)

# exports accession numbers to upload to Uniprot
fwrite(data.frame(ratio_combined$`Master_Protein_Accessions.x`), "Phospho_Accession.csv", sep = ",")

# ============== 4. Combines Uniprot Data To Combined Matrix =====
# reads in Gene Symbol table downloaded from Uniprot
gene_symbol <- fread("Phospho_Protein_Accession_Map.csv",sep=',')

# splits gene symbol by break
gene_symbol_map <- data.frame(str_split_fixed(gene_symbol$`From	To`, '\t',2))
colnames(gene_symbol_map) <- c("Master_Protein_Accessions", "Gene Symbol") 

# merges gene symbol column to main df
ratio_combined_no_na <- ratio_combined %>%
  
  # merges gene symbol column to main df
  left_join(gene_symbol_map,
            by = "Master_Protein_Accessions") %>%

  # relocates columns and removes NAs
  relocate(c(`Modifications`,
             `Master_Protein_Accessions`,
             `Gene Symbol`),
           .before = `S2_LI_0hr`) %>%
  na.omit() %>%
  
  # adds number to the end of duplicate gene symbols (ie Sptbn1-2)
  group_by(`Gene Symbol`) %>%
  mutate(`GS_count` = 1:n()) %>%
  mutate(`Gene Symbol` = ifelse(`GS_count` == 1, 
                                `Gene Symbol`, 
                                paste0(`Gene Symbol`, "-", `GS_count`))) %>%
  
  # removes unused columns
  select(-c(`GS_count`,
            `Modifications`,
            `Master_Protein_Accessions`,
            `Modifications.y`,
            `Master_Protein_Accessions.y`
            )) 
  
# exports combined grouped abundance matrix to csv
fwrite(ratio_combined_no_na, "phospho_abund_grouped_combined.csv", sep = ",")

# ============== 5. Selects Columns From Main Grouped Matrix =========
grouped_combined_GS <- ratio_combined_no_na

# == selects S1 for fuzz
fuzz_S1_LI <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S1_LI_0hr`,	`S1_LI_1hr`,	`S1_LI_6hr`,	`S1_LI_9hr`,	`S1_LI_D1`,	`S1_LI_D3`,	`S1_LI_D7`,	`S1_LI_D14`) %>%
  # replaces 0 with NA
  na_if(0) %>%
  # removes NAs
  na.omit() %>%
  # set rownames as `Gene Symbol`
  column_to_rownames(., var = "Gene Symbol")

fuzz_S1_NL <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S1_NL_0hr`,	`S1_NL_1hr`,	`S1_NL_6hr`,	`S1_NL_9hr`,	`S1_NL_D1`,	`S1_NL_D3`,	`S1_NL_D7`,	`S1_NL_D14`) %>%
  # replaces 0 with NA
  na_if(0) %>%
  # removes NAs
  na.omit() %>%
  # set rownames as `Gene Symbol`
  column_to_rownames(., var = "Gene Symbol")

# == selects S2 for fuzz
fuzz_S2_LI <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S2_LI_0hr`,	`S2_LI_1hr`,	`S2_LI_6hr`,	`S2_LI_9hr`,	`S2_LI_D1`,	`S2_LI_D3`,	`S2_LI_D7`,	`S2_LI_D14`) %>%
  # replaces 0 with NA
  na_if(0) %>%
  # removes NAs
  na.omit() %>%
  # set rownames as `Gene Symbol`
  column_to_rownames(., var = "Gene Symbol")

fuzz_S2_NL <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S2_NL_0hr`,	`S2_NL_1hr`,	`S2_NL_6hr`,	`S2_NL_9hr`,	`S2_NL_D1`,	`S2_NL_D3`,	`S2_NL_D7`,	`S2_NL_D14`) %>%
  # replaces 0 with NA
  na_if(0) %>%
  # removes NAs
  na.omit() %>%
  # set rownames as `Gene Symbol`
  column_to_rownames(., var = "Gene Symbol")

# == selects S3 for fuzz
fuzz_S3_LI <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S3_LI_0hr`,	`S3_LI_1hr`,	`S3_LI_6hr`,	`S3_LI_9hr`,	`S3_LI_D1`,	`S3_LI_D3`,	`S3_LI_D7`,	`S3_LI_D14`) %>%
  # replaces 0 with NA
  na_if(0) %>%
  # removes NAs
  na.omit() %>%
  # set rownames as `Gene Symbol`
  column_to_rownames(., var = "Gene Symbol")

fuzz_S3_NL <- grouped_combined_GS %>%
  select(`Gene Symbol`, `S3_NL_0hr`,	`S3_NL_1hr`,	`S3_NL_6hr`,	`S3_NL_9hr`,	`S3_NL_D1`,	`S3_NL_D3`,	`S3_NL_D7`, `S3_NL_D14`) %>%
  # replaces 0 with NA
  na_if(0) %>%
  # removes NAs
  na.omit() %>%
  # set rownames as `Gene Symbol`
  column_to_rownames(., var = "Gene Symbol")


# == Combines all LI grouped abundance ==
{
  fuzz_S1_LI_with_GS <- grouped_combined_GS %>%
    select(`Gene Symbol`, `S1_LI_0hr`,	`S1_LI_1hr`,	`S1_LI_6hr`,	`S1_LI_9hr`,	`S1_LI_D1`,	`S1_LI_D3`,	`S1_LI_D7`,	`S1_LI_D14`) %>%
    # replaces 0 with NA
    na_if(0) %>%
    # removes NAs
    na.omit()
  
  fuzz_S2_LI_with_GS <- grouped_combined_GS %>%
    select(`Gene Symbol`, `S2_LI_0hr`,	`S2_LI_1hr`,	`S2_LI_6hr`,	`S2_LI_9hr`,	`S2_LI_D1`,	`S2_LI_D3`,	`S2_LI_D7`,	`S2_LI_D14`) %>%
    # replaces 0 with NA
    na_if(0) %>%
    # removes NAs
    na.omit()
  
  fuzz_S3_LI_with_GS <- grouped_combined_GS %>%
    select(`Gene Symbol`, `S3_LI_0hr`,	`S3_LI_1hr`,	`S3_LI_6hr`,	`S3_LI_9hr`,	`S3_LI_D1`,	`S3_LI_D3`,	`S3_LI_D7`,	`S3_LI_D14`) %>%
    # replaces 0 with NA
    na_if(0) %>%
    # removes NAs
    na.omit()
  
  fuzz_combined_LI <- left_join(fuzz_S3_LI_with_GS, fuzz_S2_LI_with_GS, by = 'Gene Symbol') %>%
    left_join(fuzz_S1_LI_with_GS, by = 'Gene Symbol') %>%
    na.omit() %>%
    column_to_rownames(., var = "Gene Symbol")
  }


# == Combines all NL grouped abundance ==
{
  fuzz_S1_NL_with_GS <- grouped_combined_GS %>%
    select(`Gene Symbol`, `S1_NL_0hr`,	`S1_NL_1hr`,	`S1_NL_6hr`,	`S1_NL_9hr`,	`S1_NL_D1`,	`S1_NL_D3`,	`S1_NL_D7`,	`S1_NL_D14`) %>%
    # replaces 0 with NA
    na_if(0) %>%
    # removes NAs
    na.omit()
  
  fuzz_S2_NL_with_GS <- grouped_combined_GS %>%
    select(`Gene Symbol`, `S2_NL_0hr`,	`S2_NL_1hr`,	`S2_NL_6hr`,	`S2_NL_9hr`,	`S2_NL_D1`,	`S2_NL_D3`,	`S2_NL_D7`,	`S2_NL_D14`) %>%
    # replaces 0 with NA
    na_if(0) %>%
    # removes NAs
    na.omit()
  
  fuzz_S3_NL_with_GS <- grouped_combined_GS %>%
    select(`Gene Symbol`, `S3_NL_0hr`,	`S3_NL_1hr`,	`S3_NL_6hr`,	`S3_NL_9hr`,	`S3_NL_D1`,	`S3_NL_D3`,	`S3_NL_D7`,	`S3_NL_D14`) %>%
    # replaces 0 with NA
    na_if(0) %>%
    # removes NAs
    na.omit()
  
  fuzz_combined_NL <- left_join(fuzz_S3_NL_with_GS, fuzz_S2_NL_with_GS, by = 'Gene Symbol') %>%
    left_join(fuzz_S1_NL_with_GS, by = 'Gene Symbol') %>%
    na.omit() %>%
    column_to_rownames(., var = "Gene Symbol")
}

# combines all 3 LI sets and calculates average
LI_average <- fuzz_combined_LI %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  mutate('LI_0hr_mean' = rowMeans(subset(., select = c(`S1_LI_0hr`,`S2_LI_0hr`,`S3_LI_0hr`)))) %>%
  mutate('LI_1hr_mean' = rowMeans(subset(., select = c(`S1_LI_1hr`,`S2_LI_1hr`,`S3_LI_1hr`)))) %>%
  mutate('LI_6hr_mean' = rowMeans(subset(., select = c(`S1_LI_6hr`,`S2_LI_6hr`,`S3_LI_6hr`)))) %>%
  mutate('LI_9hr_mean' = rowMeans(subset(., select = c(`S1_LI_9hr`,`S2_LI_9hr`,`S3_LI_9hr`)))) %>%
  mutate('LI_D1_mean' = rowMeans(subset(., select = c(`S1_LI_D1`,`S2_LI_D1`,`S3_LI_D1`)))) %>%
  mutate('LI_D3_mean' = rowMeans(subset(., select = c(`S1_LI_D3`,`S2_LI_D3`,`S3_LI_D3`)))) %>%
  mutate('LI_D7_mean' = rowMeans(subset(., select = c(`S1_LI_D7`,`S2_LI_D7`,`S3_LI_D7`)))) %>%
  mutate('LI_D14_mean' = rowMeans(subset(., select = c(`S1_LI_D14`,`S2_LI_D14`,`S3_LI_D14`)))) %>%
  select(`LI_0hr_mean`, `LI_1hr_mean`, `LI_6hr_mean`, `LI_9hr_mean`, `LI_D1_mean`, `LI_D3_mean`, `LI_D7_mean`, `LI_D14_mean`)

# combines all 3 NL sets and calculates average
NL_average <- fuzz_combined_NL %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  mutate('NL_0hr_mean' = rowMeans(subset(., select = c(`S1_NL_0hr`,`S2_NL_0hr`,`S3_NL_0hr`)))) %>%
  mutate('NL_1hr_mean' = rowMeans(subset(., select = c(`S1_NL_1hr`,`S2_NL_1hr`,`S3_NL_1hr`)))) %>%
  mutate('NL_6hr_mean' = rowMeans(subset(., select = c(`S1_NL_6hr`,`S2_NL_6hr`,`S3_NL_6hr`)))) %>%
  mutate('NL_9hr_mean' = rowMeans(subset(., select = c(`S1_NL_9hr`,`S2_NL_9hr`,`S3_NL_9hr`)))) %>%
  mutate('NL_D1_mean' = rowMeans(subset(., select = c(`S1_NL_D1`,`S2_NL_D1`,`S3_NL_D1`)))) %>%
  mutate('NL_D3_mean' = rowMeans(subset(., select = c(`S1_NL_D3`,`S2_NL_D3`,`S3_NL_D3`)))) %>%
  mutate('NL_D7_mean' = rowMeans(subset(., select = c(`S1_NL_D7`,`S2_NL_D7`,`S3_NL_D7`)))) %>%
  mutate('NL_D14_mean' = rowMeans(subset(., select = c(`S1_NL_D14`,`S2_NL_D14`,`S3_NL_D14`)))) %>%
  select(`NL_0hr_mean`, `NL_1hr_mean`, `NL_6hr_mean`, `NL_9hr_mean`, `NL_D1_mean`, `NL_D3_mean`, `NL_D7_mean`, `NL_D14_mean`)

# ============== 6. Creates Timepoints and Binds to Original Dataframe ====

# function to create timepoints, convert to eset
create_timepoints <- function(x) {
  # creates timepoints
  timepoint <- data.frame(t(c(0,1,6,9,24,72,168,336)))
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

# runs create_timepoints function
S1_LI_eSet <- create_timepoints(fuzz_S1_LI)
S1_NL_eSet <- create_timepoints(fuzz_S1_NL)

S2_LI_eSet <- create_timepoints(fuzz_S2_LI)
S2_NL_eSet <- create_timepoints(fuzz_S2_NL)

S3_LI_eSet <- create_timepoints(fuzz_S3_LI)
S3_NL_eSet <- create_timepoints(fuzz_S3_NL)

LI_average_eSet <- create_timepoints(LI_average)
NL_average_eSet <- create_timepoints(NL_average)

# ============== 7. Scales Data ========================

# scales data
S1_LI_eSet <- standardise(S1_LI_eSet)
S1_NL_eSet <- standardise(S1_NL_eSet)

S2_LI_eSet <- standardise(S2_LI_eSet)
S2_NL_eSet <- standardise(S2_NL_eSet)

S3_LI_eSet <- standardise(S3_LI_eSet)
S3_NL_eSet <- standardise(S3_NL_eSet)

# normalizes LI average
LI_average_eSet <- standardise(LI_average_eSet)
# normalizes LI average
NL_average_eSet <- standardise(NL_average_eSet)

# ============== 8. Estimates Fuzzifier (ie m1) ================

m1_S1_LI <- mestimate(S1_LI_eSet)
m1_S1_NL <- mestimate(S1_NL_eSet)
m1_S2_LI <- mestimate(S2_LI_eSet)
m1_S2_NL <- mestimate(S2_NL_eSet)
m1_S3_LI <- mestimate(S3_LI_eSet)
m1_S3_NL <- mestimate(S3_NL_eSet)

# estimates fuzzifier for LI and NL (all sets)
m1_average_LI <- mestimate(LI_average_eSet)
m1_average_NL <- mestimate(NL_average_eSet)

# plots scree plot - determine no of centroids (12 is the best)
# Dmin(LI_average_eSet, m=m1_average_LI, crange=seq(5,20,1), repeats=3, visu=TRUE)

# ============== 9. Plots Mfuzz Plots ======================

plot_mfuzz <- function(x,m1) {
  cl <- mfuzz(x,c=12,m=m1)
  mfuzz.plot2(x,
              cl=cl,
              mfrow=c(4,3),
              time.labels = c("0hr", "1hr", "6hr", "9hr", "D1", "D3", "D7", "D14"),
              col.main = ,
              ylab = "Abundance Changes",
              ylim.set=c(-3,3),
              min.mem=0.5,
              x11 = FALSE # popup appears when TRUE
  )
}

# plots mfuzz plots
{
  # # plots mfuzz for Set 1
  # png(file = "mfuzz_S1_LI_phospho.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S1_LI_eSet, m1_S1_LI)
  # dev.off()
  # 
  # png(file="mfuzz_S1_NL_phospho.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S1_NL_eSet, m1_S1_NL)
  # dev.off()
  # 
  # # plots mfuzz for Set 2
  # png(file="mfuzz_S2_LI_phospho.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S2_LI_eSet, m1_S2_LI)
  # dev.off()
  # 
  # png(file="mfuzz_S2_NL_phospho.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S2_NL_eSet, m1_S2_NL)
  # dev.off()
  # 
  # # plots mfuzz for Set 3
  # png(file="mfuzz_S3_LI_phospho.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S3_LI_eSet, m1_S3_LI)
  # dev.off()
  # 
  # png(file="mfuzz_S3_NL_phospho.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S3_NL_eSet, m1_S3_NL)
  # dev.off()
  # 
  # # plots mfuzz for LI (all sets)
  # png(file="mfuzz_LI_average_phospho.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(LI_average_eSet, m1_average_LI)
  # dev.off()
  # 
  # # plots mfuzz for NL (all sets)
  # png(file="mfuzz_NL_average_phospho.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(NL_average_eSet, m1_average_NL)
  # dev.off()
}

# ============== 10. Validates and Evaluates Mfuzz Model ==========

# creates correlation matrix between cluster centroids
# (no more than 0.85)
# correlation_matrix <- data.frame(cor(t(cl[[1]])))

# ============== 11. Extracts Gene Lists From Clusters ==========
# creates function to extract genes (acore list) in each cluster
get_genes <- function(x) {
  acore_x <- acore(x,cl,min.acore=0)
  do.call(rbind, lapply(seq_along(acore_x), 
                        function(i){ data.frame(Cluster=i, 
                                                acore_x[[i]])}))
}

# uses get_genes function to extract acore list
S1_LI_acore_list <- get_genes(S1_LI_eSet)
S1_NL_acore_list <- get_genes(S1_NL_eSet)

S2_LI_acore_list <- get_genes(S2_LI_eSet)
S2_NL_acore_list <- get_genes(S2_NL_eSet)

S3_LI_acore_list <- get_genes(S3_LI_eSet)
S3_NL_acore_list <- get_genes(S3_NL_eSet)

# extracts acore list for combined sets (LI and NL)
LI_acore_list <- LI_average_eSet %>%
  get_genes() %>%
  na.omit() %>%
  rename("Gene Symbol" = "NAME")  

# extracts acore list for combined sets (LI and NL)
NL_acore_list <- NL_average_eSet %>%
  get_genes() %>%
  na.omit() %>%
  rename("Gene Symbol" = "NAME")  

# exports lists of genes in clusters
fwrite(LI_acore_list, "LI_cluster_genes.csv", sep = ",")
fwrite(NL_acore_list, "NL_cluster_genes.csv", sep = ",")

# creates function to join acore list to abundance by gene symbol
combine_acore <- function(abundance_df, acore_list) {
  acore_combined <- left_join(rownames_to_column(abundance_df), 
                              acore_list, 
                              by=c("rowname" = "NAME"))
  names(acore_combined)[names(acore_combined) == 'rowname'] <- 'Gene Symbol'
  return(acore_combined)
}

# joins acore list with abundance ratios 
# joins S1
S1_LI_acore_list_combined <- combine_acore(fuzz_S1_LI, S1_LI_acore_list)
S1_NL_acore_list_combined <- combine_acore(fuzz_S1_NL, S1_NL_acore_list)
# joins S2
S2_LI_acore_list_combined <- combine_acore(fuzz_S2_LI, S2_LI_acore_list)
S2_NL_acore_list_combined <- combine_acore(fuzz_S2_NL, S2_NL_acore_list)
# joins S3
S3_LI_acore_list_combined <- combine_acore(fuzz_S3_LI, S3_LI_acore_list)
S3_NL_acore_list_combined <- combine_acore(fuzz_S3_NL, S3_NL_acore_list)

# joins the combined sets (LI and NL)
LI_acore_list_combined <- combine_acore(LI_average, LI_acore_list)
NL_acore_list_combined <- combine_acore(NL_average, NL_acore_list)

# ============== 12. Creates Venn Diagram for Set Overlap ======

# splits accession by ";" delimiter (ie "Q9JHU4-1; Q9JHU4" --> "Q9JHU4-1")
Retina_WP_S1_grouped$Master_Protein_Accessions <- sapply(strsplit(Retina_WP_S1_grouped$Master_Protein_Accessions,";"), `[`, 1)
Retina_WP_S2_grouped$Master_Protein_Accessions <- sapply(strsplit(Retina_WP_S2_grouped$Master_Protein_Accessions,";"), `[`, 1)
Retina_WP_S3_grouped$Master_Protein_Accessions <- sapply(strsplit(Retina_WP_S3_grouped$Master_Protein_Accessions,";"), `[`, 1)


# creates function to prepare venn diagram lists
plot_venn_diag <- function(Set_1, Set_2, Set_3) {
  # Creates list of proteins in each set
  protein_list <- list(
    Set_1 = Set_1$Master_Protein_Accessions, 
    Set_2 = Set_2$Master_Protein_Accessions, 
    Set_3 = Set_3$Master_Protein_Accessions
  )
  # Sets names to protein list
  names(protein_list) <- c("Set 1","Set 2","Set 3")
  # plots venn diagram
  ggvenn(
    protein_list, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, 
    text_size = 2,
    set_name_size = 3
  )
}

# plots venn diagram
plot_venn_diag(Retina_WP_S1_grouped, 
               Retina_WP_S2_grouped, 
               Retina_WP_S3_grouped)

# exports venn diagram
ggsave(
  "Phospho_Protein_Venn.png",
  plot = last_plot(),
  bg = 'white',
  width = 5, 
  height = 5
)
