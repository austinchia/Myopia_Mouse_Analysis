# About
# This script is used for Whole Protein
# This script reads in raw data > preprocesses data > to form a data matrix
# Data matrix is used for further stats analysis

# ============== Clears Memory ======
# clears all objects includes hidden objects
rm(list = ls(all.names = TRUE)) 

# frees up memory and reports the memory usage.
gc() 

# ============== Loads Packages =======
library(readxl)
library(dplyr)
library(data.table)
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


# ====== A) Prepares Data Matrix for Abundance Ratio =====

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

# ============== 2. Manipulates Data  ===============
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

# ============== 3. Combines 3 Sets And Cleans Data ======

# combines all 3 sets into 1 dataset 
# (left joins to S3 because it has most no of proteins)
ratio_combined <- left_join(Retina_WP_S3, Retina_WP_S2, by = 'Accession') %>%
  left_join(Retina_WP_S1, by = 'Accession') %>%
  na.omit()

# convert blanks to NA
ratio_combined <- ratio_combined %>% 
  mutate_all(na_if,"")

# splits accession number (ie Q9JHU4-1)
ratio_combined$Accession <- sapply(strsplit(ratio_combined$Accession,"-"), `[`, 1)

# exports accession numbers to upload to Uniprot
fwrite(data.frame(ratio_combined$Accession), "Whole_Protein_Accession.csv", sep = ",")

# ============== 4. Combines Uniprot Data To Combined Matrix ======
# reads in Gene Symbol table downloaded from Uniprot
gene_symbol <- fread("Whole_Protein_Accession_Map.csv",sep=',')

# splits gene symbol by break
gene_symbol_map <- data.frame(str_split_fixed(gene_symbol$`From	To`, '\t',2))
colnames(gene_symbol_map) <- c("Accession", "Gene Symbol") 

# merges gene symbol column to main df
ratio_combined_no_na <- left_join(ratio_combined, 
                                  gene_symbol_map, 
                                  by="Accession") %>%
  relocate(`Gene Symbol`, .after = `Accession`) %>%
  na.omit() %>%
  
  # adds number to the end of duplicate gene symbols (ie Sptbn1-2)
  group_by(`Gene Symbol`) %>%
  mutate(`GS_count` = 1:n()) %>%
  mutate(`Gene Symbol` = ifelse(`GS_count` == 1, 
                                `Gene Symbol`, 
                                paste0(`Gene Symbol`, "-", `GS_count`))) %>%
  select(-`GS_count`)
  
# exports combined abundance ratio matrix to csv
fwrite(ratio_combined_no_na, "abund_ratio_combined_GS.csv", sep = ",")


# ===================== B) Using Grouped Abundance



# ====== B) Prepares Data Matrix for Grouped Abundance =====

# ============== 1. Reads Raw Data ===================

# reads S1 raw data
Retina_WP_S1_grouped <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S1', na = c("", "NA")) %>%
  select(c(`Accession`,`Abundances (Grouped)`)) %>%
  na.omit()

# reads S2 raw data
Retina_WP_S2_grouped <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S2', na = c("", "NA")) %>%
  select(c(`Accession`,`Abundances (Grouped)`)) %>%
  na.omit()

# reads S3 raw data
Retina_WP_S3_grouped <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S3', na = c("", "NA")) %>%
  select(c(`Accession`,`Abundances (Grouped)`)) %>%
  na.omit()

# ============== 2. Manipulates Data  ===============
# S1 Data manipulation
{
  # splits string into columns
  Retina_WP_S1_grouped_split <- str_split_fixed(as.character(Retina_WP_S1_grouped$`Abundances (Grouped)`), ';',15)
  
  # adds columns to original data
  Retina_WP_S1_grouped <- cbind(Retina_WP_S1_grouped, Retina_WP_S1_grouped_split)
  colnames(Retina_WP_S1_grouped) <- c('Accession', 'Abundances (Grouped)', 'S1_LI_1hr','S1_LI_6hr','S1_LI_9hr','S1_LI_D1','S1_LI_D14','S1_LI_D3','S1_LI_D7','S1_NL_0hr','S1_NL_1hr','S1_NL_6hr','S1_NL_9hr','S1_NL_D1','S1_NL_D14','S1_NL_D3','S1_NL_D7')
  
  # removes abundance column
  Retina_WP_S1_grouped <- select(Retina_WP_S1_grouped, -c(`Abundances (Grouped)`)) %>%
    mutate(S1_LI_0hr = S1_NL_0hr) %>%
    relocate(S1_LI_0hr, .after = `Accession`) %>%
    relocate(S1_LI_D14, .after = `S1_LI_D7`) %>%
    relocate(S1_NL_D14, .after = `S1_NL_D7`)
  
}

# S2 Data manipulation
{
  # splits string into columns
  Retina_WP_S2_grouped_split <- str_split_fixed(as.character(Retina_WP_S2_grouped$`Abundances (Grouped)`), ';',15)
  
  # adds columns to original data
  Retina_WP_S2_grouped <- cbind(Retina_WP_S2_grouped, Retina_WP_S2_grouped_split)
  colnames(Retina_WP_S2_grouped) <- c('Accession', 'Abundances (Grouped)', 'S2_LI_1hr','S2_LI_6hr','S2_LI_9hr','S2_LI_D1','S2_LI_D14','S2_LI_D3','S2_LI_D7','S2_NL_0hr','S2_NL_1hr','S2_NL_6hr','S2_NL_9hr','S2_NL_D1','S2_NL_D14','S2_NL_D3','S2_NL_D7')
  
  # removes abundance column
  Retina_WP_S2_grouped <- select(Retina_WP_S2_grouped, -c(`Abundances (Grouped)`)) %>%
    mutate(S2_LI_0hr = S2_NL_0hr) %>%
    relocate(S2_LI_0hr, .after = `Accession`) %>%
    relocate(S2_LI_D14, .after = `S2_LI_D7`) %>%
    relocate(S2_NL_D14, .after = `S2_NL_D7`)
  
}

# S3 Data manipulation
{
  # splits string into columns
  Retina_WP_S3_grouped_split <- str_split_fixed(as.character(Retina_WP_S3_grouped$`Abundances (Grouped)`), ';',15)
  
  # adds columns to original data
  Retina_WP_S3_grouped <- cbind(Retina_WP_S3_grouped, Retina_WP_S3_grouped_split)
  colnames(Retina_WP_S3_grouped) <- c('Accession', 'Abundances (Grouped)', 'S3_LI_1hr','S3_LI_6hr','S3_LI_9hr','S3_LI_D1','S3_LI_D14','S3_LI_D3','S3_LI_D7','S3_NL_0hr','S3_NL_1hr','S3_NL_6hr','S3_NL_9hr','S3_NL_D1','S3_NL_D14','S3_NL_D3','S3_NL_D7')
  
  # removes abundance column
  Retina_WP_S3_grouped <- select(Retina_WP_S3_grouped, -c(`Abundances (Grouped)`)) %>%
    mutate(S3_LI_0hr = S3_NL_0hr) %>%
    relocate(S3_LI_0hr, .after = `Accession`) %>%
    relocate(S3_LI_D14, .after = `S3_LI_D7`) %>%
    relocate(S3_NL_D14, .after = `S3_NL_D7`)
  
}

# ============== 3. Selects Columns From Main Grouped Matrix =========
grouped_combined_GS <- fread("grouped_combined_GS_accounted.csv",sep=',')

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
  mutate('NL_0hr_mean' = rowMeans(subset(., select = c(`S1_NL_0hr`,`S2_NL_0hr`,`S3_NL_0hr`)))) %>%
  mutate('NL_1hr_mean' = rowMeans(subset(., select = c(`S1_NL_1hr`,`S2_NL_1hr`,`S3_NL_1hr`)))) %>%
  mutate('NL_6hr_mean' = rowMeans(subset(., select = c(`S1_NL_6hr`,`S2_NL_6hr`,`S3_NL_6hr`)))) %>%
  mutate('NL_9hr_mean' = rowMeans(subset(., select = c(`S1_NL_9hr`,`S2_NL_9hr`,`S3_NL_9hr`)))) %>%
  mutate('NL_D1_mean' = rowMeans(subset(., select = c(`S1_NL_D1`,`S2_NL_D1`,`S3_NL_D1`)))) %>%
  mutate('NL_D3_mean' = rowMeans(subset(., select = c(`S1_NL_D3`,`S2_NL_D3`,`S3_NL_D3`)))) %>%
  mutate('NL_D7_mean' = rowMeans(subset(., select = c(`S1_NL_D7`,`S2_NL_D7`,`S3_NL_D7`)))) %>%
  mutate('NL_D14_mean' = rowMeans(subset(., select = c(`S1_NL_D14`,`S2_NL_D14`,`S3_NL_D14`)))) %>%
  select(`NL_0hr_mean`, `NL_1hr_mean`, `NL_6hr_mean`, `NL_9hr_mean`, `NL_D1_mean`, `NL_D3_mean`, `NL_D7_mean`, `NL_D14_mean`)

# ============== 4. Creates Timepoints and Binds to Original Dataframe ====

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

# ============== 5. Scales Data ========================

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

# ============== 6. Estimates Fuzzifier (ie m1) ================

m1_S1_LI <- mestimate(S1_LI_eSet)
m1_S1_NL <- mestimate(S1_NL_eSet)
m1_S2_LI <- mestimate(S2_LI_eSet)
m1_S2_NL <- mestimate(S2_NL_eSet)
m1_S3_LI <- mestimate(S3_LI_eSet)
m1_S3_NL <- mestimate(S3_NL_eSet)

# estimates fuzzifier for LI and NL (all sets)
m1_average_LI <- mestimate(LI_average_eSet)
m1_average_NL <- mestimate(NL_average_eSet)

# plots scree plot - determine no of centroids
# Dmin(LI_average_eSet, m=m1, crange=seq(5,20,1), repeats=3, visu=TRUE)

# ============== 7. Plots Mfuzz Plots ======================

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
  # png(file = "mfuzz_S1_LI.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S1_LI_eSet, m1_S1_LI)
  # dev.off()
  # 
  # png(file="mfuzz_S1_NL.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S1_NL_eSet, m1_S1_NL)
  # dev.off()
  # 
  # # plots mfuzz for Set 2
  # png(file="mfuzz_S2_LI.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S2_LI_eSet, m1_S2_LI)
  # dev.off()
  # 
  # png(file="mfuzz_S2_NL.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S2_NL_eSet, m1_S2_NL)
  # dev.off()
  # 
  # # plots mfuzz for Set 3
  # png(file="mfuzz_S3_LI.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S3_LI_eSet, m1_S3_LI)
  # dev.off()
  # 
  # png(file="mfuzz_S3_NL.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(S3_NL_eSet, m1_S3_NL)
  # dev.off()
  # 
  # # plots mfuzz for LI (all sets)
  # png(file="mfuzz_LI_average.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(LI_average_eSet, m1_average_LI)
  # dev.off()
  # 
  # # plots mfuzz for NL (all sets)
  # png(file="mfuzz_NL_average.png",
  #     width = 1000,
  #     height = 1000,)
  # plot_mfuzz(NL_average_eSet, m1_average_NL)
  # dev.off()
}

# ============== 8. Validates and Evaulates Mfuzz Model ==========

# creates correlation matrix between cluster centroids
# (no more than 0.85)
# correlation_matrix <- data.frame(cor(t(cl[[1]])))

# ============== 9. Extracts Gene Lists From Clusters ==========
# creates function to extract genes (acore list) in each cluster
get_genes <- function(x, m1) {
  cl <- mfuzz(x, c = 12,m = m1)
  acore_x <- acore(x,cl,min.acore=0)
  do.call(rbind, lapply(seq_along(acore_x), 
                        function(i){ data.frame(Cluster=i, 
                                                acore_x[[i]])}))
}

# uses get_genes function to extract acore list
S1_LI_acore_list <- get_genes(S1_LI_eSet, m1_S1_LI)
S1_NL_acore_list <- get_genes(S1_NL_eSet, m1_S1_NL)

S2_LI_acore_list <- get_genes(S2_LI_eSet, m1_S2_LI)
S2_NL_acore_list <- get_genes(S2_NL_eSet, m1_S2_NL)

S3_LI_acore_list <- get_genes(S3_LI_eSet, m1_S3_LI)
S3_NL_acore_list <- get_genes(S3_NL_eSet, m1_S3_NL)

# extracts acore list for combined sets (LI)
LI_acore_list <- LI_average_eSet %>%
  get_genes(., m1_average_LI) %>%
  na.omit() %>%
  rename("Gene Symbol" = "NAME")  

# extracts acore list for combined sets (NL)
NL_acore_list <- NL_average_eSet %>%
  get_genes(., m1_average_NL) %>%
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
  names(acore_combined)[names(acore_combined) == 'rowname'] <- 'NAME'
  return(acore_combined)
}

# creates function to join list with abundance (only for average combined sets)
combine_acore_GS <- function(abundance_df, acore_list) {
  acore_combined <- left_join(rownames_to_column(abundance_df), 
                              acore_list, 
                              by=c("rowname" = "Gene Symbol"))
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
LI_acore_list_combined <- combine_acore_GS(LI_average, LI_acore_list)
NL_acore_list_combined <- combine_acore_GS(NL_average, NL_acore_list)

# ============== 10. Creates Venn Diagram for Set Overlap ======

# == selects raw data for venn (identified) ==
{
  # reads S1 raw data
  Retina_WP_S1_identified <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S1', na = c("", "NA")) %>%
    select(c(`Accession`,`Abundances (Grouped)`))
  
  # reads S2 raw data
  Retina_WP_S2_identified <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S2', na = c("", "NA")) %>%
    select(c(`Accession`,`Abundances (Grouped)`))
  
  # reads S3 raw data
  Retina_WP_S3_identified <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S3', na = c("", "NA")) %>%
    select(c(`Accession`,`Abundances (Grouped)`))
}

# creates function to plot venn diagram (proteins) (using Accession No)
plot_venn_protein <- function(Set_1, Set_2, Set_3) {
  # Creates list of proteins in each set
  protein_list <- list(
    Set_1 = Set_1$Accession, 
    Set_2 = Set_2$Accession, 
    Set_3 = Set_3$Accession
  )
  # Sets names to protein list
  names(protein_list) <- c("Set 1","Set 2","Set 3")
  # plots venn diagram
  ggvenn(
    protein_list, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, 
    text_size = 3,
    set_name_size = 3
  )
}

# # plots venn diagram (quantifiable, protein)
# {
#   # plots venn
#   plot_venn_protein(Retina_WP_S1_grouped, 
#                     Retina_WP_S2_grouped, 
#                     Retina_WP_S3_grouped
#   )
#   
#   # exports venn diagrams (quantifiable, protein)
#   ggsave(
#     "Whole_Protein_Venn_Protein_Quantifiable.png",
#     plot = last_plot(),
#     bg = 'white',
#     width = 5, 
#     height = 5
#   )
# }
# 
# # plots venn diagram (identified, protein)
# {
#   # plots venn
#   plot_venn_protein(Retina_WP_S1_identified, 
#                     Retina_WP_S2_identified, 
#                     Retina_WP_S3_identified
#   )
#   
#   # exports venn diagrams (identified proteins, whole protein)
#   ggsave(
#     "Whole_Protein_Venn_Protein_Identified.png",
#     plot = last_plot(),
#     bg = 'white',
#     width = 5, 
#     height = 5
#   )
#   
# }



# ============== 11. Exports Protein Lists from Mfuzz Clusters (LI) ====
# cleans main matrix before combining
grouped_combined_GS_export <- grouped_combined_GS %>%
  select(`Gene Symbol`, `Accession`, `S1_LI_0hr`,	`S1_LI_1hr`,	`S1_LI_6hr`,	`S1_LI_9hr`,	`S1_LI_D1`,	`S1_LI_D3`,	`S1_LI_D7`,	`S1_LI_D14`) %>%
  # replaces 0 with NA
  na_if(0) %>%
  # removes NAs
  na.omit()

# merges acore list with main matrix
LI_acore_list_data <- LI_acore_list_combined %>%
  # merges list with LI data
  left_join(.,
            grouped_combined_GS_export,
            by= 'Gene Symbol') %>%
  
  # merges with NL data
  left_join(.,
            NL_acore_list_combined,
            by= 'Gene Symbol') %>%
  
  # calculates means for LI & NL + log2FC
  mutate(., LI_mean = rowMeans(select(.,
                                      LI_0hr_mean:LI_D14_mean), na.rm = TRUE)) %>%
  mutate(., NL_mean = rowMeans(select(.,
                                      NL_0hr_mean:NL_D14_mean), na.rm = TRUE)) %>%
  mutate(., log2FC = log2(LI_mean/NL_mean)) %>%
  
  # selects useful columns
  select(`Accession`, `Cluster.x`, `log2FC`)


# creates function to filter out by cluster
filter_cluster <- function(dataframe, cluster_no) {
  dataframe %>%
    filter(`Cluster.x` == cluster_no) %>%
    select(-c(`Cluster.x`))
} 

# filters acore list by cluster
{
  LI_acore_list_cl_1 <- filter_cluster(LI_acore_list_data, 1)
  LI_acore_list_cl_4 <- filter_cluster(LI_acore_list_data, 4)
  LI_acore_list_cl_5 <- filter_cluster(LI_acore_list_data, 5)
  LI_acore_list_cl_7 <- filter_cluster(LI_acore_list_data, 7)
  LI_acore_list_cl_8 <- filter_cluster(LI_acore_list_data, 8)
  LI_acore_list_cl_10 <- filter_cluster(LI_acore_list_data, 10)
  LI_acore_list_cl_11 <- filter_cluster(LI_acore_list_data, 11)
  LI_acore_list_cl_12 <- filter_cluster(LI_acore_list_data, 12)
  
}

# exports acore list for clusters 1, 4, 5, 7, 8, 10, 11, 12
{
  fwrite(LI_acore_list_cl_1, "C:/Users/austi/Documents/Data Projects/R Projects/SERIDataAnalysis/Myopia_Mouse_Analysis/Output/Whole_Protein_LI_Cluster_1.csv", sep = ",")
  fwrite(LI_acore_list_cl_4, "C:/Users/austi/Documents/Data Projects/R Projects/SERIDataAnalysis/Myopia_Mouse_Analysis/Output/Whole_Protein_LI_Cluster_4.csv", sep = ",")
  fwrite(LI_acore_list_cl_5, "C:/Users/austi/Documents/Data Projects/R Projects/SERIDataAnalysis/Myopia_Mouse_Analysis/Output/Whole_Protein_LI_Cluster_5.csv", sep = ",")
  fwrite(LI_acore_list_cl_7, "C:/Users/austi/Documents/Data Projects/R Projects/SERIDataAnalysis/Myopia_Mouse_Analysis/Output/Whole_Protein_LI_Cluster_7.csv", sep = ",")
  fwrite(LI_acore_list_cl_8, "C:/Users/austi/Documents/Data Projects/R Projects/SERIDataAnalysis/Myopia_Mouse_Analysis/Output/Whole_Protein_LI_Cluster_8.csv", sep = ",")
  fwrite(LI_acore_list_cl_10, "C:/Users/austi/Documents/Data Projects/R Projects/SERIDataAnalysis/Myopia_Mouse_Analysis/Output/Whole_Protein_LI_Cluster_10.csv", sep = ",")
  fwrite(LI_acore_list_cl_11, "C:/Users/austi/Documents/Data Projects/R Projects/SERIDataAnalysis/Myopia_Mouse_Analysis/Output/Whole_Protein_LI_Cluster_11.csv", sep = ",")
  fwrite(LI_acore_list_cl_12, "C:/Users/austi/Documents/Data Projects/R Projects/SERIDataAnalysis/Myopia_Mouse_Analysis/Output/Whole_Protein_LI_Cluster_12.csv", sep = ",")
}



# ============== 12. Runs Gene Set Variation Analysis (GSVA) ====
library(GSVA)
library(GSVAdata)
library(org.Mm.eg.db)
library(limma)

# selects D7 data for GSVA
gsva_D7_mat <- grouped_combined_GS %>%
  # selects needed columns
  dplyr::select(`Gene Symbol`,	`S1_LI_D7`,	`S2_LI_D7`, `S3_LI_D7`,`S1_NL_D7`,	`S2_NL_D7`, `S3_NL_D7`) %>%
  
  # replaces 0 with NA
  na_if(0) %>%
  
  # removes NAs
  na.omit() %>%
  
  # set rownames as `Gene Symbol`
  column_to_rownames(., var = "Gene Symbol") %>%
  
  # converts dataframe to matrix
  data.matrix()

# gets gene symbols and entrez ID from mouse GO database
go_annot <- AnnotationDbi::select(org.Mm.eg.db, keys=keys(org.Mm.eg.db), columns = c("SYMBOL"))

# # splits by entrez ID and GO and converts to list
# genes_by_symbol <- split(go_annot$`SYMBOL`,go_annot$ENTREZID)
# 
# mouse_gene_set <- data.frame(unlist(readRDS("Mm.c2.cp.kegg.v7.1.entrez.rds"))) %>%
#   rownames_to_column(., var = "Pathway") %>%
#   relocate(`Pathway`, .after = last_col()) %>%
#   rename("ENTREZID" = "unlist.readRDS..Mm.c2.cp.kegg.v7.1.entrez.rds...")
# 
# mouse_set_combined <- left_join(go_annot, mouse_gene_set, "ENTREZID") %>%
#   na.omit() %>%
#   distinct(`ENTREZID`, .keep_all= TRUE)

# calculates GSVA enrichment score estimates
gsva_estimate <- gsva(gsva_D7_mat, genes_by_symbol)

mod <- model.matrix(~ factor(c("LI", "LI", "LI", "NL", "NL", "NL")))
colnames(mod) <- c("LI", "NL")
fit <- lmFit(gsva_estimate, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.05)
summary(res)
# LI   NL
# Down   5633    0
# NotSig  149 5784
# Up        2    0

# plots volcano plot
tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.05]
plot(tt$logFC, -(tt$P.Value), pch=".", cex=4, col=grey(0.75),
     main="", xlab="GSVA enrichment score difference", ylab=expression(-log[10]~~Raw~P-value))
abline(h=max(tt$P.Value[tt$adj.P.Val <= 0.01]), col=grey(0.5), lwd=1, lty=2)
points(tt$logFC[match(DEpwys, rownames(tt))],
       tt$P.Value[match(DEpwys, rownames(tt))], pch=".", cex=5, col="darkred")
text(max(tt$logFC)*0.85, max(tt$P.Value[tt$adj.P.Val <= 0.01]), "1% FDR", pos=3)

# plots heatmap
DEpwys_es <- DEpwys
colorLegend <- c("darkred", "darkblue")
names(colorLegend) <- c("LI", "NL")
sample.color.map <- colorLegend[(gsva_D7_mat)]
names(sample.color.map) <- colnames(gsva_D7_mat)
sampleClustering <- hclust(as.dist(1-cor(gsva_D7_mat, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(gsva_D7_mat), method="pearson")),
                            method="complete")
heatmap(gsva_D7_mat, ColSideColors=sample.color.map, xlab="samples",
        ylab="Pathways", margins=c(2, 20),
        labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "",
                                          rownames(gsva_D7_mat))), 1, 35),
        labCol="", scale="row", Colv=as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering))
legend("topleft", names(colorLegend), fill=colorLegend, inset=0.01, bg="white")

# ============== sample =========
data(leukemia)
leukemia_eset
data(c2BroadSets)

leu_mat <- data.matrix(leukemia_eset)

leukemia_es <- gsva(leukemia_eset, c2BroadSets, min.sz=10, max.sz=500)
leukemia_es$subtype

mod <- model.matrix(~ factor(leukemia_es$subtype))
colnames(mod) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.01)
summary(res)

tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]
plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75),
     main="", xlab="GSVA enrichment score difference", ylab=expression(-log[10]~~Raw~P-value))
abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])), col=grey(0.5), lwd=1, lty=2)
points(tt$logFC[match(DEpwys, rownames(tt))],
       -log10(tt$P.Value[match(DEpwys, rownames(tt))]), pch=".", cex=5, col="darkred")
text(max(tt$logFC)*0.85, -log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])), "1% FDR", pos=3)

DEpwys_es <- exprs(leukemia_es[DEpwys, ])
colorLegend <- c("darkred", "darkblue")
names(colorLegend) <- c("ALL", "MLL")
sample.color.map <- colorLegend[pData(leukemia_es)[, "subtype"]]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                            method="complete")
heatmap(DEpwys_es, ColSideColors=sample.color.map, xlab="samples",
        ylab="Pathways", margins=c(2, 20),
        labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "",
                                          rownames(DEpwys_es))), 1, 35),
        labCol="", scale="row", Colv=as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering))
legend("topleft", names(colorLegend), fill=colorLegend, inset=0.01, bg="white")
