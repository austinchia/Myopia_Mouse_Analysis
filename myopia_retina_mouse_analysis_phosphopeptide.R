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
ratio_combined <- left_join(Retina_WP_S2_grouped, Retina_WP_S1_grouped, by = 'Annotated_Sequence') %>%
  left_join(Retina_WP_S3_grouped, by = 'Annotated_Sequence') %>%
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
  
  # adds number to the end of duplicate phosphopeptide (ie Sptbn1-2)
  group_by(`Annotated_Sequence`) %>%
  mutate(`seq_count` = 1:n()) %>%
  mutate(`Annotated_Sequence` = ifelse(`seq_count` == 1, 
                                `Annotated_Sequence`, 
                                paste0(`Annotated_Sequence`, "-", `seq_count`))) %>%
  
  # removes unused columns
  select(-c(`GS_count`,
            `Modifications.x`,
            `Master_Protein_Accessions.x`,
            `Modifications.y`,
            `Master_Protein_Accessions.y`
            ))

# exports combined grouped abundance matrix to csv
fwrite(ratio_combined_no_na, "phospho_abund_grouped_combined.csv", sep = ",")


