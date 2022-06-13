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
Retina_WP_S1_grouped <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S1', na = c("", "NA")) %>%
  select(c(`Accession`,`Abundances (Grouped)`))

# reads S2 raw data
Retina_WP_S2_grouped <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S2', na = c("", "NA")) %>%
  select(c(`Accession`,`Abundances (Grouped)`))

# reads S3 raw data
Retina_WP_S3_grouped <- read_excel('Myopia_retina_Whole Proteome results_3 sets.xlsx', sheet = 'Retina_WP_S3', na = c("", "NA")) %>%
  select(c(`Accession`,`Abundances (Grouped)`))