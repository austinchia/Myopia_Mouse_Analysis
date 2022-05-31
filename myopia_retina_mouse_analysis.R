library("readxl")

# Set 1
abund_s1 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S1_abundance_ratios"))                                                                            
abund_log_s1 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S1_abundance_log2"))                                                                            
abund_group_s1 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S1_abundance_grouped"))                                                                            

# Set 2
abund_s2 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S2_abundance_ratios"))                                                                            
abund_log_s2 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S2_abundance_log2"))                                                                            
abund_group_s2 <- data.frame(read_excel("Myopia_retina_Whole Proteome results_3 sets.xlsx",sheet="Retina_WP_S2_abundance_grouped"))                                                                            

library("MetaboAnalystR")
# initializing object
mSet<-InitDataObjects("conc", "stat", FALSE)
# loading in data
mSet<-Read.TextData(mSet, "abundance_S1_test.csv", "colu", "disc");
# data check
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);

# filtering features
mSet<-FilterVariable(mSet, "none", "F", 25)

# median normalization, log10 transformation, pareto scaling
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

# plot volcano plot
mSet<-Volcano.Anal(mSet, FALSE, 1.5, 0, F, 0.05, TRUE, "raw")
mSet <-PlotVolcano(mSet, "volcano_1_",1, 0, "png", 72, width=NA)

# save transformed dataset
mSet<-SaveTransformedData(mSet)
