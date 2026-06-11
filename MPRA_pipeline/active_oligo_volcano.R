#25.08.2024
#This script uses the humanMPRA quantitative data in comb_df_combined_fdr_small.csv to create volcano plot

library(ggplot2)
library(dplyr)
library(ggExtra)


setwd('/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/comparative_analysis_combined')

comp_data = read_csv('humanMPRA_with_seq_final2.csv',col_names = TRUE)

head(comp_data)


filt_data = comp_data%>%
  filter(DNA_counts_raw_ancestral>5)%>%
  filter(DNA_counts_raw_derived>5)%>%
  ggplot(aes(x = 1, y = log10(DNA_counts_raw_ancestral)))+
           geom_violin()
  
test = visual_data%>%
  filter(log_fdr==0)%>%
  filter(logFC_derived_vs_ancestral>1)

visual_data = comp_data %>%
  filter(!is.na(differential_activity))%>%
  filter(DNA_counts_raw_ancestral > 5)%>%
  filter(DNA_counts_raw_derived > 5)%>%
  mutate(libraryNum = substr(oligo, nchar(oligo) - 1, nchar(oligo)))%>%
  mutate(log_fdr = -log10(differential_activity_fdr))%>%
  mutate(effect = case_when(differential_activity_fdr >= 0.050 | abs(logFC_derived_vs_ancestral)< 0.263 ~ "no difference",
                            differential_activity_fdr < 0.05 & logFC_derived_vs_ancestral>0.263 ~ "increase",
                            differential_activity_fdr < 0.05 & logFC_derived_vs_ancestral<0.263 ~ "decrease"))

min(visual_data$logFC_derived_vs_ancestral,na.rm = TRUE)

  
ggplot(data = visual_data,aes(x = logFC_derived_vs_ancestral,y = log_fdr))+
  #geom_point(aes(color = effect))+  
  geom_point(aes(color = log(DNA_counts_raw_derived)))+  
  labs(y = "-log(Differential activity FDR)", x ="LFC") +
  geom_vline(xintercept = -0.263, linetype = "dashed", color = "blue") +  # vertical line at x = 0
  geom_vline(xintercept = 0.263, linetype = "dashed", color = "blue") +  # vertical line at x = 0
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +  # horizontal line at y = 0.05
  scale_y_continuous(limits = c(0, 10))+
  theme_minimal()


ggplot(data = visual_data,aes(x = DNA_counts_raw_derived,y = log_fdr))+
  #geom_point(aes(color = effect))+  
  geom_point()+  
  scale_y_continuous(limits = c(0, 10))+
  theme_minimal()


  