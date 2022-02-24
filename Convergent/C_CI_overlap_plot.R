
library(tidyverse)
library(data.table)
library(grid)
library(gridExtra)
options(scipen = 9999)
fread("/Users/mojtabajahani/Downloads/C_CI_overlap_merged_P_LD0.9_1cM",sep = "\t",header = T) %>% 
  mutate(porportion_number_original_cluster_overlap_CI_sig=ifelse(porportion_number_original_cluster_overlap_CI_P >= 0.95 ,
                                                                  "*","")) %>%
  mutate(porportion_length_C_overlap_CI_sig=ifelse(porportion_length_C_overlap_CI_P >= 0.95 ,
                                                   "*","")) %>%
  left_join(.,
data.frame(analysis=c("baypass","spearman","GWAS_corrected","GWAS_uncorrected"),
           group=c("Baypass+Uncorrected","Spearman+Corrected","Spearman+Corrected","Baypass+Uncorrected"))) -> data_C_CI_overlap

P1<-ggplot(data = data_C_CI_overlap, mapping = aes(x = variable,
                               y = direction_C,
                               fill = porportion_length_C_overlap_CI,
                               label=porportion_length_C_overlap_CI_sig)) +
  geom_tile()+
  geom_text()+
  facet_grid(group~data,scales = "free_x",space = "free_x")+
  scale_fill_gradientn(name = "Overlap Length Proportion  ",
                       colors = c("blue3","darkturquoise","green4","darkorange2","red"),
                       breaks = c(0.2,0.4,0.6,0.8,1), labels = c(0.2,0.4,0.6,0.8,1))+
  #guides(fill = guide_legend(title = "Length Proportion",reverse = TRUE))+
  theme(axis.text.x = element_text(size = 8, angle = 90,hjust=1,vjust=0.5),
        #strip.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )+
  xlab(label = "Variables") +
  ylab(label = "Comprisons")

#
P2<-ggplot(data = data_C_CI_overlap, mapping = aes(x = variable,
                                               y = direction_C,
                                               fill = porportion_length_C_overlap_CI_P)) +
  geom_tile()+
  facet_grid(group~data,scales = "free_x",space = "free_x")+
  scale_fill_gradientn(name = "Overlap Length Proportion P"
                       ,colors = c("blue3","darkturquoise","green4","darkorange2","red"),
                       breaks = c(0.2,0.4,0.6,0.8,1), labels = c(0.2,0.4,0.6,0.8,1))+
  #guides(fill = guide_legend(title = "Length Proportion",reverse = TRUE))+
  theme(axis.text.x = element_text(size = 8, angle = 90,hjust=1,vjust=0.5),
        #strip.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )+
  xlab(label = "Variables") +
  ylab(label = "Comprisons")

grid.arrange(P1, P2, nrow = 2)


P3<-ggplot(data = data_C_CI_overlap, mapping = aes(x = variable,
                                                   y = direction_C,
                                                   fill = porportion_number_original_cluster_overlap_CI,
                                                   label=porportion_number_original_cluster_overlap_CI_sig)) +
  geom_tile()+
  geom_text()+
  facet_grid(group~data,scales = "free_x",space = "free_x")+
  scale_fill_gradientn(name = "Overlap Number Proportion  ",
                       colors = c("blue3","darkturquoise","green4","darkorange2","red"),
                       breaks = seq(min(data_C_CI_overlap$porportion_number_original_cluster_overlap_CI),
                                    max(data_C_CI_overlap$porportion_number_original_cluster_overlap_CI),
                                    (max(data_C_CI_overlap$porportion_number_original_cluster_overlap_CI)-min(data_C_CI_overlap$porportion_number_original_cluster_overlap_CI))/4), labels = c(0.2,0.4,0.6,0.8,1))+
  #guides(fill = guide_legend(title = "Length Proportion",reverse = TRUE))+
  theme(axis.text.x = element_text(size = 8, angle = 90,hjust=1,vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )+
  xlab(label = "Variables") +
  ylab(label = "Comprisons")

#
P4<-ggplot(data = data_C_CI_overlap, mapping = aes(x = variable,
                                                   y = direction_C,
                                                   fill = porportion_number_original_cluster_overlap_CI_P)) +
  geom_tile()+
  facet_grid(group~data,scales = "free_x",space = "free_x")+
  scale_fill_gradientn(name = "Overlap Number Proportion P"
                       ,colors = c("blue3","darkturquoise","green4","darkorange2","red"),
                       breaks = c(0.2,0.4,0.6,0.8,1), labels = c(0.2,0.4,0.6,0.8,1))+
  #guides(fill = guide_legend(title = "Length Proportion",reverse = TRUE))+
  theme(axis.text.x = element_text(size = 8, angle = 90,hjust=1,vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )+
  xlab(label = "Variables") +
  ylab(label = "Comprisons")

grid.arrange(P3, P4, nrow = 2)


range(as.numeric(data_C_CI_overlap$porportion_number_original_cluster_overlap_CI))
data_C_CI_overlap %>% glimpse
