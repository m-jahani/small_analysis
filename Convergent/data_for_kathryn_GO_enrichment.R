library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)

convert <-data.frame(Variable=c(c("OM","P1","P2","BICARB","K","MG","CA","Sodium",'PH',"CEC",
                                  "PERCENT_K","PERCENT_MG","PERCENT_CA","PERCENT_NA","SOL_SALTS"),
                                c("latitude","longitude","llevation","MAT","MWMT","MCMT",
                                  "TD","MAP","MSP","AHM","SHM","DD_0",
                                  "DD5","DD_18","DD18","NFFD","bFFP","eFFP",
                                  "FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH")),type=c(rep("soil",15),rep("climate",26)))


add<-c("Baypass","GWAS_corrected","GWAS_uncorrected","Spearman")
files<-c("baypass_nullw_result_recom_bin_FDR","GWAScorrected_nullw_result_recom_bin_FDR","GWASuncorrected_nullw_result_recom_bin_FDR","spearman_nullw_result_recom_bin_FDR")
analy<-c("baypass","GWAS_corrected","GWAS_uncorrected","spearman")
datt<-c("env","gwa","gwa","env")
two_bin <- NULL

for (j in 1:length(files)) {
  read.table(paste0("/moonriseNFS/Null_W_Result/null_w_topcandidate_binning/",add[j],"/",files[j]), header = TRUE,sep="\t") %>%
    mutate(analysis=analy[j],data=datt[j]) %>% 
    rbind(.,two_bin)-> two_bin
}


two_bin %>%
  left_join(.,convert) %>%
  mutate(type=as.character(type)) %>%
  mutate(type=replace_na(type,"phenotype")) %>%
  mutate(comparisonn=paste0(Taxa1,"_",Taxa2)) %>%
  left_join(.,
            data.frame(comparisonn=c("Annuus_Argophyllus",
                                    "Argophyllus_Annuus",
                                    "Annuus_petfal",
                                    "petfal_Annuus",
                                    "Annuus_petpet",
                                    "petpet_Annuus",
                                    "Argophyllus_petfal",
                                    "petfal_Argophyllus",
                                    "Argophyllus_petpet",
                                    "petpet_Argophyllus",
                                    "petfal_petpet",
                                    "petpet_petfal"),
                       comparison_merge=c("Annuus_Argophyllus",
                                          "Annuus_Argophyllus",
                                          "Annuus_petfal",
                                          "Annuus_petfal",
                                          "Annuus_petpet",
                                          "Annuus_petpet",
                                          "Argophyllus_petfal",
                                          "Argophyllus_petfal",
                                          "Argophyllus_petpet",
                                          "Argophyllus_petpet",
                                          "petfal_petpet",
                                          "petfal_petpet"))) %>%
  separate(Window,into=c("chr","S_E"),sep=":",extra="warn",remove = F) %>%
  separate(S_E,into=c("start_window","end_window"),sep="-",extra="warn",remove = T) %>%
  mutate(start_window=as.numeric(start_window)) %>%
  mutate(end_window=as.numeric(end_window)) %>% 
  filter(!Variable %in% c("Guides_3_petals",
           "Guides_3_petals_mod",
           "Guides_individual_petal",
           "Guides_individual_petal_mod")) %>%
  select(data=type,analysis,comparison=comparison_merge,variable=Variable,window=window_name,chr,start_window,end_window,Emp_p) -> two_binn
fwrite(two_binn,"/moonriseNFS/Null_W_Result/files/convergent_2binning_GOenrichmnent_emp_P",sep = "\t")
