library(data.table)
library(tidyverse)

fread("/DATA/home/mjahani/CONVERG/all_id_gene_window_number_maf.03") %>%
  select(-Gene)  %>%
  group_by(Taxa,window,window5) %>%
  distinct(ID)  %>%
  ungroup() -> WIN_SNP_NEW
 

data.frame(V1=c("Annuus","Argophyllus"),
           V2=c("Annuus","petpet"),
           V3=c("Annuus","petfal"),
           V4=c("Argophyllus","petpet"),
           V5=c("Argophyllus","petfal"),
           V6=c("petpet","petfal")) -> TAXA1_TAXA2
           
  WIN_SNP_NEW %>% 
    group_by(Taxa,window5,window) %>% 
    tally() %>%
    ungroup() %>%
    spread(key = Taxa, value = n) -> WIN_SNP_COUNT

  
  for (j in 1:ncol(TAXA1_TAXA2)) {
    WIN_SNP_NEW %>%
      filter(Taxa %in% pull(TAXA1_TAXA2,paste0("V",j))) %>%
      select(-Taxa) %>% 
      group_by(window5,window,ID) %>% 
      tally() %>%
      ungroup() %>%
      filter(n == 2) %>%
      group_by(window5,window) %>%
      tally() %>%
      ungroup() %>%
      rename(!!paste0(pull(TAXA1_TAXA2,paste0("V",j))[1],"_",pull(TAXA1_TAXA2,paste0("V",j))[2]) := n) %>%
      right_join(.,WIN_SNP_COUNT) -> WIN_SNP_COUNT
  }
  
  WIN_SNP_COUNT %>%
    select(window5,
           window,
           Annuus,
           Argophyllus,
           petpet,
           petfal,
           Annuus_Argophyllus,
           Annuus_petpet,
           Annuus_petfal,
           Argophyllus_petpet,
           Argophyllus_petfal,
           petpet_petfal) -> WIN_SNP_COUNT
  
  WIN_SNP_COUNT[is.na(WIN_SNP_COUNT)] <- 0
  
  fwrite(WIN_SNP_COUNT,"/DATA/home/mjahani/CONVERG/SNP_count.csv",
         sep = ",",
         quote = F,
         col.names = T)
  