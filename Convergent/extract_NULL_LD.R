library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)


address <- "/home/shsoudi/scratch/output/split_files/"
address_save <- "/home/shsoudi/scratch/output/final_NULL/"

args = commandArgs(trailingOnly = TRUE)

info <- expand.grid(
  data   = c("env"),
  Taxa   = c("Annuus", "Argophyllus", "petpet", "petfal"),
  chr  =  c(paste0("Ha412HOChr0",seq(1,9,1)),paste0("Ha412HOChr",seq(10,17,1)))
) %>% 
  rbind(.,
        data.frame(
          data   = c("gwa"),
          Taxa   = c("Annuus"),
          chr  =  c(paste0("Ha412HOChr0",seq(1,9,1)),paste0("Ha412HOChr",seq(10,17,1))))
  ) %>%
  slice(as.numeric(args[1]):as.numeric(args[1]))
        
list.files(path = address) -> file_list

registerDoParallel(cores=32)
mid_data <- foreach(i=1:length(file_list), .combine ='rbind', .errorhandling='stop') %dopar% {
  fread(paste0(address,file_list[i]),
        header = F,showProgress=F) %>%
    rename(data=V1,Taxa=V2,window_A=V3,window_B=V4,LD=V5,snp_count=V6) %>%
    filter(data==as.character(info[1,1])) %>%
    filter(Taxa==as.character(info[1,2])) %>%
    filter(grepl(as.character(info[1,3]),window_A)) }

mid_data %>%
  distinct(data,Taxa,window_A,window_B,.keep_all = T) -> data
rm(mid_data)

if (as.character(info[1,1]=="gwa")) {
   fwrite(data,file = paste0(address_save,"GWA_Annuus"),col.names = F , append = T , sep = "\t")
   fwrite(info,file = paste0(address_save,"GWA_Annuus_info"),col.names = F , append = T , sep = "\t")
} else {
  if (as.character(info[1,2]=="Annuus")) {
    fwrite(data,file = paste0(address_save,"ENV_Annuus"),col.names = F , append = T , sep = "\t")
    fwrite(info,file = paste0(address_save,"ENV_Annuus_info"),col.names = F , append = T , sep = "\t")
    } else if (as.character(info[1,2])=="Argophyllus") {
    fwrite(data,file = paste0(address_save,"ENV_Argophyllus"),col.names = F , append = T , sep = "\t")
    fwrite(info,file = paste0(address_save,"ENV_Argophyllus_info"),col.names = F , append = T , sep = "\t")
    } else if (as.character(info[1,2])=="petpet") {
    fwrite(data,file = paste0(address_save,"ENV_petpet"),col.names = F , append = T , sep = "\t")
    fwrite(info,file = paste0(address_save,"ENV_petpet_info"),col.names = F , append = T , sep = "\t")
    } else if (as.character(info[1,2])=="petfal") {
    fwrite(data,file = paste0(address_save,"ENV_petfal"),col.names = F , append = T , sep = "\t")
    fwrite(info,file = paste0(address_save,"ENV_petfal_info"),col.names = F , append = T , sep = "\t")}
  }
  