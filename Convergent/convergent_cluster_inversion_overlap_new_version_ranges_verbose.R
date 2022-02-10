library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(IRanges)

address <- "/data/users/mjahani/convergent_cluster_inversion_overlap_new_version/"
fread("/data/users/mjahani/inversion.bed",sep = "\t",header = T) -> bed
convert<-data.frame(species=c("Annuus","Argophyllus","petpet","petfal"),
                    spe=c("ann","arg","pet","fal"))
bed %>%
  mutate(sv_new=gsub("ann","",sv_name)) %>%
  mutate(sv_new=gsub("arg","",sv_new)) %>%
  mutate(sv_new=gsub("pet","",sv_new)) %>%
  left_join(.,convert) %>%
  mutate(sv_name_new=paste0(spe,sv_new)) %>%
  select(chr,start,end,species,sv_name,sv_name_new)-> bed

fread(paste0("/data/users/mjahani/LD_cluster/cluster/convergent_clustering_0.9_1_cM_summary"),sep = "\t",header = T) %>%
  left_join(.,
            data.frame(direction=c("Annuus_Argophyllus",
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
                       comparison=c("Annuus_Argophyllus",
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
  distinct(data,analysis,variable,chromosome,range,size,N_convergent,direction,comparison) %>%
  separate(comparison, into = c("Taxa1","Taxa2"),sep = "_",remove = F) %>%
  separate(range,into = c("start","end"),sep = ":",remove = F) %>%
  mutate(start=as.numeric(start),end=as.numeric(end))-> convergent_cluster


#overlap between inversion and conversion regions(inclusing partial overlap)
registerDoParallel(cores=48)
window_inversion<-foreach(i=1:nrow(bed), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(convergent_cluster,
                       chromosome == as.character(bed[i,1]) &
                         (Taxa1 == as.character(bed[i,4]) | Taxa2 == as.character(bed[i,4])) &
                         ((start >= as.numeric(bed[i,2]) & start <= as.numeric(bed[i,3])) | 
                            (end >= as.numeric(bed[i,2]) & end <= as.numeric(bed[i,3])) | 
                            (start < as.numeric(bed[i,2]) & end > as.numeric(bed[i,3])))),
                data,analysis,variable,direction,comparison,chromosome,range,start,end,size,N_convergent),
         inversion_species=as.character(bed[i,4]),
         inversion_ID=as.character(bed[i,6]),
         inversion_start=as.numeric(bed[i,2]),
         inversion_end=as.numeric(bed[i,3]))}


window_inversion %>%
  mutate(overlap_start=ifelse(start>=inversion_start,start,inversion_start)) %>%
  mutate(overlap_end=ifelse(end>=inversion_end,inversion_end,end)) %>%
  mutate(overlap_size=(overlap_end-overlap_start)+1) %>%
  mutate(inversion_size=(inversion_end-inversion_start)+1) %>%
  select(data,analysis,variable,direction,comparison,chromosome,cluster_start=start,cluster_end=end,cluster_size=size,N_convergent,inversion_species,inversion_ID,inversion_start,
         inversion_end,inversion_size,overlap_start,overlap_end,overlap_size) %>%
  mutate(is_cluster_overlapping="YES") %>%
  mutate(selection_variable=paste0(data,"__",analysis,"__",variable,"__",comparison,"__",chromosome,"__",cluster_start,"__",cluster_end,"__",N_convergent)) -> window_inversion

window_inversion %>%
  distinct(selection_variable) %>%
  pull(selection_variable) -> SELECT

convergent_cluster %>%
  mutate(selection_variable=paste0(data,"__",analysis,"__",variable,"__",comparison,"__",chromosome,"__",start,"__",end,"__",N_convergent)) %>%
  filter(!selection_variable %in% SELECT) %>%
  mutate(is_cluster_overlapping="NO",
         inversion_species="NA",
         inversion_ID="NA",
         inversion_start="NA",
         inversion_end="NA",
         inversion_size="NA",
         overlap_start="NA",
         overlap_end="NA",
         overlap_size=0) %>%
  select(data,analysis,variable,direction,comparison,chromosome,cluster_start=start,cluster_end=end,cluster_size=size,N_convergent,inversion_species,inversion_ID,inversion_start,
         inversion_end,inversion_size,overlap_start,overlap_end,overlap_size,is_cluster_overlapping) %>%
  rbind(.,select(window_inversion,-selection_variable))  -> non_overlap

fwrite(non_overlap,paste0(address,"convergent_cluster_inversion_overlap_LD0.9_1cM_ranges_verbose"), sep = "\t",col.names = T)
