#directions have been merged (ann-arg and arg_ann were merged to ann_arg)
#note: if ann_arg and arg_ann have the same convergent island we consider it only once
#for cluster number convergent , it has been count based on original convergent cluster numbers not C_resions
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(IRanges)
address <- "/data/users/mjahani/C_CI_overlap/"
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
  select(-direction) %>%
  dplyr::rename(direction=comparison) %>%
  distinct(data,analysis,variable,chromosome,range,size,N_convergent,direction) %>%
  separate(direction, into = c("Taxa1","Taxa2"),sep = "_",remove = F) %>%
  separate(range,into = c("start","end"),sep = ":",remove = F) %>%
  mutate(start=as.numeric(start),end=as.numeric(end))-> convergent_cluster


#overlap between inversion and conversion regions(including partial overlap)
registerDoParallel(cores=48)
window_inversion<-foreach(i=1:nrow(bed), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(convergent_cluster,
                       chromosome == as.character(bed[i,1]) &
                         (Taxa1 == as.character(bed[i,4]) | Taxa2 == as.character(bed[i,4])) &
                         ((start >= as.numeric(bed[i,2]) & start <= as.numeric(bed[i,3])) | 
                            (end >= as.numeric(bed[i,2]) & end <= as.numeric(bed[i,3])) | 
                            (start < as.numeric(bed[i,2]) & end > as.numeric(bed[i,3])))),
                data,analysis,variable,direction,chromosome,range,start,end,size,N_convergent),
         inversion=as.character(bed[i,6]),
         inversion_start=as.numeric(bed[i,2]),
         inversion_end=as.numeric(bed[i,3]))}

#regions of convergent clusters that ovelaps with inversion in the same species contrast (one of two species) 
window_inversion %>%
  mutate(cluster=paste0(chromosome,":",range)) %>%
  mutate(overlap_start=ifelse(start>=inversion_start,start,inversion_start)) %>%
  mutate(overlap_end=ifelse(end>=inversion_end,inversion_end,end)) %>%
  select(-inversion_start,-inversion_end) %>%  
  select(data,analysis,variable,direction,cluster,chromosome,overlap_start,overlap_end,inversion)-> CI_regions  #overlap_regions
#convergent clusters that completely are not ovelaps with inversion in the same species contrast (one of two species) 
window_inversion %>% 
  mutate(selection=paste0(data,analysis,variable,direction,chromosome,range,start,end,size,N_convergent)) %>% pull(selection)-> OVER_WIN
convergent_cluster %>%
  mutate(selection=paste0(data,analysis,variable,direction,chromosome,range,start,end,size,N_convergent)) %>%
  filter(!selection %in% OVER_WIN) %>%
  mutate(cluster=paste0(chromosome,":",range)) %>%
  select(data,analysis,variable,direction,cluster,chromosome,no_overlap_start=start,no_overlap_end=end) -> zero_overlap
rm(OVER_WIN)

window_inversion %>%
  mutate(cluster=paste0(chromosome,":",range)) %>%
  mutate(overlap_start=ifelse(start>=inversion_start,start,inversion_start)) %>%
  mutate(overlap_end=ifelse(end>=inversion_end,inversion_end,end)) %>%
  mutate(selection=paste0(data,"__",analysis,"__",variable,"__",direction,"__",cluster,"__",chromosome)) %>%
  select(data,analysis,variable,direction,cluster,chromosome,start,end,overlap_start,overlap_end,selection)-> data

data %>% distinct(selection) %>% pull(selection) -> SELEC
partial_non_overlap <- NULL
for (i in 1:length(SELEC)) {
  data %>%
    filter(selection==SELEC[i]) -> tmp1
  #regions in convergent clusters that are not ovelaps with inversion in the same species contrast (one of two species)
  gaps(IRanges(pull(tmp1,overlap_start),pull(tmp1,overlap_end)),start=pull(distinct(tmp1,start),start),end = pull(distinct(tmp1,end),end)) %>%
    data.frame() %>%
    rbind(data.frame(start=c("tmp"),end=c("tmp"),width=c("tmp"))) %>%
    mutate(selection=SELEC[i]) %>%
    separate(selection,into = c("data","analysis","variable","direction","cluster","chromosome"),sep = "__",remove = T) %>%
    select(data,analysis,variable,direction,cluster,chromosome,no_overlap_start=start,no_overlap_end=end) %>%
    filter(no_overlap_start != "tmp" ) %>%
    rbind(.,partial_non_overlap) -> partial_non_overlap
  rm(tmp1)}

rbind(zero_overlap,partial_non_overlap)  %>% 
  distinct(data,analysis,variable,direction,cluster,chromosome,no_overlap_start,no_overlap_end) -> C_regions  #nonoverlap_regions

#overlaps between non-overlab regions(convergent only) and overlap regions(conversion_inversion) 
registerDoParallel(cores=48)
non_overlaps<-foreach(i=1:nrow(C_regions), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(CI_regions,
                       analysis == as.character(C_regions[i,2]) &
                         variable == as.character(C_regions[i,3]) &
                         direction != as.character(C_regions[i,4]) &
                         chromosome == as.character(C_regions[i,6]) &
                         ((overlap_start >= as.numeric(C_regions[i,7]) & overlap_start <= as.numeric(C_regions[i,8])) | 
                            (overlap_end >= as.numeric(C_regions[i,7]) & overlap_end <= as.numeric(C_regions[i,8])) | 
                            (overlap_start < as.numeric(C_regions[i,7]) & overlap_end > as.numeric(C_regions[i,8])))),
                data,analysis,variable,chromosome,direction_CI=direction,start_CI=overlap_start,end_CI=overlap_end,inversion),
         direction_C=as.character(C_regions[i,4]),
         start_C=as.character(C_regions[i,7]),
         end_C=as.character(C_regions[i,8]),
         original_cluster=as.character(C_regions[i,5]))}


non_overlaps %>%
  mutate(start_C=as.numeric(start_C),
         end_C=as.numeric(end_C),
         start_CI=as.numeric(start_CI),
         end_CI=as.numeric(end_CI)) %>%
  select(data,analysis,variable,chromosome,original_cluster,direction_C,start_C,end_C,start_CI,end_CI,direction_CI,inversion) -> non_overlaps


non_overlaps %>%
  mutate(overlap_start=ifelse(start_C>=start_CI,start_C,start_CI)) %>%
  mutate(overlap_end=ifelse(end_C>=end_CI,end_CI,end_C)) %>%
  mutate(selection=paste0(data,"__",analysis,"__",variable,"__",chromosome,"__",original_cluster,"__",direction_C,"__",start_C,"__",end_C)) %>%
  select(data,analysis,variable,chromosome,original_cluster,direction_C,start_C,end_C,direction_CI,start_CI,end_CI,overlap_start,overlap_end,inversion,selection)-> nonoverlap_regions_overlap

nonoverlap_regions_overlap %>%
  distinct(selection) %>%
  pull(selection) -> SELECT

C_regions %>% 
  mutate(selection=paste0(data,"__",analysis,"__",variable,"__",chromosome,"__",cluster,"__",direction,"__",no_overlap_start,"__",no_overlap_end)) %>%
  filter(!selection %in% SELECT) %>%
  mutate(direction_CI="NA",
         start_CI="NA",
         end_CI="NA",
         overlap_start="NA",
         overlap_end="NA",
         inversion="NA") %>%
  select(data,
         analysis,
         variable,
         chromosome,
         original_cluster=cluster,
         direction_C=direction,
         start_C=no_overlap_start,
         end_C=no_overlap_end,
         direction_CI,
         start_CI,
         end_CI,
         overlap_start,
         overlap_end,
         inversion,
         selection) %>%
  rbind(.,nonoverlap_regions_overlap) %>%
  select(data,
         analysis,
         variable,
         chromosome,
         original_cluster,
         comparison_C=direction_C,
         start_C,
         end_C,
         comparison_CI=direction_CI,
         start_CI,
         end_CI,
         overlap_start,
         overlap_end,
         inversion) -> C_CI_LD0.9_1cM_ranges_verbose
  fwrite(C_CI_LD0.9_1cM_ranges_verbose,"/data/users/mjahani/C_CI_overlap/C_CI_LD0.9_1cM_ranges_verbose", sep = "\t",col.names = T)




