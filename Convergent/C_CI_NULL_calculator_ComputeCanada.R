library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(IRanges)

args = commandArgs(trailingOnly = TRUE)

address <- "~/C_CI_input/"

fread(paste0(address,"convergent_clustering_0.9_1_cM_summary"),sep = "\t",header = T) %>%
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

perm <- seq(1,25,1)
Taxa <- c("Annuus","Argophyllus","petpet","petfal")

for (x in 1:length(perm)) {
  bed <- NULL
  for (i in 1:length(Taxa)) {
    system(paste0("bedtools shuffle -i ",address,Taxa[i],".bed -g ",address,"Ha412.genome -noOverlapping > ",address,Taxa[i],"_perm",perm[x],args[1],"_new.bed")) 
    
    fread(paste0(address,Taxa[i],"_perm",perm[x],args[1],"_new.bed"),sep = "\t",header = F) %>%
      mutate(species=Taxa[i]) %>%
      select(chr=V1,start=V2,end=V3,species,sv_name_new=V4) %>%
      rbind(.,bed)-> bed
    system(paste0("rm ",address,Taxa[i],"_perm",perm[x],args[1],"_new.bed"))
  }
  
  registerDoParallel(cores=48)
  window_inversion<-foreach(i=1:nrow(bed), .combine='rbind', .errorhandling='stop') %dopar% { 
    mutate(select(filter(convergent_cluster,
                         chromosome == as.character(bed[i,1]) &
                           (Taxa1 == as.character(bed[i,4]) | Taxa2 == as.character(bed[i,4])) &
                           ((start >= as.numeric(bed[i,2]) & start <= as.numeric(bed[i,3])) | 
                              (end >= as.numeric(bed[i,2]) & end <= as.numeric(bed[i,3])) | 
                              (start < as.numeric(bed[i,2]) & end > as.numeric(bed[i,3])))),
                  data,analysis,variable,direction,chromosome,range,start,end,size,N_convergent),
           inversion=as.character(bed[i,5]),
           inversion_start=as.numeric(bed[i,2]),
           inversion_end=as.numeric(bed[i,3]))}
  rm(bed)
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
  
  rbind(zero_overlap,partial_non_overlap) -> C_regions  #nonoverlap_regions
 
  rm(zero_overlap,partial_non_overlap,data,SELEC)
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
  
  rm(overlap_regions)
  non_overlaps %>%
    mutate(start_C=as.numeric(start_C),
           end_C=as.numeric(end_C),
           start_CI=as.numeric(start_CI),
           end_CI=as.numeric(end_CI)) %>%
    select(data,analysis,variable,chromosome,original_cluster,direction_C,start_C,end_C,start_CI,end_CI,direction_CI,inversion) -> non_overlaps
  
  
  non_overlaps %>%
    mutate(overlap_start=ifelse(start_C>=start_CI,start_C,start_CI)) %>%
    mutate(overlap_end=ifelse(end_C>=end_CI,end_CI,end_C)) %>%
    select(-start_CI,-end_CI) %>% 
    mutate(selection=paste0(data,"__",analysis,"__",variable,"__",chromosome,"__",original_cluster,"__",direction_C,"__",start_C,"__",end_C)) %>%
    select(data,analysis,variable,chromosome,original_cluster,direction_C,start_C,end_C,overlap_start,overlap_end,direction_CI,inversion,selection)-> nonoverlap_regions_overlap
  
  rm(non_overlaps) 
  nonoverlap_regions_overlap %>% distinct(selection) %>% pull(selection) -> SELEC_NON
  overlap_length <- NULL
  for (i in 1:length(SELEC_NON)) {
    nonoverlap_regions_overlap %>%
      filter(selection==SELEC_NON[i]) -> tmp2
    
    reduce(IRanges(pull(tmp2,overlap_start),pull(tmp2,overlap_end))) %>%
      data.frame() %>%
      summarise(overlap_length=sum(width)) %>%
      mutate(selection=SELEC_NON[i]) %>%
      separate(selection,into = c("data","analysis","variable","chromosome","original_cluster","direction_C","start_C","end_C"),sep = "__",remove = F) %>%
      mutate(start_C=as.numeric(start_C),
             end_C=as.numeric(end_C)) %>%
      mutate(length_C=(end_C-start_C)+1) %>%
      mutate(overlap_porportion=(overlap_length/length_C)) %>% #what porportion of a convergent region(not overlab with inversion) is overlap with inversion in other contrasts
      select(data,analysis,variable,chromosome,original_cluster,start_C,end_C,direction_C,length_C,overlap_length,overlap_porportion,selection) %>%
      rbind(.,overlap_length) -> overlap_length
    rm(tmp2)}
  
  #only for NULL make a facke data frame which contains  all data,analysis,variable,directoins and fake chromosomem,cluster,no_overlap_start,no_overlap_end
  
  overlap_length %>% 
    distinct(selection) %>%
    pull(selection) -> selection_variables
  #add zero overlap between convergent reions(C) and convergent_inversion regions(CI) to the result
  C_regions %>%
    mutate(selection=paste0(data,"__",analysis,"__",variable,"__",chromosome,"__",cluster,"__",direction,"__",no_overlap_start,"__",no_overlap_end)) %>%
    filter(!selection %in% selection_variables) %>%  
    mutate(length_C=(as.numeric(no_overlap_end)-as.numeric(no_overlap_start))+1) %>%
    mutate(overlap_length=0,overlap_porportion=0) %>%
    select(data,analysis,variable,chromosome,original_cluster=cluster,
           start_C=no_overlap_start,end_C=no_overlap_end,direction_C=direction,
           length_C,overlap_length,overlap_porportion,selection) %>%
    rbind(.,overlap_length) -> all_nonoverlap_regions_overlaps
  rm(selection_variables)
  
  
  #to fill zero overlap data
  convergent_cluster %>% 
    mutate(original_cluster=paste0(chromosome,":",range)) %>%
    group_by(data,analysis,variable,direction) %>%
    summarise(N_original_cluster=n_distinct(original_cluster)) %>%
    ungroup() %>%
    mutate(length_C=0,
           N_original_cluster_overlab_CI=0,
           length_C_overlap_CI=0,
           porportion_number_original_cluster_overlap_CI=0,
           porportion_length_C_overlap_CI=0) %>%
    dplyr::rename(direction_C=direction) %>%
    mutate(selection=paste0(data,"__",analysis,"__",variable,"__",direction_C)) -> FILL_direction
  
  all_nonoverlap_regions_overlaps %>%
    group_by(data,analysis,variable,direction_C) %>%
    summarise(N_original_cluster=n_distinct(original_cluster),#Number of original convergent clusters
              length_C=sum(as.numeric(length_C))) %>% #sum of all the convergent(not overlap with inversion) regions
    ungroup(.) %>% 
    full_join(.,
              all_nonoverlap_regions_overlaps %>%
                filter(as.numeric(overlap_length) > 0) %>%
                group_by(data,analysis,variable,direction_C) %>%
                summarise(N_original_cluster_overlab_CI=n_distinct(original_cluster),#Number of original convergent clusters which (either fully or partially) overlap with conversion_inversion regions in the other contrasts 
                          length_C_overlap_CI=sum(as.numeric(overlap_length))#length of the portion of all the convergent(not overlap with inversion) regions which overlap with conversion inversion in other contrasts
                )) %>%
    mutate(N_original_cluster_overlab_CI=replace_na(N_original_cluster_overlab_CI,0)) %>%
    mutate(length_C_overlap_CI=replace_na(length_C_overlap_CI,0)) %>%
    mutate(porportion_number_original_cluster_overlap_CI=(N_original_cluster_overlab_CI/N_original_cluster),#what porportion of original clusters overlap with conversion_inversion regions in the other contrasts (either fully or partially)
           porportion_length_C_overlap_CI=(length_C_overlap_CI/length_C)#what porportion of convergent(not overlap with inversion) regions overlap with conversion inversion in other contrasts
    ) %>%
    mutate(selection=paste0(data,"__",analysis,"__",variable,"__",direction_C)) -> C_CI_overlap_direction
  
  C_CI_overlap_direction %>%
    distinct(selection) %>%
    pull(selection)-> SELEC_direction
  
  FILL_direction %>%
    filter(!selection %in% SELEC_direction) %>%
    rbind(.,
          C_CI_overlap_direction) %>%
    select(-selection) %>%
    fwrite(paste0(address,"C_CI_overlap_LD0.9_1cM_direction_NULL"), sep = "\t",col.names = F, append = T)
  rm(FILL_direction,C_CI_overlap_direction,SELEC_direction)
  
  convergent_cluster %>% 
    mutate(original_cluster=paste0(chromosome,":",range)) %>%
    group_by(data,analysis,variable) %>%
    summarise(N_original_cluster=n_distinct(original_cluster)) %>%
    ungroup() %>%
    mutate(length_C=0,
           N_original_cluster_overlab_CI=0,
           length_C_overlap_CI=0,
           porportion_number_original_cluster_overlap_CI=0,
           porportion_length_C_overlap_CI=0) %>%
    mutate(selection=paste0(data,"__",analysis,"__",variable)) -> FILL
  
  all_nonoverlap_regions_overlaps %>%
    group_by(data,analysis,variable) %>%
    summarise(N_original_cluster=n_distinct(original_cluster),#Number of original convergent clusters
              length_C=sum(as.numeric(length_C))) %>% #sum of all the convergent(not overlap with inversion) regions
    ungroup(.) %>% 
    full_join(.,
              all_nonoverlap_regions_overlaps %>%
                filter(as.numeric(overlap_length) > 0) %>%
                group_by(data,analysis,variable) %>%
                summarise(N_original_cluster_overlab_CI=n_distinct(original_cluster),#Number of original convergent clusters which (either fully or partially) overlap with conversion_inversion regions in the other contrasts 
                          length_C_overlap_CI=sum(as.numeric(overlap_length))#length of the portion of all the convergent(not overlap with inversion) regions which overlap with conversion inversion in other contrasts
                )) %>%
    mutate(N_original_cluster_overlab_CI=replace_na(N_original_cluster_overlab_CI,0)) %>%
    mutate(length_C_overlap_CI=replace_na(length_C_overlap_CI,0)) %>%
    mutate(porportion_number_original_cluster_overlap_CI=(N_original_cluster_overlab_CI/N_original_cluster),#what porportion of original clusters overlap with conversion_inversion regions in the other contrasts (either fully or partially)
           porportion_length_C_overlap_CI=(length_C_overlap_CI/length_C)#what porportion of convergent(not overlap with inversion) regions overlap with conversion inversion in other contrasts
    )  %>%
    mutate(selection=paste0(data,"__",analysis,"__",variable)) -> C_CI_overlap
  
  C_CI_overlap %>%
    distinct(selection) %>%
    pull(selection) -> SELEC
  
  FILL %>%
    filter(!selection %in% SELEC) %>%
    rbind(.,
          C_CI_overlap) %>%
    select(-selection) %>%
    fwrite(paste0(address,"C_CI_overlap_LD0.9_1cM_NULL"), sep = "\t",col.names = F, append = T)
  
  
  rm(SELEC,C_CI_overlap,FILL,nonoverlap_regions,nonoverlap_regions_overlap,SELEC_NON,overlap_length,all_nonoverlap_regions_overlaps)}

