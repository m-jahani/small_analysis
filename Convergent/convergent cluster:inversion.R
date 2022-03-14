library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(gtools)


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

LD <- c(0.8,0.85,0.9,0.95,0.99)
cM <- c(0.1,0.25,0.5,1,2,5,10)

for (z in 1:length(LD)) {
  for (y in 1:length(cM)) {
fread(paste0("/data/users/mjahani/LD_cluster/cluster/convergent_clustering_",LD[z],"_",cM[y],"_cM_summary"),sep = "\t",header = T) %>%
  separate(direction, into = c("Taxa1","Taxa2"),sep = "_",remove = F) %>%
  separate(range,into = c("start","end"),sep = ":",remove = F) %>%
  mutate(start=as.numeric(start),end=as.numeric(end))-> convergent_cluster


registerDoParallel(cores=40)
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

window_inversion %>%
  mutate(overlap_start=ifelse(start>=inversion_start,start,inversion_start)) %>%
  mutate(overlap_end=ifelse(end>=inversion_end,inversion_end,end)) %>%
  mutate(overlap_size=(overlap_end-overlap_start)+1) %>%
  mutate(overlap_percentage=((as.numeric(overlap_size)-1)/(end-start))*100) %>%
  select(-inversion_start,-inversion_end) %>%
  mutate(selection_variable=paste0(data,analysis,variable,direction,chromosome,range,size)) -> window_inversion

window_inversion %>%
  distinct(selection_variable) %>%
  pull(selection_variable) -> loop_list

overlap_percentage_data <- NULL
for (i in 1:length(loop_list)) {
  window_inversion %>%
    filter(selection_variable==loop_list[i]) %>%
    select(-selection_variable) %>% 
    arrange(overlap_start) -> temporary_data

  if (1<nrow(temporary_data)){
    combination_row <- as.data.frame(permutations(n=nrow(temporary_data),r=2,v=seq(1,nrow(temporary_data),1),repeats.allowed=T)) %>%
      mutate(test=(as.numeric(V1)<=as.numeric(V2))) %>% filter(test!=TRUE) %>% select(-test)
    overlap_combination <- NULL
    for (j in 1:nrow(combination_row)) {
      data.frame(row1=combination_row[j,1],
                 row2=combination_row[j,2],
                 size=temporary_data[combination_row[j,1],14],
                 overlap_raw=(min(
                   c(
                     as.numeric(temporary_data[combination_row[j,1],13]),
                     as.numeric(temporary_data[combination_row[j,2],13])
                     )#vector of to overlap end
                   ) - as.numeric(temporary_data[combination_row[j,1],12])+1)
      ) %>%
        mutate(overlap=if_else(as.numeric(overlap_raw)>=0,as.numeric(overlap_raw),0)) %>%
        mutate(non_overlap=as.numeric(size)-as.numeric(overlap)) %>%
        select(-overlap_raw) %>%
        rbind(.,overlap_combination)-> overlap_combination
    }#combination_row
    overlap_combination %>% 
      group_by(row1) %>%
      summarise(a=min(non_overlap)) %>%
      ungroup() %>%
      summarise(overall_overlap_percentage=(100*((sum(a)+as.numeric(temporary_data[1,14]))/as.numeric(temporary_data[1,9]))),
                overall_overlap=(sum(a)+as.numeric(temporary_data[1,14]))) %>%
      mutate(data=temporary_data[1,1],
             analysis=temporary_data[1,2],
             variable=temporary_data[1,3],
             direction=temporary_data[1,4],
             chromosome=temporary_data[1,5],
             range=temporary_data[1,6],
             size=temporary_data[1,9]) %>%
      select(data,analysis,variable,direction,chromosome,range,size,overall_overlap_percentage,overall_overlap) %>%
      rbind(.,overlap_percentage_data)-> overlap_percentage_data
  }else{temporary_data %>%
      select(data,analysis,variable,direction,chromosome,range,size,overall_overlap_percentage=overlap_percentage,overall_overlap=overlap_size) %>%
      rbind(.,overlap_percentage_data)-> overlap_percentage_data
  }
}


overlap_percentage_data %>%
  mutate(selection_variable=paste0(data,analysis,variable,direction,chromosome,range,size)) %>%
  distinct(selection_variable) %>%
  pull(selection_variable) -> list_convergent
  
convergent_cluster %>%
  mutate(selection_variable=paste0(data,analysis,variable,direction,chromosome,range,size)) %>%
  filter(!selection_variable %in% list_convergent) %>%
  mutate(overall_overlap_percentage=0,overall_overlap=0) %>%
  select(data,analysis,variable,direction,chromosome,range,size,overall_overlap_percentage,overall_overlap) %>%
  rbind(.,overlap_percentage_data) -> full_data

full_data %>%
  group_by(data,analysis,variable) %>%
  summarise(N_cluster=n_distinct(range),
            cluster_length=sum(size)) %>%
  ungroup(.) %>% 
  full_join(.,
  full_data %>%
  filter(overall_overlap>0) %>%
  group_by(data,analysis,variable) %>%  
  summarise(overlap_length=sum(overall_overlap),
            N_overlap_cluster=n_distinct(range)) %>%
  ungroup(.)) %>%
  full_join(.,
  convergent_cluster %>% 
  distinct(data,analysis,variable,direction,chromosome,range,.keep_all = T) %>%
  group_by(data,analysis,variable) %>% 
  summarise(N_convergent_window=sum(N_convergent)) %>%
  ungroup(.)) %>%
  full_join(.,
  window_inversion %>%
  group_by(data,analysis,variable) %>%
  summarise(N_overlap_inversion_fragment=n(),
            N_overlap_inversion=n_distinct(inversion)) %>%
  ungroup(.)) %>% 
  mutate(overlap_length=replace_na(overlap_length,0),
         N_overlap_cluster=replace_na(N_overlap_cluster,0),
         N_overlap_inversion_fragment=replace_na(N_overlap_inversion_fragment,0),
         N_overlap_inversion=replace_na(N_overlap_inversion,0)) %>%
  select(data,analysis,variable,N_overlap_inversion_fragment,N_overlap_inversion,
         N_cluster,N_convergent_window,overlap_length,cluster_length,N_overlap_cluster) -> final
  
final %>%
  mutate(Number_overlap_cluster=(N_overlap_cluster/N_cluster),
         length_overlap=(overlap_length/cluster_length)) %>%
  fwrite(paste0("/data/users/mjahani/LD_cluster/inv_con/conversion_inversion_LD",LD[z],"_",cM[y],"cM"), sep = "\t",col.names = F)
  

data.frame(LD=LD[z],
           cM=cM[y],
           overall_cluster_overlap=sum(final$N_overlap_cluster)/sum(final$N_cluster),
           overall_length_overlap=sum(final$overlap_length)/sum(final$cluster_length)) %>%
  fwrite("/data/users/mjahani/LD_cluster/inv_con/conversion_inversion_cutoffs_summary", sep = "\t",col.names = F,append = T)

rm(list_convergent,full_data,final,convergent_cluster,window_inversion,loop_list,temporary_data,overlap_combination,overlap_percentage_data)
  }
}

# library(tidyverse)
# library(data.table)
# LD <- c(0.8,0.85,0.9,0.95,0.99)
# cM <- c(0.1,0.25,0.5,1,2,5,10)
# 
# for (z in 1:length(LD)) {
#   for (y in 1:length(cM)) {
# read.table(paste0("/data/users/mjahani/LD_cluster/inv_con/conversion_inversion_LD",LD[z],"_",cM[y],"cM"), sep = "\t",header = F) %>%
#       select(data=V1,analysis=V2,variable=V3,N_overlap_inversion_fragment=V4,N_overlap_inversion=V5,
#     N_cluster=V6,N_convergent_window=V7,overlap_length=V8,cluster_length=V9,N_overlap_cluster=V10,Number_overlap_cluster=V11,overall_length_overlap=V12) %>%
#       filter(!analysis %in% c("spearman","GWAS_uncorrected"))-> final
#    
# data.frame(LD=LD[z],
#            cM=cM[y],
#            overall_cluster_overlap_proportion=sum(as.numeric(final$N_overlap_cluster))/sum(as.numeric(final$N_cluster)),
#            overall_length_overlap_proportion=sum(as.numeric(final$overlap_length))/sum(as.numeric(final$cluster_length)),
#            overall_cluster_overlap=sum(as.numeric(final$N_overlap_cluster)),
#            N_clusters=sum(as.numeric(final$N_cluster)),
#            overall_length_overlap=sum(as.numeric(final$overlap_length)),
#            total_length=sum(as.numeric(final$cluster_length))
#            ) %>%
#   fwrite("/data/users/mjahani/LD_cluster/inv_con/conversion_inversion_cutoffs_corrected_summary", sep = "\t",col.names = F,append = T)}}




########null distribution
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(gtools)

data.frame(chr=seq(1,17,1),
           end=c(159217232,
                 184765313,
                 182121094,
                 221926752,
                 187413619,
                 157516673,
                 162054088,
                 177936428,
                 200254352,
                 190919061,
                 199477616,
                 179008902,
                 175682153,
                 189916084,
                 186609792,
                 218409126,
                 205797941)) %>%
  mutate(chrom=ifelse(chr<=9,"Ha412HOChr0","Ha412HOChr")) %>%
  mutate(chr=paste0(chrom,chr)) %>%
  select(chr,end) %>%
  fwrite(paste0(address,"Ha412.genome"),sep = "\t",col.names = F)

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(gtools)


"/data/users/mjahani/LD_cluster/inv_con/temprary_run/" -> address  
fread("/data/users/mjahani/inversion.bed",sep = "\t",header = T) -> bed_temp
convert<-data.frame(species=c("Annuus","Argophyllus","petpet","petfal"),
                    spe=c("ann","arg","pet","fal"))
bed_temp %>%
  mutate(sv_new=gsub("ann","",sv_name)) %>%
  mutate(sv_new=gsub("arg","",sv_new)) %>%
  mutate(sv_new=gsub("pet","",sv_new)) %>%
  left_join(.,convert) %>%
  mutate(sv_name_new=paste0(spe,sv_new)) %>%
  select(chr,start,end,species,sv_name_new)-> bed_temp 



library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(gtools)
perm <- seq(9001,10000,1)


"/data/users/mjahani/LD_cluster/inv_con/temprary_run/" -> address
Taxa <- c("Annuus","Argophyllus","petpet","petfal")
for (x in 1:length(perm)) {
bed <- NULL
for (i in 1:length(Taxa)) {
  # bed_temp %>%
  #   filter(species==Taxa[i]) %>%
  #   select(chr,start,end,sv_name_new) %>%
  #   fwrite(paste0(address,Taxa[i],".bed"),sep = "\t",col.names = F)

  system(paste0("bedtools shuffle -i ",address,Taxa[i],".bed -g ",address,"Ha412.genome -noOverlapping > ",address,Taxa[i],"_perm",perm[x],"_new.bed")) 
  
  fread(paste0(address,Taxa[i],"_perm",perm[x],"_new.bed"),sep = "\t",header = F) %>%
    mutate(species=Taxa[i]) %>%
    select(chr=V1,start=V2,end=V3,species,sv_name_new=V4) %>%
    rbind(.,bed)-> bed
  system(paste0("rm ",address,Taxa[i],"_perm",perm[x],"_new.bed"))
}



# LD <- c(0.8,0.85,0.9,0.95,0.99)
# cM <- c(0.1,0.25,0.5,1,2,5,10)

LD <- c(0.9)
cM <- c(1)

for (z in 1:length(LD)) {
  for (y in 1:length(cM)) {
    fread(paste0("/data/users/mjahani/LD_cluster/cluster/convergent_clustering_",LD[z],"_",cM[y],"_cM_summary"),sep = "\t",header = T) %>%
      separate(direction, into = c("Taxa1","Taxa2"),sep = "_",remove = F) %>%
      separate(range,into = c("start","end"),sep = ":",remove = F) %>%
      mutate(start=as.numeric(start),end=as.numeric(end))-> convergent_cluster
    
    
    registerDoParallel(cores=40)
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
    
    window_inversion %>%
      mutate(overlap_start=ifelse(start>=inversion_start,start,inversion_start)) %>%
      mutate(overlap_end=ifelse(end>=inversion_end,inversion_end,end)) %>%
      mutate(overlap_size=(overlap_end-overlap_start)) %>%
      mutate(overlap_percentage=(overlap_size/(end-start))*100) %>%
      select(-inversion_start,-inversion_end) %>%
      mutate(selection_variable=paste0(data,analysis,variable,direction,chromosome,range,size)) -> window_inversion
    
    window_inversion %>%
      distinct(selection_variable) %>%
      pull(selection_variable) -> loop_list
    
    overlap_percentage_data <- NULL
    for (i in 1:length(loop_list)) {
      window_inversion %>%
        filter(selection_variable==loop_list[i]) %>%
        select(-selection_variable) %>% 
        arrange(overlap_start) -> temporary_data
      
      if (1<nrow(temporary_data)){
        combination_row <- as.data.frame(permutations(n=nrow(temporary_data),r=2,v=seq(1,nrow(temporary_data),1),repeats.allowed=T)) %>%
          mutate(test=(as.numeric(V1)<=as.numeric(V2))) %>% filter(test!=TRUE) %>% select(-test)
        overlap_combination <- NULL
        for (j in 1:nrow(combination_row)) {
          data.frame(row1=combination_row[j,1],
                     row2=combination_row[j,2],
                     size=temporary_data[combination_row[j,1],14],
                     overlap_raw=min(
                       c(
                         as.numeric(temporary_data[combination_row[j,1],13]),
                         as.numeric(temporary_data[combination_row[j,2],13])
                       )#vector of to overlap end
                     ) - as.numeric(temporary_data[combination_row[j,1],12])
          ) %>%
            mutate(overlap=if_else(as.numeric(overlap_raw)>=0,as.numeric(overlap_raw),0)) %>%
            mutate(non_overlap=as.numeric(size)-as.numeric(overlap)) %>%
            select(-overlap_raw) %>%
            rbind(.,overlap_combination)-> overlap_combination
        }#combination_row
        overlap_combination %>% 
          group_by(row1) %>%
          summarise(a=min(non_overlap)) %>%
          ungroup() %>%
          summarise(overall_overlap_percentage=(100*((sum(a)+as.numeric(temporary_data[1,14]))/as.numeric(temporary_data[1,9]))),
                    overall_overlap=(sum(a)+as.numeric(temporary_data[1,14]))) %>%
          mutate(data=temporary_data[1,1],
                 analysis=temporary_data[1,2],
                 variable=temporary_data[1,3],
                 direction=temporary_data[1,4],
                 chromosome=temporary_data[1,5],
                 range=temporary_data[1,6],
                 size=temporary_data[1,9]) %>%
          select(data,analysis,variable,direction,chromosome,range,size,overall_overlap_percentage,overall_overlap) %>%
          rbind(.,overlap_percentage_data)-> overlap_percentage_data
      }else{temporary_data %>%
          select(data,analysis,variable,direction,chromosome,range,size,overall_overlap_percentage=overlap_percentage,overall_overlap=overlap_size) %>%
          rbind(.,overlap_percentage_data)-> overlap_percentage_data
      }
    }
    
    
    overlap_percentage_data %>%
      mutate(selection_variable=paste0(data,analysis,variable,direction,chromosome,range,size)) %>%
      distinct(selection_variable) %>%
      pull(selection_variable) -> list_convergent
    
    convergent_cluster %>%
      mutate(selection_variable=paste0(data,analysis,variable,direction,chromosome,range,size)) %>%
      filter(!selection_variable %in% list_convergent) %>%
      mutate(overall_overlap_percentage=0,overall_overlap=0) %>%
      select(data,analysis,variable,direction,chromosome,range,size,overall_overlap_percentage,overall_overlap) %>%
      rbind(.,overlap_percentage_data) -> full_data
    
    full_data %>%
      group_by(data,analysis,variable) %>%
      summarise(N_cluster=n_distinct(range),
                cluster_length=sum(size)) %>%
      ungroup(.) %>% 
      full_join(.,
                full_data %>%
                  filter(overall_overlap>0) %>%
                  group_by(data,analysis,variable) %>%  
                  summarise(overlap_length=sum(overall_overlap),
                            N_overlap_cluster=n_distinct(range)) %>%
                  ungroup(.)) %>%
      full_join(.,
                convergent_cluster %>% 
                  distinct(data,analysis,variable,direction,chromosome,range,.keep_all = T) %>%
                  group_by(data,analysis,variable) %>% 
                  summarise(N_convergent_window=sum(N_convergent)) %>%
                  ungroup(.)) %>%
      full_join(.,
                window_inversion %>%
                  group_by(data,analysis,variable) %>%
                  summarise(N_overlap_inversion_fragment=n(),
                            N_overlap_inversion=n_distinct(inversion)) %>%
                  ungroup(.)) %>% 
      mutate(overlap_length=replace_na(overlap_length,0),
             N_overlap_cluster=replace_na(N_overlap_cluster,0),
             N_overlap_inversion_fragment=replace_na(N_overlap_inversion_fragment,0),
             N_overlap_inversion=replace_na(N_overlap_inversion,0)) %>%
      select(data,analysis,variable,N_overlap_inversion_fragment,N_overlap_inversion,
             N_cluster,N_convergent_window,overlap_length,cluster_length,N_overlap_cluster) -> final
    
    final %>%
      mutate(Number_overlap_cluster=(N_overlap_cluster/N_cluster),
             length_overlap=(overlap_length/cluster_length),
             permutation=perm[x]) %>%
      select(permutation,data,analysis,variable,N_overlap_inversion_fragment,N_overlap_inversion,
             N_cluster,N_convergent_window,overlap_length,cluster_length,N_overlap_cluster,Number_overlap_cluster,length_overlap) %>%
      fwrite(paste0("/data/users/mjahani/LD_cluster/inv_con/conversion_inversion_LD",LD[z],"_",cM[y],"cM_permutation"), sep = "\t",col.names = F,append = T)
    
    
    data.frame(permutation=perm[x],
               LD=LD[z],
               cM=cM[y],
               overall_cluster_overlap=sum(final$N_overlap_cluster)/sum(final$N_cluster),
               overall_length_overlap=sum(final$overlap_length)/sum(final$cluster_length)) %>%
      fwrite("/data/users/mjahani/LD_cluster/inv_con/conversion_inversion_cutoffs_summary_permutation", sep = "\t",col.names = F,append = T)
    
    rm(list_convergent,full_data,final,convergent_cluster,window_inversion,loop_list,temporary_data,overlap_combination,overlap_percentage_data)
  }
}
}



fread(paste0("/data/users/mjahani/LD_cluster/inv_con/conversion_inversion_LD0.9_1cM"),sep = "\t",header = F) %>%
  select(data=V1,analysis=V2,variable=V3,N_overlap_inversion_fragment=V4,N_overlap_inversion=V5,
         N_cluster=V6,N_convergent_window=V7,overlap_length=V8,cluster_length=V9,N_overlap_cluster=V10,Number_overlap_cluster=V11,length_overlap=V12) %>%
  mutate(dav=paste0(data,analysis,variable)) -> data

fread(paste0("/data/users/mjahani/LD_cluster/inv_con/conversion_inversion_LD0.9_1cM_permutation"),sep = "\t",header = F) %>%
  select(permutation=V1,data=V2,analysis=V3,variable=V4,N_overlap_inversion_fragment=V5,N_overlap_inversion=V6,
         N_cluster=V7,N_convergent_window=V8,overlap_length=V9,cluster_length=V10,N_overlap_cluster=V11,Number_overlap_cluster=V12,length_overlap=V13) %>%
  mutate(dav=paste0(data,analysis,variable)) -> data_perm

data %>% 
  distinct(dav) %>%
  pull(dav)-> list_dav
final_result <- NULL

for (i in 1:length(list_dav)) {
  data %>% 
    filter(dav==list_dav[i]) -> a
  
  data_perm %>% 
    filter(dav==list_dav[i]) %>%
    pull(Number_overlap_cluster) -> NOC_list
  
  data_perm %>% 
    filter(dav==list_dav[i]) %>%
    pull(Number_overlap_cluster) -> LO_list
  

  a[1,14] <- mean(NOC_list)
  a[1,15] <- data_perm %>% 
    filter(dav==list_dav[i]) %>%
    filter(Number_overlap_cluster< a[1,11]) %>%
    nrow()
  

  a[1,16] <- mean(LO_list)
  a[1,17] <- data_perm %>% 
    filter(dav==list_dav[i]) %>%
    filter(Number_overlap_cluster < a[1,12]) %>%
    nrow()
  
  a %>% rbind(.,final_result) -> final_result
  rm(a,NOC_list,LO_list)
}

final_result %>% 
  select(data,analysis,variable,Number_overlap_cluster,Number_overlap_cluster_random_mean=V14,
         Number_overlap_cluster_P=V15,length_overlap,length_overlap_random_mean=V16,length_overlap_P=V17) %>%
  mutate(Number_overlap_cluster_P=Number_overlap_cluster_P/10000) %>%
  mutate(length_overlap_P=length_overlap_P/10000) %>%
  mutate(Number_overlap_cluster_effectsize=Number_overlap_cluster/Number_overlap_cluster_random_mean) %>%
  mutate(length_overlap_effectsize=length_overlap/length_overlap_random_mean) %>%
  select(data,analysis,variable,
         Number_overlap_cluster,Number_overlap_cluster_P,Number_overlap_cluster_effectsize,Number_overlap_cluster_random_mean,
         length_overlap,length_overlap_P,length_overlap_effectsize,length_overlap_random_mean) %>%
  fwrite("/data/users/mjahani/LD_cluster/inv_con/conversion_inversion_LD0.9_1cM_final", sep = "\t",col.names = T)
  
         
  

