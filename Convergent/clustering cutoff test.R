library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)

fread("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/10k_random_all_limitless",header = F) %>%
  select(Taxa=V1,distance=V2,count=V3)-> tenk_random_all
cM <- c(0.1,0.25,0.5,1,2,5,10)
LD <- c(0.8,0.85,0.90,0.95,0.99)
TAXA <- c("Annuus","Argophyllus","petpet","petfal")
for (j in 1:length(LD)) {
for (i in 1:length(TAXA) ) {
  fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/",TAXA[i],"_distance_meanLD_limitless"),header = F) %>%
    select(distance=V1,mean_LD=V2) %>%
    mutate(Taxa=TAXA[i])  %>%
    left_join(.,tenk_random_all) %>%
    filter(count>=10000) %>%
    group_by(Taxa,distance) %>%
    sample_n(10000) %>%
    ungroup() %>%
    group_by(Taxa,distance) %>%  
    summarise(threshold=quantile(mean_LD, LD[j])) %>%
    ungroup() %>%
    fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/Taxa_distance_threshold_limitless_",LD[j]),
           sep = "\t",col.names= F,append = T)
  
  fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/",TAXA[i],"_distance_meanLD_limitless"),header = F) %>%
    select(distance=V1,mean_LD=V2) %>%
    mutate(Taxa=TAXA[i])  %>%
    left_join(.,tenk_random_all) %>%
    filter(count<10000) %>%
    group_by(Taxa,distance) %>%  
    summarise(threshold=quantile(mean_LD, LD[j])) %>%
    ungroup() %>%
    fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/Taxa_distance_threshold_limitless_",LD[j]),
           sep = "\t",col.names= F,append = T)}}


fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/temp/convergent_result"),header = F) %>%
  select(Taxa=V1,window_A=V2,window_B=V3,mean_LD=V4) %>%
  separate(window_A,into = c("chrom_A","start_end_A"), sep=":", remove = F) %>%
  separate(start_end_A,into = c("start_A","end_A"), sep="-", remove = T) %>%
  separate(window_B,into = c("chrom_B","start_end_B"), sep=":", remove = F) %>%
  separate(start_end_B,into = c("start_B","end_B"), sep="-", remove = T) %>%
  mutate(distance=abs(as.numeric(start_A)-as.numeric(start_B))) %>%
  select(Taxa,window_A,window_B,mean_LD,distance)-> conv_res

registerDoParallel(cores=20) 
  foreach(j = 1:length(LD),
          .combine = 'rbind',
          .errorhandling = 'stop') %dopar% {
fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/Taxa_distance_threshold_limitless_",LD[j]),header = F) %>%
  select(Taxa=V1,distance=V2,threshold=V3) %>%
  left_join(conv_res,.) %>%
  mutate(significant_LD=(mean_LD>=threshold)) %>%
  fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/convergent_result_sig_limitless_",LD[j]),sep = "\t")}
  
  rm(conv_res,tenk_random_all)
            
            
fread("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/temp/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") -> map
            map %>%
              group_by(chr) %>%
              arrange(pos,.by_group = TRUE) %>%
              mutate(delta_cM=cM-lag(cM, default = first(cM)),
                     delta_bp=pos-lag(pos, default = first(pos)),
                     start=lag(pos)) %>% 
              ungroup() %>%
              mutate(recomb_rate=delta_cM/delta_bp) %>%
              select(chr,start,end=pos,recomb_rate) %>%
              filter(!is.na(start)) %>%
              left_join(.,rename(map,start=pos)) -> recombination
            rm(map)
  
for (j in 1:length(LD)) {
fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/convergent_result_sig_limitless_",LD[j]),sep = "\t")-> convergent_result                      
distinct(rbind(rename(distinct(convergent_result,window_A),window=window_A),rename(distinct(convergent_result,window_B),window=window_B)),window) %>%
              separate(window,into = c("chr","start_end"), sep=":", remove = F) %>%
              separate(start_end,into = c("start","end"), sep="-", remove = T) %>% 
              mutate(start=as.numeric(start),end=as.numeric(end)) %>%
              mutate(pos=(((end-start)/2)+start)) %>%
              select(window,chr,pos) -> conv_pos #window  chr pos
registerDoParallel(cores=40)
window_recombination<-foreach(i=1:nrow(recombination), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(conv_pos,
                       chr == as.character(recombination[i,1])  &
                         pos >= as.numeric(recombination[i,2])&
                         pos <= as.numeric(recombination[i,3])),
                window,pos),
         start=as.numeric(recombination[i,2]),
         recomb_rate=as.numeric(recombination[i,4]),
         cM=as.numeric(recombination[i,5]))}
rm(conv_pos)
window_recombination %>% 
  mutate(G_position=(((pos-start)*recomb_rate)+cM)) %>%
  select(window,G_position) -> window_recombination # window G_position

registerDoParallel(cores=20) 
foreach(k = 1:length(cM),
        .combine = 'rbind',
        .errorhandling = 'stop') %dopar% {
convergent_result %>%
  left_join(.,
            rename(window_recombination,
                   window_A = window,
                   G_position_A = G_position)) %>%
  left_join(.,
            rename(window_recombination,
                   window_B = window,
                   G_position_B = G_position)) %>%
  mutate(Genetic_distance = abs(G_position_A-G_position_B)) %>%
  mutate(GD_significant = ifelse(Genetic_distance <=  cM[k] ,
                                 "TRUE",
                                 "FALSE")) %>%
  select(-G_position_A,-G_position_B) %>%
  filter(distance != 0) %>%
  mutate(LD_GD_sig = ifelse(significant_LD == "TRUE" & GD_significant == "TRUE",
                            "TRUE",
                            "FALSE")) %>%
  fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/temporary_convergent_result_sig_limitless_",LD[j],"_",cM[k]),sep = "\t")}}

  rm(convergent_result)
            
 for (j in 1:length(LD)) {
   for (k in 1:length(cM)) {
     fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/temporary_convergent_result_sig_limitless_",LD[j],"_",cM[k]),sep = "\t",header = T) -> convergent_result                      
   rbind(convergent_result,
         select(convergent_result,
                Taxa,window_A=window_B,window_B=window_A,mean_LD,distance,threshold,significant_LD,Genetic_distance,GD_significant,LD_GD_sig)) %>%
     separate(window_A,into = c("chrom","start_end_A"), sep=":", remove = F) %>%
     select(Taxa,chrom,window_A,window_B,LD_GD_sig) %>% 
     fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/convergent_result_",LD[j],"_",cM[k],"_cM"),sep = "\t")
   rm(convergent_result)  
   
 }  }         
            

  library(tidyverse)
  library(data.table)
  library(foreach)
  library(doParallel)  

  cM <- c(0.1)#1
  LD <- c(0.8,0.85)
  
  cM <- c(0.25)#2
  LD <- c(0.8,0.85)
  
  cM <- c(0.5) #3
  LD <- c(0.8,0.85)
  
  cM <- c(1)#4
  LD <- c(0.8,0.85)
  
  cM <- c(2)#5
  LD <- c(0.8,0.85)
  
  cM <- c(5)#6
  LD <- c(0.8,0.85)
  
  cM <- c(10)#t7
  LD <- c(0.8,0.85)
  
  cM <- c(0.1)#8
  LD <- c(0.90,0.95,0.99)
  
  cM <- c(0.25)#9
  LD <- c(0.90,0.95,0.99)
  
  cM <- c(0.5) #10
  LD <- c(0.90,0.95,0.99)
  
  cM <- c(1)#11
  LD <- c(0.90,0.95,0.99)
  
  cM <- c(2)#12
  LD <- c(0.90,0.95,0.99)
  
  cM <- c(5)#13
  LD <- c(0.90,0.95,0.99)
  
  cM <- c(10)#t14
  LD <- c(0.90,0.95,0.99)
  
  LD <- c(0.8,0.85,0.90,0.95,0.99)
  cM <- c(0.1,0.25,0.5,1,2,5,10)
  
  # data_thresholds <- NULL
  # for (i in 1:length(LD)) {
  #   for (j in 1:length(cM)) {
  #     fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/convergent_result_",LD[i],"_",cM[j],"_cM"),
  #           sep = "\t",header = T) %>%
  #       mutate(LD=LD[i],cM=cM[j]) %>%
  #       select(LD,cM,Taxa,chrom,window_A,window_B,LD_GD_sig) %>%
  #       rbind(.,data_thresholds)-> data_thresholds}}
  library(tidyverse)
  library(data.table)
  library(foreach)
  library(doParallel)  
  
  
  fread("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/temp/convergent_2binning",sep = "\t",header=T) %>%
    select(type,analysis,Taxa1,Taxa2,variable,window5,Emp_p) -> convergent_window
  
  convergent_window %>% distinct(type) %>% pull(type) -> TYPE
  for (i in 1:length(TYPE)) {
    convergent_window %>%
      filter(type==TYPE[i]) %>%
      select(-type) -> middle_data0
    middle_data0 %>% distinct(analysis) %>% pull(analysis) -> ANALAY
    for (j in 1:length(ANALAY)) {
      middle_data0 %>%
        filter(analysis==ANALAY[j]) %>%
        select(-analysis) -> middle_data
      middle_data %>% distinct(variable) %>% pull(variable) -> VARI
      for (k in 1:length(VARI)) {
        middle_data %>%
          filter(variable==VARI[k]) %>%
          select(-variable)-> middle_data1 
        middle_data1 %>% distinct(Taxa2) %>% pull(Taxa2) -> TAXA2
        for (l in 1:length(TAXA2)) {
          middle_data1 %>% 
            filter(Taxa2==TAXA2[l]) %>%
            select(-Taxa2)-> middle_data2 
          middle_data2 %>% distinct(Taxa1) %>% pull(Taxa1) -> TAXA1
          for (m in 1:length(TAXA1)) {
            middle_data2 %>%
              filter(Taxa1==TAXA1[m]) %>%
              select(-Taxa1) %>%
              separate(window5,into = c("chrom","start_end"), sep=":", remove = F) %>%
              separate(start_end,into = c("start","end"), sep="-", remove = T) %>%
              mutate(start=as.numeric(start)) %>%
              select(-end) ->  middle_data3
            middle_data3 %>% distinct(chrom) %>% pull(chrom) -> CHR
            for (n in 1:length(CHR)) {
              middle_data3 %>%
                filter(chrom==CHR[n]) %>%
                arrange(start) %>%
                select(-chrom,-start) -> middle_data4
              
              for (q in 1:length(LD)) {
                for (r in 1:length(cM)) {
              fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/convergent_result_",LD[q],"_",cM[r],"_cM"),
                        sep = "\t",header = T) %>%     
                filter(Taxa==TAXA1[m]) %>%
                filter(chrom==CHR[n]) %>%
                select(window_A,window_B,LD_GD_sig) -> con_sig
              cluster <- NULL
              
              if (1<nrow(middle_data4)){
                for (o in 1:(nrow(middle_data4)-1)) {
                  rbind(data.frame(Taxa=TAXA1[m],window_A=middle_data4[o,1],window_B=middle_data4[o+1,1],cycle=o),cluster)-> cluster
                }#rows
              }else{rbind(data.frame(Taxa=TAXA1[m],window_A=middle_data4[1,1],window_B="-",cycle=1),cluster)-> cluster}
              left_join(cluster,con_sig) %>%
                mutate(LD_GD_sig=replace_na(LD_GD_sig,"FALSE")) %>% arrange(cycle) %>% mutate(clustering=0) -> cluster
              
              for (p in 1:nrow(cluster)) {
                cluster[p,6] <- ifelse(cluster[p,4]==1,1,
                                       cluster[p,6] <- ifelse(cluster[p,6]== 0 & cluster[p-1,5]=="TRUE" & cluster[p,5]=="TRUE" ,cluster[p-1,6],as.numeric(cluster[p-1,6])+1))
              }#cluster rows
              cluster %>% filter(LD_GD_sig=="TRUE") %>% gather(windows,window_name,2:3) %>% select(window_name,clustering) %>% distinct(window_name,.keep_all = T) -> true_cluster
              middle_data4 %>%
                filter(!window5 %in% pull(distinct(true_cluster,window_name),window_name)) %>%
                mutate(clustering=window5) %>%
                select(window_name=window5,clustering) %>%
                rbind(.,true_cluster) %>%
                group_by(clustering) %>%
                mutate(cluster=group_indices()) %>%
                ungroup() %>%
                mutate(type=TYPE[i],analysis=ANALAY[j],variable=VARI[k],Taxa1=TAXA1[m],Taxa2=TAXA2[l],chrom=CHR[n]) %>%
                select(type,analysis,variable,Taxa1,Taxa2,chrom,window_name,cluster) %>%
                fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/convergent_clustering_",LD[q],"_",cM[r],"_cM"),sep = "\t",col.names= F,append = T)
                }#LD cuttoff
              }#cM
            }#chrom
          }#taxa1
        }#taxa2
      }#variable
    }#analysis
  }#type
  
  for (q in 1:length(LD)) {
    for (r in 1:length(cM)) {  
  fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/convergent_clustering_",LD[q],"_",cM[r],"_cM"),header=F) %>%
    select(type=V1,analysis=V2,variable=V3,Taxa1=V4,Taxa2=V5,chrom=V6,window_name=V7,cluster=V8) %>%
    separate(window_name,into = c("chr","start_end"), sep=":", remove = F) %>%
    select(-chr) %>%
    separate(start_end,into = c("start","end"), sep="-", remove = T) %>% 
    mutate(start=as.numeric(start),end=as.numeric(end)) %>%
    group_by(type,analysis,variable,Taxa1,Taxa2,chrom,cluster) %>%
    summarise(N_convergent=n(),cluster_start=min(start),cluster_end=max(end)) %>%
    ungroup() %>%
    mutate(range=paste0(as.character(cluster39_start),":",as.character(cluster_end))) %>%
    mutate(size=(as.numeric(cluster_end)-as.numeric(cluster_start))+1) %>%
    mutate(direction=paste0(Taxa1,"_",Taxa2)) %>%
    select(data=type,analysis,variable,direction,chromosome=chrom,range,size,N_convergent) %>%
    fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/convergent_clustering_",LD[q],"_",cM[r],"_cM_summary"),sep = "\t",col.names= T,append = T)}}
  
            
  
  x= 2265592-868987
  t= 5605
  (((1396640*t)/(x-422))/60)-(t/60)
  (x/1396640)*100
  

#   
#   fread(paste0("/data/users/mjahani/LD_cluster/cluster/convergent_clustering_0.95_5cM_summary"),header=T) -> a
#   #vector for chromosomes with more than one convergent window
#   a %>%
#     group_by(data,analysis,variable,direction,chromosome) %>% 
#     summarise(count=sum(N_convergent)) %>%
#     ungroup() %>%
#     filter(count>1) %>%
#     mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
#     distinct(davdc) %>%
#     pull(davdc) -> vect
#   
#   #filter for the vector
#   a %>%
#     mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
#     filter(davdc %in% vect) %>%
#     select(-davdc) -> b
#   
# test <- NULL  
#   b %>% summarise(total_window=sum(N_convergent)) -> test[1,1] #maximum posible number of clusters
#   b %>% distinct(data,analysis,variable,direction,chromosome) %>% nrow -> test[1,2] #minimum posible number of clusters
#   b %>% summarise(N_cluster=n()) -> test[1,3] #number of clusters
#   b %>% summarise(max_size=max(size)) -> test[1,4] #maximum cluster size
#   b %>% group_by(data,analysis,variable,direction,chromosome) %>% tally() %>% ungroup() %>% 
#     summarise(max_cluster=max(n)) -> test[1,5] # maximum cluster number in one test
#   b %>% group_by(data,analysis,variable,direction,chromosome) %>% tally() %>% ungroup() %>% 
#     summarise(median_cluster=median(n)) -> test[1,6]# median of cluster numbers in one test
#   b %>% group_by(data,analysis,variable,direction,chromosome) %>% tally() %>% ungroup() %>% 
#     summarise(mean_cluster=mean(n)) -> test[1,7]#average cluster number in one test
#   b %>% group_by(data,analysis,variable,direction,chromosome) %>% tally() %>% ungroup() %>% 
#     summarise(sd_cluster=sd(n)) -> test[1,8]# standard deviation of cluster numbers
#   
#   #test with more than one cluster
#   b %>% group_by(data,analysis,variable,direction,chromosome) %>% tally() %>% ungroup() %>%
#     filter(n>1) %>% mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>% distinct(davdc) %>%
#     pull(davdc) -> vect2
#   
#   b %>% mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
#     filter(davdc %in% vect2) %>%
#     select(-davdc) %>%
#     separate(range,into = c("start","end"), sep=":", remove = F) %>%
#     group_by(data,analysis,variable,direction,chromosome) %>%
#     arrange(as.numeric(start)) %>%
#     mutate(distance=lead(as.numeric(start))-as.numeric(end)) %>% 
#     ungroup() %>%
#     filter(!is.na(distance)) %>%
#     select(distance)-> c 
#   c %>% summarise(max_distance=max(distance)) -> test[1,9]
#   c %>% summarise(min_distance=min(distance)) -> test[1,10]
#   c %>% summarise(sd_distance=sd(distance)) -> test[1,11]
#   c %>% summarise(average_distance=mean(distance)) -> test[1,12]
#   c %>% summarise(median_distance=median(distance)) -> test[1,13]
# 
#   
#       separate(start_end,into = c("start","end"), sep="-", remove = T) %>% 


  library(tidyverse)
  library(data.table)
  
LD <- c(0.8,0.85,0.90,0.95,0.99)
cM <- c(0.1,0.25,0.5,1,2,5,10) 

#to calculate distance between clusters
for (i in 1:length(LD)) {
  for (j in 1:length(cM)) {
    fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/convergent_clustering_"
                 ,LD[i],"_",cM[j],"_cM_summary"),header=T) -> a 
    a %>%
      group_by(data,analysis,variable,direction,chromosome) %>% 
      summarise(count=sum(N_convergent)) %>%
      ungroup() %>%
      filter(count>1) %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      distinct(davdc) %>%
      pull(davdc) -> vect
    
    #filter for the vector
    a %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      filter(davdc %in% vect) %>%
      select(-davdc) -> b

    rm(a)
    
    b %>% group_by(data,analysis,variable,direction,chromosome) %>% tally() %>% ungroup() %>%
      filter(n>1) %>% mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>% distinct(davdc) %>%
      pull(davdc) -> vect2
    
    b %>% mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      filter(davdc %in% vect2) %>%
      select(-davdc) %>%
      separate(range,into = c("start","end"), sep=":", remove = F) %>%
      group_by(data,analysis,variable,direction,chromosome) %>%
      arrange(as.numeric(start)) %>%
      mutate(distance=lead(as.numeric(start))-as.numeric(end)) %>% 
      ungroup() %>%
      filter(!is.na(distance)) %>%
      mutate(LD=LD[i],cM=cM[j]) %>%
      select(LD,cM,distance) %>% 
      fwrite("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/clustering_distance",sep = "\t",col.names= T,append = T)
    rm(b)
  }}

#to calculate cluster size
for (i in 1:length(LD)) {
  for (j in 1:length(cM)) {
    fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/convergent_clustering_"
                 ,LD[i],"_",cM[j],"_cM_summary"),header=T) -> a 
    a %>%
      group_by(data,analysis,variable,direction,chromosome) %>% 
      summarise(count=sum(N_convergent)) %>%
      ungroup() %>%
      filter(count>1) %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      distinct(davdc) %>%
      pull(davdc) -> vect
    
    #filter for the vector
    a %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      filter(davdc %in% vect) %>%
      select(-davdc) %>%
      filter(N_convergent>1) %>%
      mutate(LD=LD[i],cM=cM[j]) %>%
      select(LD,cM,size) %>%
      fwrite("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/clustering_size",sep = "\t",col.names= T,append = T)
      rm(a)
  }}

for (i in 1:length(LD)) {
  for (j in 1:length(cM)) {
    fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/convergent_clustering_"
                 ,LD[i],"_",cM[j],"_cM_summary"),header=T) -> a 
    a %>%
      group_by(data,analysis,variable,direction,chromosome) %>% 
      summarise(count=sum(N_convergent)) %>%
      ungroup() %>%
      filter(count>1) %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      distinct(davdc) %>%
      pull(davdc) -> vect
    
    #filter for the vector
    a %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      filter(davdc %in% vect) %>%
      select(-davdc) %>%
      group_by(data,analysis,variable,direction,chromosome) %>%
      tally() %>%
      ungroup() %>% 
      mutate(LD=LD[i],cM=cM[j]) %>%
      select(LD,cM,N_cluster=n) %>%
      fwrite("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/clustering_number",sep = "\t",col.names= T,append = T)
    rm(a)
  }}

#############uncorrected and baypass
library(tidyverse)
library(data.table)

LD <- c(0.8,0.85,0.90,0.95,0.99)
cM <- c(0.1,0.25,0.5,1,2,5,10) 

#to calculate distance between clusters
for (i in 1:length(LD)) {
  for (j in 1:length(cM)) {
    fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/convergent_clustering_"
                 ,LD[i],"_",cM[j],"_cM_summary"),header=T) %>% filter(analysis %in% c("baypass","GWAS_uncorrected"))-> a 
    a %>%
      group_by(data,analysis,variable,direction,chromosome) %>% 
      summarise(count=sum(N_convergent)) %>%
      ungroup() %>%
      filter(count>1) %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      distinct(davdc) %>%
      pull(davdc) -> vect
    
    #filter for the vector
    a %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      filter(davdc %in% vect) %>%
      select(-davdc) -> b
    
    rm(a)
    
    b %>% group_by(data,analysis,variable,direction,chromosome) %>% tally() %>% ungroup() %>%
      filter(n>1) %>% mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>% distinct(davdc) %>%
      pull(davdc) -> vect2
    
    b %>% mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      filter(davdc %in% vect2) %>%
      select(-davdc) %>%
      separate(range,into = c("start","end"), sep=":", remove = F) %>%
      group_by(data,analysis,variable,direction,chromosome) %>%
      arrange(as.numeric(start)) %>%
      mutate(distance=lead(as.numeric(start))-as.numeric(end)) %>% 
      ungroup() %>%
      filter(!is.na(distance)) %>%
      mutate(LD=LD[i],cM=cM[j]) %>%
      select(LD,cM,distance) %>% 
      fwrite("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/clustering_uncorrected_distance",sep = "\t",col.names= T,append = T)
    rm(b)
  }}

#to calculate cluster size
for (i in 1:length(LD)) {
  for (j in 1:length(cM)) {
    fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/convergent_clustering_"
                 ,LD[i],"_",cM[j],"_cM_summary"),header=T) %>% filter(analysis %in% c("baypass","GWAS_uncorrected")) -> a 
    a %>%
      group_by(data,analysis,variable,direction,chromosome) %>% 
      summarise(count=sum(N_convergent)) %>%
      ungroup() %>%
      filter(count>1) %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      distinct(davdc) %>%
      pull(davdc) -> vect
    
    #filter for the vector
    a %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      filter(davdc %in% vect) %>%
      select(-davdc) %>%
      filter(N_convergent>1) %>%
      mutate(LD=LD[i],cM=cM[j]) %>%
      select(LD,cM,size) %>%
      fwrite("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/clustering_uncorrected_size",sep = "\t",col.names= T,append = T)
    rm(a)
  }}

for (i in 1:length(LD)) {
  for (j in 1:length(cM)) {
    fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/convergent_clustering_"
                 ,LD[i],"_",cM[j],"_cM_summary"),header=T) %>% filter(analysis %in% c("baypass","GWAS_uncorrected")) -> a 
    a %>%
      group_by(data,analysis,variable,direction,chromosome) %>% 
      summarise(count=sum(N_convergent)) %>%
      ungroup() %>%
      filter(count>1) %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      distinct(davdc) %>%
      pull(davdc) -> vect
    
    #filter for the vector
    a %>%
      mutate(davdc=paste0(data,analysis,variable,direction,chromosome)) %>%
      filter(davdc %in% vect) %>%
      select(-davdc) %>%
      group_by(data,analysis,variable,direction,chromosome) %>%
      tally() %>%
      ungroup() %>% 
      mutate(LD=LD[i],cM=cM[j]) %>%
      select(LD,cM,N_cluster=n) %>%
      fwrite("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/combination/result/clustering_uncorreced_number",sep = "\t",col.names= T,append = T)
    rm(a)
  }}

#MAC
##########
library(tidyverse)
library(data.table)
fread("/Users/mojtabajahani/Downloads/clustering_distance",header = T,sep = "\t") %>% filter(LD!="LD") ->  clustering_distance
clustering_distance %>%
  mutate(distance=as.numeric(distance)) %>%
  group_by(LD,cM) %>%
  mutate(distance=as.numeric(distance)) %>%
  filter(distance<10000000) %>%
  tally() %>%
  ungroup() %>%
  #mutate(n=(n*-1)) %>%
  rename(under1=n) -> a

clustering_distance %>%
  mutate(distance=as.numeric(distance)) %>%
  group_by(LD,cM) %>%
  filter(distance>1000000) %>%
  tally() %>%
  ungroup() %>%
  rename(upper10=n)-> b

full_join(a,b) -> c

fread("/Users/mojtabajahani/Downloads/clustering_number",header = T,sep = "\t") %>%
  mutate(distance=as.numeric(N_cluster)) %>%
  filter(LD!="LD") %>%
  group_by(LD,cM) %>% 
  summarise(NC=sum(as.numeric(N_cluster))) %>% full_join(.,c) -> d


fread("/Users/mojtabajahani/Downloads/clustering_size",header = T,sep = "\t") %>%
  filter(LD!="LD") %>%
  mutate(distance=as.numeric(size)) -> clustering_size
  group_by(LD,cM) %>% 
  filter(size<10000000) %>%



%>% mutate(Not_clust=NC-5905) %>% mutate(ratunder=abs(under1)/Not_clust,ratup=upper10/Not_clust) -> d

####
clustering_distance %>% mutate(distance=as.numeric(distance),LD=as.numeric(LD),cM=as.character(cM)) -> clustering_distance
ggplot(clustering_distance, aes(LD, (distance/1000000),colour=cM )) + geom_smooth(method = lm, se = FALSE) 

clustering_distance %>% mutate(distance=as.numeric(distance),cM=as.numeric(cM),LD=as.character(LD)) -> clustering_distance
ggplot(clustering_distance, aes(cM, (distance/1000000),colour=LD )) + geom_smooth(method = lm, se = FALSE)


cor(clustering_distance$LD,clustering_distance$distance)

clustering_distance %>% filter(LD=="LD")

clustering_distance %>% mutate(cutoff=paste0("LD=",LD,"_","cM=",cM)) -> clustering_distance
clustering_distance %>% group_by(cutoff) %>% sample_n(30) %>% ungroup()  -> a 
a %>% distinct(cutoff) %>% sample_n(3) %>% pull(cutoff) -> cutt
a %>% filter(cutoff %in% cutt) -> b
ggplot(clustering_distance, aes(cutoff, (as.numeric(distance)/1000000))) + geom_boxplot() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


##########################
fread("/Users/mojtabajahani/Downloads/clustering_uncorrected_size",header = T,sep = "\t") %>%
  filter(LD!="LD") %>%
  mutate(size=as.numeric(size)) %>%
  mutate(cutoff=paste0("LD=",LD,"_","cM=",cM)) -> clustering_size
ggplot(clustering_size, aes(cutoff, (as.numeric(size)/1000000))) + 
  geom_boxplot() +
  labs(title="Cluster Size", y ="Mbp", x = "Cutoff") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


fread("/Users/mojtabajahani/Downloads/clustering_uncorrected_distance",header = T,sep = "\t") %>% 
  filter(LD!="LD") %>%
  mutate(distance=as.numeric(distance)) %>%
  mutate(cutoff=paste0("LD=",LD,"_","cM=",cM)) ->  clustering_distance
ggplot(clustering_distance, aes(cutoff, (as.numeric(distance)/1000000))) + 
  geom_boxplot() +
  labs(title="Distance Between Clusters", y ="Mbp", x = "Cutoff") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


fread("/Users/mojtabajahani/Downloads/clustering_uncorreced_number",header = T,sep = "\t") %>%
  filter(LD!="LD") %>%
  mutate(N_cluster=as.numeric(N_cluster)) %>%
  mutate(cutoff=paste0("LD=",LD,"_","cM=",cM)) ->  clustering_number
ggplot(clustering_number, aes(cutoff,N_cluster)) + 
  geom_boxplot() +
  labs(title="Number of Clusters per chromosome", y ="Count", x = "Cutoff") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
 

clustering_distance %>%
  mutate(gr=cut(distance, breaks= c(0,1000000,5000000,10000000,50000000,500000000))) %>%
  group_by(LD,cM,gr) %>%
  tally() %>%
  ungroup() %>%
  mutate(cutoff=paste0("LD=",LD,"_","cM=",cM)) -> underMB
underMB$gr <- factor(underMB$gr, levels = rev(levels(underMB$gr)))
ggplot(underMB,aes(x=cutoff,y=n,fill=gr)) +
  geom_bar(stat="identity") +
  labs(title="Count of distance between clusters", y ="Count", x = "Cutoff") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_discrete(name = "Distance Category", labels = c( ">50 Mbp","10-50 Mbp", "5-10 Mbp","1-5 Mbp","<1 Mbp"))

clustering_size %>%
  mutate(gr=cut(size, breaks= c(0,1000000,5000000,10000000,50000000,500000000))) %>%
  group_by(LD,cM,gr) %>%
  tally() %>%
  ungroup() %>%
  mutate(cutoff=paste0("LD=",LD,"_","cM=",cM)) -> underMB_size
underMB_size$gr <- factor(underMB_size$gr, levels = rev(levels(underMB_size$gr)))
ggplot(underMB_size,aes(x=cutoff,y=n,fill=gr)) +
  geom_bar(stat="identity") +
  labs(title="Count of clusters based on size", y ="Count", x = "Cutoff") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_discrete(name = "Distance Category", labels = c( ">50 Mbp","10-50 Mbp", "5-10 Mbp","1-5 Mbp","<1 Mbp"))

fread("/Users/mojtabajahani/Downloads/clustering_uncorrected_size",header = T,sep = "\t") %>%
  filter(LD!="LD") %>%
  mutate(size=as.numeric(size))-> a

cor(as.numeric(a$cM),as.numeric(a$size))

fread("/Users/mojtabajahani/Downloads/clustering_uncorreced_number",header = T,sep = "\t") %>%
  filter(LD!="LD") %>%
  mutate(N_cluster=as.numeric(N_cluster))-> b 
cor(as.numeric(b$LD),as.numeric(b$N_cluster))
