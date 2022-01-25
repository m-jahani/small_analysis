################################################################################################################################################################
#########################################################LD base window clustering. May 2020####################################################################
################################################################################################################################################################


# The clustering procedure conducts a stepwise process as follows: at the initial step, it takes the first convergent window from each chromosome and estimates
# the Linkage Disequilibrium (LD) with the second convergent window. At the following step, the procedure moves a window forward and calculates LD between the
# second and third windows. The process continues moving forward step by step as long as these conditions are met in each step; 1) the calculated LD is
# significantly different (P < .05) from the null distribution, and 2) the distance between the two windows is not more than 5 cM. Immediately after a step has
# not met either of the two conditions, a cluster forms with all previous steps’ windows. These windows are then subsetted from the data, and the clustering
# algorithm starts from the first step with the remaining windows. This process continues until all the convergent windows have been clustered.
# In each step of the procedure, PLINK (reference) was used to calculate LD (r2) among all pairwise combinations of the SNPs (MAF < 3%) within the two windows.
# Later on, all the calculated r2 values were summarized by averaging to estimate the LD level among the two windows of the step. The null distribution in each
# step was constructed with LD estimations among 10,000 random, non-convergent window pairs. The distance between each window pair matches the physical distance
# of the two convergent windows of the corresponding step.
# 
# 
# Mojtaba Jahani 
# mojtaba.jahani@hotmail.com
#################################################################Softwares######################################################################################

#VCFtools
#PLINK

#############################################################Files and directories##############################################################################
directory = "/data/mojtaba/test/"#The main directory to save results and temporary files
window_IDs  = "/data/mojtaba/example/window_IDs" #List of the SNPs co-locating with the windows need to be clustered  (MAF < 3 or 5%)
BG_IDs  = "/data/mojtaba/example/background_IDs" #The list for the rest of SNPs in the VCF (MAF < 3 or 5%), which do not exist in the first file ("window_IDs")
VCF = "/data/mojtaba/example/example.vcf" #The VCF
map_address = "/data/mojtaba/example/map.txt" #Genetic map
window_list = "/data/mojtaba/example/list_window" #The list of windows for LD clustering
ID_window_file = "/data/mojtaba/example/ID_window"#all the SNPs and the corresponding windows

################################################################Theresholds######################################################################################

threshold_cM = 5 #Distance threshold
quant = 0.95 #LD cutoff
nthreads = 48 #number of threads for parallel computation 

#################################################################Packages#######################################################################################

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)

##############################################################Directories#######################################################################################

setwd(directory)
system(paste0("mkdir ","ID_list VCFs ped paiwise_LD Ld_gz pair_window cluster"))
ID_list = paste0(directory,"ID_list/")
VCF_out = paste0(directory,"VCFs/")
PLINK_in = paste0(directory,"ped/")
LD_dir = paste0(directory,"paiwise_LD/")
LD_GZ = paste0(directory,"Ld_gz/")
LD_win = paste0(directory,"pair_window/")
cluster_dir = paste0(directory,"cluster/")

##############################################################Data preparation##################################################################################
#LD calculation outputs are going to be huge files. To make the analysis, RAM and storage-friendly, all the calculations will be carried out in chromosome level.

#generates a vector of chromosome IDs
fread(paste0(window_IDs),
      header  =  F) %>% 
  rename(chr = V1,
         pos = V2) %>%
  distinct(chr) %>%
  pull(chr) -> CHR


#Divides ID list to single chromosome files
IDs <- c(window_IDs,
         BG_IDs)
names <- c("window",
           "BG")
for (i in 1:length(IDs)) {
  for (j in 1:length(CHR)) {
    fread(paste0(IDs[i]),
          header  =  F) %>%
      rename(chr = V1,
             pos = V2) %>%
      filter(chr == CHR[j]) %>%
      select(chr,
             pos) %>%
      fwrite(.,
             file  =  paste0(ID_list,
                             "IDs_",
                             names[i],
                             "_",
                             CHR[j]),
             sep  =  "\t",
             col.names  =  F)
  }
}

#breaks the main VCF to single chromosome VCFs
registerDoParallel(cores=nthreads) 
  for (i in 1:length(IDs)) {  
foreach(j = 1:length(CHR),
        .combine = 'rbind',
        .errorhandling = 'stop') %dopar% {
      system(paste0("vcftools --vcf ",
                    VCF,
                    " --positions ",
                    ID_list,
                    "IDs_",
                    names[i],
                    "_",
                    CHR[j],
                    " --recode --recode-INFO-all --out ",
                    VCF_out,
                    names[i],
                    "_",
                    CHR[j]))
        }
    }


#convert VCFs to an input format for PLINK (ped format)
registerDoParallel(cores=nthreads) 
for (i in 1:length(IDs)) {  
  foreach(j = 1:length(CHR),
          .combine = 'rbind',
          .errorhandling = 'stop') %dopar% {
      system(paste0("/data/mojtaba/example/plink --vcf ",
                    VCF_out,
                    names[i],
                    "_",
                    CHR[j],
                    ".recode.vcf  --allow-extra-chr --double-id --recode --out ",
                    PLINK_in,
                    names[i],
                    "_",
                    CHR[j]))
          }
  }

#Remove all VCFs
system(paste0("rm -r ",
              VCF_out))

#This part is not needed if the VCFs have the ID column information. This script can produce ID column in *.map files
registerDoParallel(cores=nthreads) 
for (i in 1:length(IDs)) {
    foreach(j = 1:length(CHR),
            .combine = 'rbind',
            .errorhandling = 'stop') %dopar% {
      fread(paste0(PLINK_in,
                   names[i],
                   "_",
                   CHR[j],
                   ".map"),header  =  F) %>%
        mutate(V5 = paste0(V1,
                         ":",
                         V4)) %>%
        select(V1,
               V5,
               V3,
               V4) %>%
        fwrite(paste0(PLINK_in,
                      names[i],
                      "_",
                      CHR[j],
                      ".map"),col.names  =  F,
               sep  =  "\t")
            }
  }


#load data
fread(ID_window_file,
      sep  =  "\t",header  =  F) %>%
  select(ID = V1,
         window5 = V2)-> ID_WIN_CHR

##############################################################Data preparation##################################################################################
#This section of the script calculates LD (r2) between all window combinations in each chromosome

#This part calculates LD among all pairwise combinations of the SNPs for each chromosome.
for (i in 1:length(IDs)) {  
 for (j in 1:length(CHR)) {
   
   system(paste0("/data/mojtaba/example/plink --file ",
                    PLINK_in,
                    names[i],
                    "_",
                    CHR[j],
                    " --r2 --ld-window-kb 999999999 --ld-window 999999999 --ld-window-r2 0 --allow-extra-chr --out ",
                    LD_dir,
                    names[i],
                    "_",
                    CHR[j]))
   
#This part was only designed to compare the number of SNP to test if any small part of the script failed to run.
#count number of pairwise LD calculation
system(paste0("wc -l ",
              LD_dir,
              names[i],
              "_",
              CHR[j],
              ".ld > ",
              LD_dir,
              names[i],
              "_",
              CHR[j],
              ".count "))


# In some cases, even one chromosome LD calculation file is much bigger than our RAM.
# This part of the script breaks down LD calculation too small files to be loaded on RAM.
setwd(LD_GZ)#directory to save the small files
system(paste0("split -l 50000 ",#number of line per small file
              LD_dir,
              names[i],
              "_",
              CHR[j],
              ".ld"))

#Calculates LD between two windows by averaging pairwise SNP LDs.
list.files(path  =  LD_GZ) -> file_list
registerDoParallel(cores=nthreads) 
foreach(z=1:length(file_list),
        .combine='rbind',
        .errorhandling='stop') %dopar% {
  #for (z in 1:length(file_list)) {
  fread(paste0(LD_GZ,
               file_list[z]),
        header  =  F,
        showProgress  =  F) %>%
    select(SNP_A = V3,
           SNP_B = V6,
           R2 = V7) %>%
    left_join(.,
              rename(ID_WIN_CHR,
                     SNP_A=ID,
                     WIN_A=window5)) %>%
    filter(!is.na(WIN_A)) %>%
    left_join(.,
              rename(ID_WIN_CHR,
                     SNP_B=ID,
                     WIN_B=window5)) %>%
    mutate(DW = paste(WIN_A,
                      WIN_B,
                      sep="_"),
           R2 = as.numeric(R2)) %>%
    select(DW,
           R2) %>%
    group_by(DW) %>%
    summarise(LD = mean(R2),
              snp_count = n()) %>%
    ungroup() %>%
    separate(DW,
             into  =  c("window_A",
                        "window_B"),
             sep = "_",
             remove  =  T) %>%
    select(window_A,
           window_B,
           LD,
           snp_count) %>%
    fwrite(paste0(LD_win,
                  names[i],
                  "_",
                  CHR[j],
                  ".LD_WINDOW"),
           sep="\t",
           col.names =F,
           append=T)
          }


#Remove extra files
system(paste0("rm ",
              LD_dir,
              "*.ld"))
system(paste0("rm ",
              LD_dir,
              "*.log"))
system(paste0("rm ",
              LD_dir,
              "*.nosex"))
system(paste0("rm ",
              LD_GZ,
              "*"))
rm(file_list)


#Merge window pair LD results
fread(paste0(LD_win,
             names[i],
             "_",
             CHR[j],
             ".LD_WINDOW"),header  =  F) -> LD_WINDOW 
LD_WINDOW %>%
  rename(window_A = V1,
         window_B = V2,
         LD = V3,
         snp_count = V4) %>%
  mutate(mult = (LD*snp_count)) %>%
  group_by(window_A,
           window_B) %>%
  summarise(n_snp = sum(snp_count),
            multip = sum(mult)) %>%
  ungroup() %>% 
  mutate(mean_LD = (multip/n_snp)) %>%
  select(window_A,
         window_B,
         LD = mean_LD,
         snp_count = n_snp) %>%
  fwrite(paste0(LD_win,
                names[i],
                "_result"),
         sep  =  "\t",
         col.names =  F,
         append  =  T)


#This part was only designed to compare the number of SNP to test if any small part of the script failed to run.
#The last column af the "counting" file must be always zero, otherwise something wend wrong in the calculation
read.table(paste0(LD_dir,
                  names[i],
                  "_",
                  CHR[j],
                  ".count"),
           header  =  F)-> N
data.frame(data = names[i],
           chr = CHR[j],
           count1 = N[1,1]-1,
           count2 = (sum(LD_WINDOW$V4))) %>%
  mutate(test = (count1-count2)) %>%
  select(data,
         chr,
         count1,
         count2,
         test) %>%
  fwrite(paste0(LD_win,
                "counting"),
         sep  =  "\t",
         col.names =  F,
         append  =  T)
rm(LD_WINDOW,
   N)

#delete files
system(paste0("rm ",LD_win,
              names[i],
              "_",
              CHR[j],
              ".LD_WINDOW"))

system(paste0("rm ",
              LD_dir,
              names[i],
              "_",
              CHR[j],
              ".count"))
 }
  }


##############################################################Null distribution##################################################################################
#This section of the script calculates a Null distribution based on background SNPs/windows.

#reports the distance between two windows in column 1 and LD in column 2
  fread(paste0(LD_win,
               "BG_result"),
        header  =  F) %>%
    select(window_A = V1,
           window_B = V2,
           mean_LD = V3) %>%
    separate(window_A,
             into  =  c("chrom_A",
                        "start_end_A"),
             sep = ":",
             remove  =  F) %>%
    separate(start_end_A,
             into  =  c("start_A","end_A"),
             sep = "-",
             remove  =  T) %>%
    separate(window_B,
             into  =  c("chrom_B",
                        "start_end_B"),
             sep = ":",
             remove  =  F) %>%
    separate(start_end_B,
             into  =  c("start_B",
                        "end_B"),
             sep = "-",
             remove  =  T) %>%
    mutate(distance = abs(as.numeric(start_A)-as.numeric(start_B))) %>%
    select(-window_A,
           -window_B) %>%
    select(distance,
           mean_LD) %>%
    fwrite(paste0(LD_win,
                  "BG_distance_meanLD"),
           sep  =  "\t",
           col.names =  F)

#counts the number of windows pairs for each possible distance
  fread(paste0(LD_win,
               "BG_distance_meanLD"),
        header  =  F) %>%
    select(distance = V1) %>%
    group_by(distance) %>%
    tally() %>%
    ungroup() %>%
    fwrite(paste0(LD_win,"distance_count"),
           sep  =  "\t",
           col.names =  F,
           append  =  T)
fread(paste0(LD_win,"distance_count"),header  =  F) %>%
  select(distance = V1,
         count = V2)-> tenk_random_all


#samples random 10,000 LD for each window-pair distance (when at least 10,000 window pair exist for the distance), 
#and calculates the quantile value (for example 0.95)
  fread(paste0(LD_win,
               "BG_distance_meanLD"),
        header  =  F) %>%
    select(distance = V1,mean_LD = V2) %>%
    left_join(.,tenk_random_all) %>%
    filter(count >= 10000) %>%
    group_by(distance) %>%
    sample_n(10000) %>%
    ungroup() %>%
    group_by(distance) %>%  
    summarise(threshold = quantile(mean_LD,
                                   quant)) %>%
    ungroup() %>%
    fwrite(paste0(LD_win,
                  "BG_threshold"),
           sep  =  "\t",
           col.names =  F,
           append  =  T)
  
  #Exactly like the previous step, except it uses all window-pairs for each distance (when window pairs for the distance are less than 10,000),
  #the quantile values will be added to the same file as the previous step.
  fread(paste0(LD_win,
               "BG_distance_meanLD"),
        header  =  F) %>%
    select(distance = V1,mean_LD = V2) %>%
    left_join(.,
              tenk_random_all) %>%
    filter(count<10000) %>%
    group_by(distance) %>%  
    summarise(threshold = quantile(mean_LD,
                                   quant)) %>%
    ungroup() %>%
    fwrite(paste0(LD_win,
                  "BG_threshold"),
           sep  =  "\t",
           col.names =  F,
           append  =  T)
  
  #load threshold file (result of the two last steps)
  fread(paste0(LD_win,
               "BG_threshold"),
        header  =  F) %>%
    select(distance = V1,
           threshold = V2) -> BG_threshold
  
#This step uses the threshold data to see if LD between the window pairs are significant (larger than the threshold)  
  fread(paste0(LD_win,
               "window_result"),
        header  =  F) %>%
    select(window_A = V1,
           window_B = V2,
           mean_LD = V3) %>%
    separate(window_A,
             into  =  c("chrom_A",
                        "start_end_A"),
             sep = ":",
             remove  =  F) %>%
    separate(start_end_A,
             into  =  c("start_A",
                        "end_A"),
             sep = "-",
             remove  =  T) %>%
    separate(window_B,
             into  =  c("chrom_B",
                        "start_end_B"),
             sep = ":", remove  =  F) %>%
    separate(start_end_B,into  =  c("start_B",
                                    "end_B"),
             sep = "-",
             remove  =  T) %>%
    mutate(distance = abs(as.numeric(start_A)-as.numeric(start_B))) %>%
    select(window_A,
           window_B,
           mean_LD,
           distance) %>%
    left_join(.,
              BG_threshold) %>%
    mutate(significant_LD = (mean_LD >= threshold))-> convergent_result

##########################################################################genetic map##################################################################################
#load genetic map
  fread(map_address,
        header = T) -> map

#It define bins and calculate the recombination rate
  map %>%
  group_by(chr) %>%
    arrange(pos,
            .by_group  =  TRUE) %>%
    mutate(delta_cM = cM-lag(cM,
                             default  =  first(cM)),
           delta_bp = pos-lag(pos,
                              default  =  first(pos)),
           start = lag(pos)) %>% 
    ungroup() %>%
    mutate(recomb_rate = delta_cM/delta_bp) %>%
    select(chr,
           start,
           end = pos,
           recomb_rate) %>%
    filter(!is.na(start)) %>%
    left_join(.,
              rename(map,
                     start = pos)) -> recombination
  rm(map)
#It calculates the center position of each window as "window position"
  distinct(rbind(rename(distinct(convergent_result,
                                 window_A),
                        window = window_A),rename(distinct(convergent_result,
                                                           window_B),
                                                  window = window_B)),
           window) %>%
    separate(window,
             into  =  c("chr",
                        "start_end"),
             sep = ":",
             remove  =  F) %>%
    separate(start_end,
             into  =  c("start","end"),
             sep = "-", remove  =  T) %>% 
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>%
    mutate(pos = (((end-start)/2)+start)) %>%
    select(window,
           chr,
           pos) -> conv_pos
  
  
  #It finds the corresponding recombination bin for each window, and assign the Genetic distance (cM)
  registerDoParallel(cores=nthreads) 
  window_recombination<-foreach(i = 1:nrow(recombination),
                                .combine = 'rbind',
                                .errorhandling = 'stop') %dopar% { 
    mutate(select(filter(conv_pos,
                         chr  ==  as.character(recombination[i,1])  &
                           pos >=  as.numeric(recombination[i,2])&
                           pos <=  as.numeric(recombination[i,3])),
                  window,pos),
           start = as.numeric(recombination[i,2]),
           recomb_rate = as.numeric(recombination[i,4]),
           cM = as.numeric(recombination[i,5]))}
  rm(conv_pos,recombination)

  #It calculates the genetic position based on the location between genetic map markers
  window_recombination %>% 
    mutate(G_position = (((pos-start)*recomb_rate)+cM)) %>%
    select(window,G_position) -> window_recombination
  
 #It generates a list of windows that needed to be clustered (plus chromosome and position) 
  distinct(rbind(rename(distinct(convergent_result,window_A),
                        window = window_A),
                 rename(distinct(convergent_result,
                                 window_B),
                        window = window_B)),window) %>%
    separate(window,
             into  =  c("chr",
                               "start_end"),
             sep = ":",
             remove  =  F) %>%
    separate(start_end,
             into  =  c("start",
                        "end"),
             sep = "-",
             remove  =  T) %>% 
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>%
    mutate(pos = (((end-start)/2)+start)) %>%
    select(window,
           chr,
           pos)  -> window_pos
  
 #It calculates the genetic distance between all window pairs and based on genetic distance, and LD of each window determines if two windows are linked  
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
    mutate(GD_significant = ifelse(Genetic_distance <=  threshold_cM ,
                                   "TRUE",
                                   "FALSE")) %>%
    select(-G_position_A,-G_position_B) %>%
    filter(distance != 0) %>%
    mutate(LD_GD_sig = ifelse(significant_LD == "TRUE" & GD_significant == "TRUE",
                              "TRUE",
                              "FALSE"))-> convergent_result
  
  #It generates a table that shows if two windows are linked (table of linkage)
  rbind(convergent_result,
        select(convergent_result,
               window_A = window_B,
               window_B = window_A,
               mean_LD,
               distance,
               threshold,
               significant_LD,
               Genetic_distance,
               GD_significant,
               LD_GD_sig)) %>%
    separate(window_A,
             into  =  c("chrom",
                        "start_end_A"),
             sep = ":",
             remove  =  F) %>%
    select(chrom,
           window_A,
           window_B,
           LD_GD_sig)-> convergent_result_sig_double
  rm(convergent_result)  

  ##########################################################################LD clustering##################################################################################
  
#loads the window list  
  fread(window_list,
        sep  =  "\t",
        header = F) %>%
    rename(window = V1) %>%
    separate(window,
             into  =  c("chrom",
                        "start_end"),
             sep = ":",
             remove  =  F) %>%
    separate(start_end,
             into  =  c("start","end"),
             sep = "-",
             remove  =  T) %>%
    mutate(start = as.numeric(start)) %>%
    select(chrom,
           start,
           window) -> convergent_window
  
  #list of chromosomes
  convergent_window %>% 
    distinct(chrom) %>% 
    pull(chrom)-> CHR 
  
  for (i in 1:length(CHR)) {
    #filter for chromosome i
    convergent_window %>%
      filter(chrom == CHR[i]) %>%
      arrange(start) %>%
      select(-chrom,
             -start) %>%
      distinct(window)-> temporary
    #filter the table of linkage for chromosome i
    convergent_result_sig_double %>%
      filter(chrom == CHR[i]) %>%
      select(window_A,
             window_B,
             LD_GD_sig) -> con_sig
    #empty cluster variable
    cluster <- NULL
    #fills the cluster variable up with LD clustering steps(a window pair per step)
    if (1<nrow(temporary)){
      for (o in 1:(nrow(temporary)-1)) {
        rbind(data.frame(window_A = temporary[o,1],
                         window_B = temporary[o+1,1],
                         cycle = o),
              cluster)-> cluster
      }#rows
    }else{rbind(data.frame(window_A = temporary[1,1],
                           window_B = "-",
                           cycle = 1),cluster)-> cluster}
    #calculates if each of the two windows is linked
    left_join(cluster,
              con_sig) %>%
      mutate(LD_GD_sig = replace_na(LD_GD_sig,
                                    "FALSE")) %>% arrange(cycle) %>% mutate(clustering = 0) -> cluster
    #It forms clusters
    for (p in 1:nrow(cluster)) {
      cluster[p,5] <- ifelse(cluster[p,3] == 1,1,
                             cluster[p,5] <- ifelse(cluster[p,5] ==  0 & cluster[p-1,4] == "TRUE" & cluster[p,4] == "TRUE" ,
                                                    cluster[p-1,5],
                                                    as.numeric(cluster[p-1,5])+1))
    }#cluster rows
    #save the clustering information
    cluster %>% filter(LD_GD_sig == "TRUE") %>%
      gather(windows,
             window_name,
             1:2) %>% 
      select(window_name,
             clustering) %>% distinct(window_name,
                                      .keep_all  =  T) -> true_cluster
    temporary %>%
      filter(!window %in% pull(distinct(true_cluster,
                                        window_name),
                               window_name)) %>%
      mutate(clustering = window) %>%
      select(window_name = window,
             clustering) %>%
      rbind(.,true_cluster) %>%
      mutate(cluster = group_indices(.,clustering)) %>%
      mutate(chrom = CHR[i]) %>%
      select(chrom,
             window_name,
             cluster) %>%
      fwrite(paste0(cluster_dir,
                    "cluster"),
             sep  =  "\t",
             col.names =  F,
             append  =  T)
  }#chrom
  

  #summarize cluster data
  fread(paste0(cluster_dir,
               "cluster"),
        header = F) %>%
    select(chrom = V1,
           window_name = V2,
           cluster = V3) %>%
    separate(window_name,
             into  =  c("chr",
                        "start_end"),
             sep = ":",
             remove  =  F) %>%
    select(-chr) %>%
    separate(start_end,
             into  =  c("start",
                        "end"),
             sep = "-",
             remove  =  T) %>% 
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>%
    group_by(chrom,
             cluster) %>%
    summarise(N_windows = n(),
              cluster_start = min(start),
              cluster_end = max(end)) %>%
    ungroup() %>%
    mutate(range = paste0(as.character(cluster_start),
                          ":",
                          as.character(cluster_end))) %>%
    mutate(size = (as.numeric(cluster_end)-as.numeric(cluster_start))+1) %>%
    select(chromosome = chrom,
           range,
           size,
           N_windows) %>%
    fwrite(paste0(cluster_dir,"cluster_summary"),
           sep  =  "\t",
           col.names =  T)
  
  