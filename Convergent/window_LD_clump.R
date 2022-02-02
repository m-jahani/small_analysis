library(foreach)
library(doParallel)
library(tidyverse)
library(chunked)
library(data.table)

options(scipen = 9999)
fread("/moonriseNFS/Null_W_Result/files/convergent_2binning",sep = "\t",header=T) %>%
  distinct(window5) %>%
  mutate(convergence="convergent")-> Convergent_window

fread("/data/users/mjahani/GWAS_conver/LD_cluster/all_id_gene_window_number_maf.03",sep = "\t",header = T) %>%
  filter(data== "env") %>%
  select(Taxa,ID,window5) %>%
  left_join(.,Convergent_window) %>%
  mutate(convergence=replace_na(convergence,"non_convergent"))-> ID_convergent

ID_convergent %>% distinct(Taxa) %>% pull(Taxa) -> TAXA
ID_convergent %>% distinct(convergence) %>% pull(convergence) -> conver
ID_convergent %>% separate(ID,into = c("chr","pos"), sep=":", remove = T) %>% distinct(chr) %>% pull(chr) -> CHR

registerDoParallel(cores=45)
for (i in 1:length(conver)) {
  for (j in 1:length(TAXA)) {
      foreach(k=1:length(CHR), .combine='rbind', .errorhandling='stop') %dopar% {
    ID_convergent %>%
      filter(convergence==conver[i]) %>%
      distinct(Taxa,ID) %>%
      separate(ID,into = c("chrom","pos"), sep=":", remove = T) %>%
      filter(Taxa==TAXA[j]) %>%
      filter(chrom==CHR[k]) %>%
      select(chrom,pos) %>%
      fwrite(.,file = paste0("/data/users/mjahani/GWAS_conver/LD_cluster/id_list/IDs_",TAXA[j],"_",conver[i],"_",CHR[k]),sep = "\t",col.names = F)
  }}}
rm(ID_convergent)
c("Annuus.tranche90.snp.env","Argophyllus.tranche90.snp.gwas","Petiolaris.tranche90.snp.petfal","Petiolaris.tranche90.snp.petpet") -> VCF

for (i in 1:length(conver)) {
  foreach(j=1:length(TAXA), .combine='rbind', .errorhandling='stop') %dopar% {
    foreach(k=1:length(CHR), .combine='rbind', .errorhandling='stop') %dopar% {
    system(paste0("vcftools --gzvcf /moonriseNFS/soudi/null_W/ld/",
                  VCF[j],
                  ".90.bi.remappedHa412HO.vcf.gz --positions /data/users/mjahani/GWAS_conver/LD_cluster/id_list/IDs_",
                  TAXA[j],
                  "_",
                  conver[i],
                  "_",
                  CHR[k],
                  " --recode --recode-INFO-all --out /data/users/mjahani/GWAS_conver/LD_cluster/VCFs/",
                  TAXA[j],
                  "_",
                  conver[i],
                  "_",
                  CHR[k]))
  }}}

for (i in 1:length(conver)) {
  foreach(j=1:length(TAXA), .combine='rbind', .errorhandling='stop') %dopar% {
    foreach(k=1:length(CHR), .combine='rbind', .errorhandling='stop') %dopar% {
  system(paste0("/home/mjahani/bin/plink --vcf /data/users/mjahani/GWAS_conver/LD_cluster/VCFs/",
                TAXA[j],
                "_",
                conver[i],
                "_",
                CHR[k],
                ".recode.vcf  --allow-extra-chr --double-id --recode --out /data/users/mjahani/GWAS_conver/LD_cluster/ped/",
                TAXA[j],
                "_",
                conver[i],
                "_",
                CHR[k]))
  }}}

system(paste0("rm /data/users/mjahani/GWAS_conver/LD_cluster/VCFs/*"))



for (i in 1:length(conver)) {
  for (j in 1:length(TAXA)) {
    foreach(k=1:length(CHR), .combine='rbind', .errorhandling='stop') %dopar% {
  fread(paste0("/data/users/mjahani/GWAS_conver/LD_cluster/ped/",
               TAXA[j],
               "_",
               conver[i],
               "_",
               CHR[k],
               ".map"),header = F) %>%
    mutate(V5=paste0(V1,":",V4)) %>%
    select(V1,V5,V3,V4) %>%
    fwrite(paste0("/data/users/mjahani/GWAS_conver/LD_cluster/ped/",
                  TAXA[j],
                  "_",
                  conver[i],
                  "_",
                  CHR[k],
                  ".map"),col.names = F,sep = "\t")
  }}}
  
# fread("/data/users/mjahani/GWAS_conver/LD_cluster/all_id_gene_window_number_maf.03",sep = "\t",header = T) %>%
#   distinct(ID,window5) %>%
#   separate(ID,into = c("chrom","pos"), sep=":", remove = F) %>%
#   select(chrom,ID,window5)-> ID_WIN


TAXA <- c("Annuus","Argophyllus","petfal","petpet")
conver <- c("convergent")
CHR <- c("Ha412HOChr01","Ha412HOChr02","Ha412HOChr03","Ha412HOChr04","Ha412HOChr05",
         "Ha412HOChr06","Ha412HOChr07","Ha412HOChr08","Ha412HOChr09","Ha412HOChr10",
         "Ha412HOChr11","Ha412HOChr12","Ha412HOChr13","Ha412HOChr14","Ha412HOChr15",
         "Ha412HOChr16","Ha412HOChr17")




registerDoParallel(cores=5)
for (i in 1:length(conver)) {
  for (j in 1:length(TAXA)) {
    foreach(k=1:length(CHR), .combine='rbind', .errorhandling='stop') %dopar% {
      fread(paste0("/data/users/mjahani/LD_cluster/ped/",
                   TAXA[j],
                   "_",
                   conver[i],
                   "_",
                   CHR[k],
                   ".map"),header = F) %>%
        mutate(V1=gsub("Ha412HOChr","",V1)) %>%
        mutate(V2=gsub("Ha412HOChr","",V2)) %>%
        fwrite(paste0("/data/users/mjahani/LD_cluster/ped/",
                      TAXA[j],
                      "_",
                      conver[i],
                      "_",
                      CHR[k],
                      ".map"),col.names = F,sep = "\t")}}}


start_time_total <- Sys.time() 
for (i in 1:length(conver)) {
  for (j in 1:length(TAXA)) {
    for (k in 1:length(CHR)) {
  system(paste0("/home/mjahani/bin/plink --file /data/users/mjahani/LD_cluster/ped/",
                TAXA[j],
                "_",
                conver[i],
                "_",
                CHR[k],
                "  /data/users/mjahani/LD_cluster/paiwise_LD/",
                TAXA[j],
                "_",
                conver[i],
                "_",
                CHR[k]))
      
      
      system(paste0("wc -l /data/users/mjahani/LD_cluster/paiwise_LD/",
                    TAXA[j],
                    "_",
                    conver[i],
                    "_",
                    CHR[k],
                    ".ld > /data/users/mjahani/LD_cluster/paiwise_LD/",
                    TAXA[j],
                    "_",
                    conver[i],
                    "_",
                    CHR[k],
                    ".count &"))
      
      
      setwd("/data/users/mjahani/LD_cluster/paiwise_LD/temp/")

      system(paste0("split -l 500000  /data/users/mjahani/LD_cluster/paiwise_LD/",
                TAXA[j],
                     "_",
                     conver[i],
                     "_",
                     CHR[k],".ld"))
      awk '{print $3,$6,$7}' /home/mjahani/scratch/pauline/Annuus.ann_fullsam.tranche90_snps_bi_AN50_AF99_MAF03_CHR10.ld | split -l 500000
      # system(paste0("pigz -dc  /data/users/mjahani/LD_cluster/paiwise_LD/",
      #               TAXA[j],
      #               "_",
      #               conver[i],
      #               "_",
      #               CHR[k],
      #               ".ld.gz | awk '{print $3,$6,$7}' | parallel --block 999M --pipe ' pigz > /data/users/mjahani/LD_cluster/paiwise_LD/temp/part_{#}.gz'"))
      
      
      
      library(foreach)
      library(doParallel)
      library(tidyverse)
      library(data.table)
      
      fread("/home/mjahani/scratch/pauline/SNP_WIN",sep = "\t",header = F) %>%
        rename(ID=V1,window5=V2) %>%
        mutate(ID=gsub("Ha412HOChr10:","",ID)) %>%
        mutate(window5=gsub("Ha412HOChr10:","",window5)) %>%
        select(ID,window5) -> ID_WIN_CHR
       
  
      
      list.files(path = "/home/mjahani/scratch/pauline/SPLIT/") -> file_list
        registerDoParallel(cores=10)
        foreach(z=1:length(file_list), .combine ='rbind', .errorhandling='stop') %dopar% {
        fread(paste0("/home/mjahani/scratch/pauline/SPLIT/",file_list[z]),
                       header = F,showProgress=F) %>%
              select(SNP_A=V1,SNP_B=V2,R2=V3) %>%
              mutate(SNP_A=gsub("Ha412HOChr10:","",SNP_A)) %>%
              mutate(SNP_B=gsub("Ha412HOChr10:","",SNP_B)) %>%
              left_join(.,rename(ID_WIN_CHR,SNP_A=ID,WIN_A=window5)) %>%
              filter(!is.na(WIN_A)) %>%
              left_join(.,rename(ID_WIN_CHR,SNP_B=ID,WIN_B=window5)) %>%
              mutate(DW=paste(WIN_A,WIN_B,sep = "_"),R2=as.numeric(R2)) %>%
              select(DW,R2) %>%
              group_by(DW) %>%
              summarise(LD=mean(R2),snp_count=n()) %>%
              ungroup() %>%
              separate(DW,into = c("window_A","window_B"), sep="_", remove = T) %>%
              select(window_A,window_B,LD,snp_count) %>%
              fwrite(paste0("/home/mjahani/scratch/pauline/CHR10.LD_WINDOW"),sep = "\t",col.names= F,append = T)}
    

result <- NULL
for (i in c("xaa","xab","xac","xad","xae","xaf")) {
fread(paste0("/home/mjahani/scratch/pauline/SPLT2/",i),header = F) %>%
  rename(window_A=V1,window_B=V2,LD=V3,snp_count=V4) %>%
  mutate(mult=(LD*snp_count)) %>%
  group_by(window_A,window_B) %>%
  summarise(n_snp=sum(snp_count),multip=sum(mult)) %>%
  ungroup() %>% 
  mutate(mean_LD=(multip/n_snp)) %>%
  select(window_A,window_B,LD=mean_LD,snp_count=n_snp) %>%
    rbind(.,result)-> result }
47306381694 data1
47306381694
  fwrite(paste0("/home/mjahani/scratch/pauline/CHR10.LD_WINDOW_FINAL"),sep = "\t",col.names= F,append = T)
  
/home/mjahani/scratch/pauline/SPLT2/xaa
/home/mjahani/scratch/pauline/SPLT2/xab
/home/mjahani/scratch/pauline/SPLT2/xac
/home/mjahani/scratch/pauline/SPLT2/xad
/home/mjahani/scratch/pauline/SPLT2/xae
/home/mjahani/scratch/pauline/SPLT2/xaf
  
  read.table("/home/mjahani/scratch/pauline/CHR10.LD_WINDOW_FINAL",header = F) -> N
  sum()
  
  
  fread("/data/users/mjahani/scaffolding/PK/REP/PK_assembly_chr_only_rep1.mapinfo") %>%
    rename(QNAME=V1,RNAME=V2,POS=V3,MAPQ=V4,RNEXT=V5,PNEXT=V6) -> data

  data %>% 
    mutate(chr=ifelse(RNAME==RNEXT,"intra","inter")) %>%
    mutate(dist=abs(PNEXT-POS)) -> data
  
    
data.frame(taxa=TAXA[j],
                          convergent=conver[i],
                          chr=CHR[k],
                          count1=N[1,1]-1,
                          count2=(sum(LD_WINDOW$V5))) %>%
  mutate(test=(count1-count2)) %>%
  select(taxa,convergent,chr,count1,count2,test) %>%
  fwrite(paste0("/data/users/mjahani/LD_cluster/pair_wind/counting"),sep = "\t",col.names= F,append = T)
  

rm(LD_WINDOW,LD_WINDOW_result,N)

system(paste0("rm /data/users/mjahani/LD_cluster/pair_wind/",
             TAXA[j],
             "_",
             conver[i],
             "_",
             CHR[k],".LD_WINDOW"))

system(paste0("rm /data/users/mjahani/LD_cluster/paiwise_LD/",
              TAXA[j],
              "_",
              conver[i],
              "_",
              CHR[k],
              ".count"))
}}}
end_time_total <- Sys.time()
end_time_total - start_time_total

library(tidyverse)
library(data.table)
options(scipen = 999999)
TAXA<-c("petpet","petfal","Annuus","Argophyllus")
limit<-c(105450000,108695000,110760000,112460000)
for (i in 1:length(TAXA) ) {
fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/",TAXA[i],"_non_convergent_full"),header = F) %>%
    select(window_A=V2,window_B=V3,mean_LD=V4) -> a
  a %>%
    mutate(window_A=gsub("Ha412HOChr\\d\\d:","",window_A),window_B=gsub("Ha412HOChr\\d\\d:","",window_B)) -> a
  a %>%
    mutate(window_A=gsub("-\\d+","",window_A),window_B=gsub("-\\d+","",window_B)) -> a
  a %>%
    mutate(distance=abs(as.numeric(window_A)-as.numeric(window_B))) %>%
    select(-window_A,-window_B) %>%
    #filter(distance<=limit[i]) %>%
    select(distance,mean_LD) %>%
    fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/",TAXA[i],"_distance_meanLD_limitless"),
           sep = "\t",col.names= F)
  fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/",TAXA[i],"_distance_meanLD_limitless"),header = F) %>%
  select(distance=V1) %>%
  mutate(Taxa=TAXA[i]) %>%
    group_by(Taxa,distance) %>%
    tally() %>%
    ungroup() %>%
  fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/10k_random_all_limitless"),
                  sep = "\t",col.names= F,append = T)}

fread("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/10k_random_all_limitless",header = F) %>%
  select(Taxa=V1,distance=V2,count=V3)-> tenk_random_all

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
    summarise(threshold=quantile(mean_LD, 0.95)) %>%
    ungroup() %>%
    fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/Taxa_distance_threshold_limitless"),
           sep = "\t",col.names= F,append = T)
  
  fread(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/",TAXA[i],"_distance_meanLD_limitless"),header = F) %>%
    select(distance=V1,mean_LD=V2) %>%
    mutate(Taxa=TAXA[i])  %>%
    left_join(.,tenk_random_all) %>%
    filter(count<10000) %>%
    group_by(Taxa,distance) %>%  
    summarise(threshold=quantile(mean_LD, 0.95)) %>%
    ungroup() %>%
    fwrite(paste0("/data/mojtaba/GWAS_conver/LD_cluster/pair_wind/final/Taxa_distance_threshold_limitless"),
           sep = "\t",col.names= F,append = T)}
  


fread(paste0("/data/users/mjahani/LD_cluster/pair_wind/Taxa_distance_threshold_limitless_0.9"),header = F) %>%
select(Taxa=V1,distance=V2,threshold=V3) -> Taxa_distance_threshold



fread(paste0("/data/users/mjahani/LD_cluster/pair_wind/convergent_result"),header = F) %>%
  select(Taxa=V1,window_A=V2,window_B=V3,mean_LD=V4) %>%
  separate(window_A,into = c("chrom_A","start_end_A"), sep=":", remove = F) %>%
  separate(start_end_A,into = c("start_A","end_A"), sep="-", remove = T) %>%
  separate(window_B,into = c("chrom_B","start_end_B"), sep=":", remove = F) %>%
  separate(start_end_B,into = c("start_B","end_B"), sep="-", remove = T) %>%
  mutate(distance=abs(as.numeric(start_A)-as.numeric(start_B))) %>%
  select(Taxa,window_A,window_B,mean_LD,distance) %>%
  left_join(.,Taxa_distance_threshold) %>%
  mutate(significant_LD=(mean_LD>=threshold))-> convergent_result

fwrite(convergent_result,"/data/users/mjahani/LD_cluster/pair_wind/convergent_result_sig_limitless_0.9",sep = "\t")
#######
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
options(scipen = 9999)
fread("/data/users/mjahani/LD_cluster/pair_wind/convergent_result_sig_limitless_0.9",sep = "\t") -> convergent_result


fread("/data/users/mjahani/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") -> map
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

distinct(rbind(rename(distinct(convergent_result,window_A),window=window_A),rename(distinct(convergent_result,window_B),window=window_B)),window) %>%
  separate(window,into = c("chr","start_end"), sep=":", remove = F) %>%
  separate(start_end,into = c("start","end"), sep="-", remove = T) %>% 
  mutate(start=as.numeric(start),end=as.numeric(end)) %>%
  mutate(pos=(((end-start)/2)+start)) %>%
  select(window,chr,pos) -> conv_pos


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
rm(conv_pos,recombination)

window_recombination %>% 
  mutate(G_position=(((pos-start)*recomb_rate)+cM)) %>%
  select(window,G_position) -> window_recombination

convergent_result %>% distinct(Taxa) %>% pull(Taxa)->TAX
taxa_conv_pos <- NULL

for (i in 1:length(TAX)) {
distinct(rbind(rename(distinct(filter(convergent_result,Taxa==TAX[i]),Taxa,window_A),window=window_A),rename(distinct(filter(convergent_result,Taxa==TAX[i]),Taxa,window_B),window=window_B)),Taxa,window) %>%
  separate(window,into = c("chr","start_end"), sep=":", remove = F) %>%
  separate(start_end,into = c("start","end"), sep="-", remove = T) %>% 
  mutate(start=as.numeric(start),end=as.numeric(end)) %>%
  mutate(pos=(((end-start)/2)+start)) %>%
  select(Taxa,window,chr,pos) %>%
  rbind(.,taxa_conv_pos) -> taxa_conv_pos}

fread("/data/users/mjahani/inversion.bed",sep = "\t") -> bed  
bed %>% filter(species=="Petiolaris") %>% mutate(Taxa="petpet") %>% select(-species) ->pet
bed %>% filter(species=="Petiolaris") %>% mutate(Taxa="petfal") %>% select(-species) ->fal
rbind(pet,fal)-> petfal
bed %>% filter(species!="Petiolaris") %>% rename(Taxa=species) %>% rbind(.,petfal) -> bed
rm(pet,fal,petfal)

registerDoParallel(cores=40)
window_inversion<-foreach(i=1:nrow(bed), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(taxa_conv_pos,
                       chr == as.character(bed[i,1]) &
                         pos >= as.numeric(bed[i,2]) &
                         pos <= as.numeric(bed[i,3]) &
                         Taxa == as.character(bed[i,5])),
                window,pos,Taxa),
         inversion=as.character(bed[i,4]))}

window_inversion %>%
  mutate(tax_win=paste0(Taxa,"_",window)) %>%
  group_by(tax_win) %>%
  tally() %>% 
  ungroup() %>% 
  filter(n>1) %>% 
  distinct(tax_win) %>% 
  pull(tax_win) -> vector

window_inversion %>%
  mutate(tax_win=paste0(Taxa,"_",window)) %>%
  filter(tax_win %in% vector) %>%
  distinct(inversion)

# inversion
# 1  pet14.01
# 2  pet14.02
# 3  pet16.01
# 4  pet16.02
# 5  pet17.03
# 6  pet17.04

window_inversion %>%
  mutate(inversion=ifelse(inversion=="pet14.02","pet14.01",inversion)) %>%
  mutate(inversion=ifelse(inversion=="pet16.02","pet16.01",inversion)) %>%
  mutate(inversion=ifelse(inversion=="pet17.04","pet17.03",inversion)) -> window_inversion
  
window_inversion %>% 
  distinct(window,pos,Taxa,inversion) %>%
  select(Taxa,window,inversion)-> window_inversion



convergent_result %>%
  left_join(.,rename(window_recombination,window_A=window,G_position_A=G_position)) %>%
  left_join(.,rename(window_recombination,window_B=window,G_position_B=G_position)) %>%
  mutate(Genetic_distance=abs(G_position_A-G_position_B)) %>%
  mutate(GD_significant=ifelse(Genetic_distance<=5,"TRUE","FALSE")) %>%
  select(-G_position_A,-G_position_B) %>%
  left_join(.,rename(window_inversion,window_A=window,inversion_A=inversion)) %>%
  left_join(.,rename(window_inversion,window_B=window,inversion_B=inversion)) %>% 
  mutate(inversion_A=replace_na(inversion_A,"not_inversion")) %>%
  mutate(inversion_B=replace_na(inversion_B,"not_inversion"))  -> convergent_result
fwrite(convergent_result,"/data/users/mjahani/LD_cluster/pair_wind/convergent_result_sig_limitless_0.9_LD_GD_INV",sep = "\t") 

# convergent_result %>%   
#   filter(inversion_A!="not_inversion") %>%
#   filter(inversion_B!="not_inversion") %>%
#   mutate(inv_iden_significant=(inversion_A==inversion_B)) %>%
#   filter(inv_iden_significant=="TRUE")-> inversion_sig
# 
# left_join(convergent_result,inversion_sig) %>%
#   mutate(inv_iden_significant=replace_na(inv_iden_significant,"FALSE"))-> convergent_result

##############################
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
options(scipen = 9999)

fread("/data/users/mjahani/LD_cluster/pair_wind/convergent_result_sig_limitless_0.95_LD_GD_INV",header = T) ->   convergent_result

convergent_result %>%
  filter(distance!=0) %>%
  mutate(LD_GD_sig=ifelse(significant_LD=="TRUE" & GD_significant=="TRUE","TRUE","FALSE")) -> data
# data %>%
#   mutate(LD_GD_INV_sig=ifelse(significant_LD=="TRUE" & inv_iden_significant=="TRUE","TRUE","FALSE"))-> data


  rbind(data,
        select(data,
               Taxa,
               window_A=window_B,
               window_B=window_A,
               mean_LD,
               distance,
               threshold,
               significant_LD,
               Genetic_distance,
               GD_significant,
               inversion_A,
               inversion_B,
               LD_GD_sig)) %>%
    separate(window_A,into = c("chrom","start_end_A"), sep=":", remove = F) %>%
    select(Taxa,chrom,window_A,window_B,LD_GD_sig)-> convergent_result_sig_double
rm(data)  


fread("/data/users/mjahani/Null_W_Result/files/convergent_2binning",sep = "\t",header=T) %>%
  select(type,analysis,Taxa1,Taxa2,variable,window5,Emp_p) -> convergent_window

# fread("/data/users/mjahani/Null_W_Result/files/convergent_2binning",sep = "\t",header=T) %>%
#   select(type,analysis,Taxa1,Taxa2,variable,window5,Emp_p) %>%
#   select(-analysis) %>%
#   distinct(type,Taxa1,Taxa2,variable,window5,.keep_all = T) %>%
#   mutate(analysis="merge") %>%
#   select(type,analysis,Taxa1,Taxa2,variable,window5,Emp_p)-> convergent_window


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
                
                 convergent_result_sig_double %>%
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
                   fwrite("/data/users/mjahani/LD_cluster/cluster/convergent_clustering_merge_0.95_5cM",sep = "\t",col.names= F,append = T)
            }#chrom
          }#taxa1
        }#taxa2
      }#variable
    }#analysis
   }#type

fread("/data/users/mjahani/LD_cluster/cluster/convergent_clustering_0.9_5cM",header=F) %>%
  select(type=V1,analysis=V2,variable=V3,Taxa1=V4,Taxa2=V5,chrom=V6,window_name=V7,cluster=V8) %>%
  separate(window_name,into = c("chr","start_end"), sep=":", remove = F) %>%
  select(-chr) %>%
  separate(start_end,into = c("start","end"), sep="-", remove = T) %>% 
  mutate(start=as.numeric(start),end=as.numeric(end)) %>%
  group_by(type,analysis,variable,Taxa1,Taxa2,chrom,cluster) %>%
  summarise(N_convergent=n(),cluster_start=min(start),cluster_end=max(end)) %>%
  ungroup() %>%
  mutate(range=paste0(as.character(cluster_start),":",as.character(cluster_end))) %>%
  mutate(size=(as.numeric(cluster_end)-as.numeric(cluster_start))+1) %>%
  mutate(direction=paste0(Taxa1,"_",Taxa2)) %>%
  select(data=type,analysis,variable,direction,chromosome=chrom,range,size,N_convergent) %>%
  fwrite("/data/users/mjahani/LD_cluster/cluster/convergent_clustering_0.9_5cM_summary",sep = "\t",col.names= T,append = T)

fread("/data/users/mjahani/LD_cluster/cluster/convergent_clustering_summary",header=T) -> a


fread("/home/mjahani/pauline/LD_calculation/Annuus.ann_fullsam.tranche90_snps_bi_AN50_AF99_MAF03_CHR10.map") %>%
  select(CHR=V1,
         POS=V4) -> SNPS

190802806/5000
38161*5000
seq(1,190802806,5000)-> a

data.table(a) -> b

b %>%
  rename(START=a) %>%
  mutate(CHR="Ha412HOChr10") %>%
  mutate(END=START+4999) %>%
  mutate(WIN_ID=paste0(CHR,":",START,"-",END)) %>% 
  select(CHR,
         START,
         END,
         WIN_ID) -> BED
  
registerDoParallel(cores=48)
SNP_GENE <- foreach(i=1:nrow(BED), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(SNPS,
                       CHR == as.character(BED[i,1]),
                       (POS >= as.numeric(BED[i,2]) &
                          POS <= as.numeric(BED[i,3]))
  ),
  CHR,POS),
  Gene_ID=as.character(BED[i,4]))}
