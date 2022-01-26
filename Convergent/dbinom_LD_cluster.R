library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)


dbinom_result_all <- NULL
    for (i in c("Baypass","GWAS_corrected","GWAS_uncorrected","Spearman")) {
      list.files(path = paste0("/data/users/mjahani/dbinom_LD_cluster/dbinom_recombination_bins_scores/",i))-> file_list
      for (j in 1:length(file_list)) { 
        tryCatch(fread(paste0("/data/users/mjahani/dbinom_LD_cluster/dbinom_recombination_bins_scores/",i,"/",file_list[j]),
                       header = T ,stringsAsFactors = FALSE), error=function(e) NULL) -> data
        if (!is.null(data)) {
          data %>%
            gather(species,dbinom,1:4) %>%
            mutate(window=paste0(chrom,":",start_window,"-",end_window)) %>%
            mutate(analysis=i) %>%
            select(analysis,variable,species,window,dbinom) -> tmp1
          tmp1 %>%
            filter(!is.na(dbinom)) %>%
            group_by(analysis,variable,species) %>%
            summarise(threshold=quantile(dbinom, 0.95)) %>%
            ungroup() %>%
            full_join(.,tmp1) %>%
            mutate(is_top_window=ifelse(dbinom > threshold,"YES", "NO")) %>%
            filter(is_top_window=="YES") %>%
            select(analysis,variable,species,window,dbinom) %>%
            rbind(.,dbinom_result_all) -> dbinom_result_all
        }
      }
    }
fwrite(dbinom_result_all,"/data/users/mjahani/dbinom_LD_cluster/dbinom_recombination_bins_scores/dbinom_result_all",sep = "\t",col.names = T)
#some of the top dibinom snps in annuus(from GWA analysis) did not exist in GEA's VCF. So I had to use to VCFs (GEA and GWA) for Annuus 
c("Argophyllus","petpet","petfal") -> SPE
fread("/data/users/mjahani/all_id_gene_window_number_maf.03",sep = "\t",header = T) -> all_ID
for (i in 1:length(SPE)) {
 all_ID %>%
  filter(Taxa == SPE[i]) %>%
  distinct(ID) %>%
    separate(ID,into = c("chrom","pos"), sep=":", remove = T) %>%
    select(chrom,pos) -> mid_data
    mid_data %>% distinct(chrom) %>% pull(chrom) -> CHR
    registerDoParallel(cores=48)
    foreach(j=1:length(CHR), .combine='rbind', .errorhandling='stop') %dopar% {
      mid_data %>%
        filter(chrom==CHR[j]) %>%
        fwrite(.,file = paste0("/data/users/mjahani/dbinom_LD_cluster/id_list/IDs_",SPE[i],"_",CHR[j],"_env"),sep = "\t",col.names = F)
    }
    rm(mid_data,CHR)
    }

c("gwa","env")-> DATA
for (i in 1:length(data)) {
  all_ID %>%
    filter(Taxa == "Annuus") %>%
    filter(data == DATA[i]) %>%
    distinct(ID) %>%  
    separate(ID,into = c("chrom","pos"), sep=":", remove = T) %>%
    select(chrom,pos) -> mid_data
  mid_data %>% distinct(chrom) %>% pull(chrom) -> CHR
  registerDoParallel(cores=48)
  foreach(j=1:length(CHR), .combine='rbind', .errorhandling='stop') %dopar% {
    mid_data %>%
      filter(chrom==CHR[j]) %>%
      fwrite(.,file = paste0("/data/users/mjahani/dbinom_LD_cluster/id_list/IDs_Annuus_",CHR[j],"_",DATA[i]),sep = "\t",col.names = F)
    
  }
  rm(mid_data,CHR)
}        


list.files(path = paste0("/data/users/mjahani/dbinom_LD_cluster/id_list/"))-> ID_file_list

c(rep("Annuus.tranche90.snp.env",17),
  rep("Argophyllus.tranche90.snp.gwas",17),
  rep("Petiolaris.tranche90.snp.petfal",17),
  rep("Petiolaris.tranche90.snp.petpet",17)) -> VCF

ID_file_list[!grepl(".gwa",ID_file_list)] -> ID_file_list_env
registerDoParallel(cores=48)
foreach(i=1:length(ID_file_list_env), .combine='rbind', .errorhandling='stop') %dopar% {
  system(paste0("vcftools --gzvcf /moonriseNFS/wild_gwas/",
                VCF[i],
                ".90.bi.remappedHa412HO.vcf.gz --positions /data/users/mjahani/dbinom_LD_cluster/id_list/",
                ID_file_list_env[i],
                " --recode --recode-INFO-all --out /data/users/mjahani/dbinom_LD_cluster/VCFs/",
                substring(ID_file_list_env[i],5,last = 100L)))
}


ID_file_list[grepl(".gwa",ID_file_list)] -> ID_file_list_gwa
registerDoParallel(cores=48)
foreach(i=1:length(ID_file_list_gwa), .combine='rbind', .errorhandling='stop') %dopar% {
  system(paste0("vcftools --gzvcf /moonriseNFS/wild_gwas/",
                "Annuus.tranche90.snp.gwas.90.bi.remappedHa412HO.vcf.gz --positions /data/users/mjahani/dbinom_LD_cluster/id_list/",
                ID_file_list_gwa[i],
                " --recode --recode-INFO-all --out /data/users/mjahani/dbinom_LD_cluster/VCFs/",
                substring(ID_file_list_gwa[i],5,last = 100L)))
}

list.files(path = paste0("/data/users/mjahani/dbinom_LD_cluster/VCFs/"))-> VCF_file_list
VCF_file_list[grepl(".vcf",VCF_file_list)]-> VCF_file_list
foreach(i=1:length(VCF_file_list), .combine='rbind', .errorhandling='stop') %dopar% {
  system(paste0("/home/mjahani/bin/plink --vcf /data/users/mjahani/dbinom_LD_cluster/VCFs/",
                VCF_file_list[i],
                " --allow-extra-chr --double-id --recode --out /data/users/mjahani/dbinom_LD_cluster/ped/",
                VCF_file_list[i]))
}

list.files(path = paste0("/data/users/mjahani/dbinom_LD_cluster/ped/"))-> ped_file_list
ped_file_list[grepl(".map",ped_file_list)]-> ped_file_list
foreach(i=1:length(ped_file_list), .combine='rbind', .errorhandling='stop') %dopar% {
  fread(paste0("/data/users/mjahani/dbinom_LD_cluster/ped/",
               ped_file_list[i]),header = F) %>%
    mutate(V5=gsub("Ha412HOChr","",V1)) %>%
    mutate(V6=paste0("C",V5)) %>%
    mutate(V7=paste0(V6,":",V4)) %>%
    select(V6,V7,V3,V4) %>%
    fwrite(paste0("/data/users/mjahani/dbinom_LD_cluster/ped/",
                  ped_file_list[i]),col.names = F,sep = "\t")
}

Annuus_Ha412HOChr02_gwa.recode.vcf.ld

list.files(path = paste0("/data/users/mjahani/dbinom_LD_cluster/ped/"))-> ped_file_list
ped_file_list[grepl(".map",ped_file_list)]-> ped_file_list
#Annuus_Ha412HOChr02_gwa.recode.vcf.map
for (i in 1:length(ped_file_list)) {
  
system(paste0("/home/mjahani/bin/plink --file /data/users/mjahani/dbinom_LD_cluster/ped/",
              gsub(".map","",ped_file_list[i]),
              " --r2 --ld-window-kb 999999999 --ld-window 999999999 --ld-window-r2 0 --allow-extra-chr --out /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
              gsub(".map","",ped_file_list[i])))


system(paste0("wc -l /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
              gsub(".map","",ped_file_list[i]),
              ".ld > /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
              gsub(".map","",ped_file_list[i]),
              ".count &"))

# system(paste0("du --block-size=1G /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#               gsub(".map","",ped_file_list[i]),
#               ".ld > /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#               gsub(".map","",ped_file_list[i]),
#               ".size"))

# fread(paste0("/data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#              gsub(".map","",ped_file_list[i]),
#              ".size")) %>%
#   pull(V1) -> size
  
# if (as.numeric(size) < 2000) {
#   setwd("/data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp/")
#   system(paste0("split -l 500000  /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#                 gsub(".map","",ped_file_list[i]),".ld"))
# }else{
#   Sys.time() -> start2
  system(paste0("pigz /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
                gsub(".map","",ped_file_list[i]),".ld"))
  
  system(paste0("pigz -dc  /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
                gsub(".map","",ped_file_list[i]),
                ".ld.gz |  parallel --block 999M --pipe ' pigz > /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp3/",gsub(".recode.vcf.map","",ped_file_list[i]),"_part_{#}.gz'"))
  

#   Sys.time() -> end2
#   end2-start2 -> duration
# }
###
# Sys.time() -> start1
# setwd("/data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp/")
# system(paste0("split -l 500000  /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#               gsub(".map","",ped_file_list[i]),".ld"))
# Sys.time() -> end1
# end1-start1 -> duration1 #Time difference of 3.113045 hours #1102GB #1102GB
# 
# 
# Sys.time() -> start2
# 
# system(paste0("awk '{print $3,$6,$7}' /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/", gsub(".map","",ped_file_list[i]),
#               ".ld | parallel --block 500M --pipe ' cat > /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp1/part_{#}'"))
# Sys.time() -> end2
# end2-start2 -> duration2 # Time difference of 7.08393 hours #500GB #1102GB
# 
# Sys.time() -> start3
# system(paste0("pigz /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#               gsub(".map","",ped_file_list[i]),".ld"))
# 
# system(paste0("pigz -dc  /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#               gsub(".map","",ped_file_list[i]),
#               ".ld.gz | awk '{print $3,$6,$7}' | parallel --block 500M --pipe ' cat > /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp2/part_{#}'"))
# Sys.time() -> end3
# end3-start3 -> duration3 #Time difference of 8.813148 hours #152GB #500
# 
# Sys.time() -> start4
# system(paste0("pigz -dc  /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#               gsub(".map","",ped_file_list[i]),
#               ".ld.gz | parallel --block 500M --pipe ' cat > /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp3/part_{#}'"))
# Sys.time() -> end4
# end4-start4 -> duration4  #Time difference of 2.818966 hours #152GB #1102GB
# 
# Sys.time() -> start5
# system(paste0("pigz -dc  /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#               gsub(".map","",ped_file_list[i]),
#               ".ld.gz | awk '{print $3,$6,$7}' | parallel --block 999M --pipe ' pigz > /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp/part_{#}.gz'"))
# Sys.time() -> end5
# end5-start5 -> duration5#Time difference of 6.616865 hours #152GB #108GB
# 
# Sys.time() -> start6
# system(paste0("pigz -dc  /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#               gsub(".map","",ped_file_list[i]),
#               ".ld.gz | parallel --block 500M --pipe ' pigz > /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp1/part_{#}.gz'"))
# Sys.time() -> end6
# end6-start6 -> duration6 #Time difference of 2.399034 hours #152GB #152GB
# 
# system("pigz -d /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/Annuus_Ha412HOChr02_gwa.recode.vcf.ld.gz")
# 
# Sys.time() -> start7
# system(paste0("pigz /data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
#               gsub(".map","",ped_file_list[i]),".ld"))
# Sys.time() -> end7
# end7-start7 -> duration7 #43 mins


fread("/data/users/mjahani/all_id_gene_window_number_maf.03",sep = "\t",header = T) %>%
  distinct(ID,window5) %>%
  separate(ID,into = c("chrom","pos"), sep=":", remove = F) %>%
  select(chrom,ID,window5) %>%
  filter(chrom==gsub("(_gwa.recode.vcf.map|_env.recode.vcf.map)","",gsub("[A-z]{6,11}_","",ped_file_list[i]))) %>%
  mutate(ID=gsub("Ha412HOChr","C",ID)) %>%
  select(ID,window5) -> ID_WIN_CHR

list.files(path = "/data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp3/") -> file_list

registerDoParallel(cores=40)
foreach(z=1:length(file_list), .combine ='rbind', .errorhandling='stop') %dopar% {
  fread(paste0("/data/users/mjahani/dbinom_LD_cluster/paiwise_LD/temp1/",file_list[z]),
        header = F,showProgress=F) %>%
    select(SNP_A=V3,SNP_B=V6,R2=V7) %>%
    left_join(.,rename(ID_WIN_CHR,SNP_A=ID,WIN_A=window5)) %>%
    filter(!is.na(WIN_A)) %>%
    left_join(.,dplyr::rename(ID_WIN_CHR,SNP_B=ID,WIN_B=window5)) %>%
    mutate(DW=paste(WIN_A,WIN_B,sep = "_"),R2=as.numeric(R2)) %>%
    select(DW,R2) %>%
    group_by(DW) %>%
    summarise(LD=mean(R2),snp_count=n()) %>%
    ungroup() %>%
    separate(DW,into = c("window_A","window_B"), sep="_", remove = T) %>%
    mutate(tmp=ped_file_list[i]) %>%
    mutate(tmp=gsub(".recode.vcf.map","",tmp)) %>%
    separate(tmp,into= c("Taxa","ch","data"), sep= "_",remove= T) %>%
    select(data,Taxa,window_A,window_B,LD,snp_count) %>%
    fwrite(paste0("/data/users/mjahani/dbinom_LD_cluster/pair_wind/",
                  gsub(".recode.vcf.map","",ped_file_list[i]),".LD_WINDOW"),sep = "\t",col.names= F,append = T)}

system(paste0("rm /data/users/mjahani/LD_cluster/paiwise_LD/*.ld"))
#system(paste0("rm /data/users/mjahani/LD_cluster/paiwise_LD/*.ld.gz"))
system(paste0("rm /data/users/mjahani/LD_cluster/paiwise_LD/*log"))
system(paste0("rm /data/users/mjahani/LD_cluster/paiwise_LD/*.nosex"))
system(paste0("rm /data/users/mjahani/LD_cluster/paiwise_LD/temp/*"))

rm(file_list,ID_WIN_CHR)


fread(paste0("/data/users/mjahani/dbinom_LD_cluster/pair_wind/",
             gsub(".recode.vcf.map","",ped_file_list[i]),".LD_WINDOW"),header = F) -> LD_WINDOW 
LD_WINDOW %>%
  rename(data=V1,Taxa=V2,window_A=V3,window_B=V4,LD=V5,snp_count=V6) %>%
  mutate(mult=(LD*snp_count)) %>%
  group_by(Taxa,window_A,window_B) %>%
  summarise(n_snp=sum(snp_count),multip=sum(mult)) %>%
  ungroup() %>% 
  mutate(mean_LD=(multip/n_snp)) %>%
  select(Taxa,window_A,window_B,LD=mean_LD,snp_count=n_snp) %>%
  fwrite(paste0("/data/users/mjahani/dbinom_LD_cluster/pair_wind/window_pairwise_LD_result"),sep = "\t",col.names= F,append = T)


              gsub(".map","",ped_file_list[i]),
".count &"

read.table(paste0("/data/users/mjahani/dbinom_LD_cluster/paiwise_LD/",
                  gsub(".map","",ped_file_list[i]),".count"),header = F)-> N


data.frame(data=ped_file_list[i],
           count1=N[1,1]-1,
           count2=(sum(LD_WINDOW$V6))) %>%
  mutate(test=(count1-count2)) %>%
  select(data,count1,count2,test) %>%
  fwrite(paste0("/data/users/mjahani/dbinom_LD_cluster/pair_wind/counting"),sep = "\t",col.names= F,append = T)


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



######################################
dbinom_result_all %>%
  mutate(data=ifelse(analysis=="Spearman" | analysis=="Baypass" ,"env","gwa")) %>%
  distinct(data,species,window) %>%
  mutate(is_convergent="YES") -> spe_win

fread("/data/users/mjahani/all_id_gene_window_number_maf.03",sep = "\t",header = T) %>%
  select(data,species=Taxa,ID,window=window5) %>%
  left_join(.,spe_win) %>%
  filter(!is.na(is_convergent)) %>%
  select(data,species,window,ID) ->  ID_convergent

for file in *,gz
do
mv "$file" "${file/,gz/.gz}"
done
