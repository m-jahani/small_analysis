library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)

inversion_directory <- "/data/users/mjahani/FST/petiolaris/"
VCF <-"/moonriseNFS/wild_gwas/Petiolaris.tranche90.snp.gwas.90.bi.remappedHa412HO.vcf.gz"
fread(paste0("/data/users/mjahani/FST/Petfal.individual.GEA.list"),header = F) %>%
  mutate(taxa="petfal") %>% 
  select(sample=V1,taxa)-> petfal_list #petfal sample list
fread(paste0("/data/users/mjahani/FST/Petpet.individual.list"),header = F) %>%
  mutate(taxa="petpet") %>% 
  select(sample=V1,taxa)-> petpet_list
rbind(petfal_list,petpet_list) -> petpet_petfal #petpet sample list
rm(petfal_list,petpet_list)

# fread("/data/users/mjahani/Ha412HO_inv.v3.inversions.regions.v1.txt") -> inversions.regions
# fread("/data/users/mjahani/Ha412HO_inv.v3.txt") %>%
#   mutate(chr=ifelse(chr>=10,paste0("Ha412HOChr",chr),paste0("Ha412HOChr0",chr))) %>% 
#   mutate(mds=paste0(direction,mds)) %>%
#   select(-species) %>%
#   left_join(inversions.regions,.) %>%
#   filter(species!="Niveus") %>%
#   select(species,chr,mds,sv_name) -> name_key

# bcftools query -l /moonriseNFS/wild_gwas/Petiolaris.tranche90.snp.gwas.90.bi.remappedHa412HO.vcf.gz > indiv_list
list.files(path = inversion_directory) -> file_list
sample_inversion <- NULL
for (i in 1:length(file_list)) {
  fread(paste0(inversion_directory,file_list[i]),sep="\t",header = T) %>%
    select(sample,triangle_genotype) %>%
    left_join(petpet_petfal,.) %>%
    mutate(inversion=file_list[i]) %>%
    mutate(inversion=gsub("Ha412HO_inv.v3.pcasites.","",inversion)) %>%
    mutate(inversion=gsub(".genotypes.txt","",inversion)) %>%
    rbind(.,sample_inversion) -> sample_inversion #genotype each inversion for each sample
  }

fread("/data/users/mjahani/Ha412HO_inv.v3.inversions.regions.v1.txt") -> inversions.regions
fread("/data/users/mjahani/Ha412HO_inv.v3.petinvfreq.txt") -> Petiolaris.frequency

inversions.regions %>%
  filter(species=="Petiolaris") %>%
  select(-species) %>%
  mutate(species="PetFal") %>%
  left_join(.,Petiolaris.frequency)%>%
  filter(as.numeric(freq) > 0) %>%
  select(-species) %>%
  mutate(species="petfal") %>%
  select(species,start,end,chr,mds) -> inversions.regions.petfal #cordination of each inversion in petfal

inversions.regions %>%
  filter(species=="Petiolaris") %>%
  select(-species) %>%
  mutate(species="PetPet") %>%
  left_join(.,Petiolaris.frequency)%>%
  filter(as.numeric(freq) > 0) %>%
  select(-species) %>%
  mutate(species="petpet") %>%
  select(species,start,end,chr,mds) -> inversions.regions.petpet #cordination of each inversion in petpet

# rbind(inversions.regions.petfal,inversions.regions.petpet) %>%
#  mutate(tinv=paste0(species,":",chr,".",mds))  %>%
#   distinct(tinv) %>%
#   pull(tinv) -> TINV
# 
# sample_inversion %>%
#   mutate(tinv=paste0(taxa,":",inversion)) %>%
#   filter(tinv %in% TINV) %>%
#   select(-tinv) -> sample_inversion

rbind(inversions.regions.petfal,inversions.regions.petpet) %>%
  mutate(inversion=paste0(chr,".",mds)) %>%
  distinct(species,inversion) %>%
  group_by(inversion) %>%
  tally() %>%
  ungroup() %>%
  rename(N_inversion=n) -> inversion_count #number of inversions in petpet and petfal (1= inversion is private to one species, 2= inversion is common to both species)

  right_join(inversion_count,sample_inversion) -> sample_inversion
  
  # inversion_count %>% 
  #   filter(N_inversion==1) %>%
  #   distinct(inversion) %>%
  #   pull(inversion) -> INV1
  #   rbind(inversions.regions.petfal,inversions.regions.petpet) %>%
  #     mutate(inversion=paste0(chr,".",mds)) %>%
  #     distinct(species,inversion) %>%
  #     filter(inversion %in% INV1)
# sample_inversion %>% 
#   spread(inversion,triangle_genotype) -> sample_inversion

#bcftools query -f '%CHROM:%POS\n' /moonriseNFS/wild_gwas/Petiolaris.tranche90.snp.gwas.90.bi.remappedHa412HO.vcf.gz >  /data/users/mjahani/FST/Petiolaris_ID_list 
    
    
vanil_50K<-read.table(file ="/data/users/mjahani/bins_good_5000",header = FALSE)
vanil_50K$V6<-NULL
vanil_50K$V1<-paste("window",vanil_50K$V1,sep = "_")
vanil_50K$V2 <- gsub("chr", "", vanil_50K$V2)
colnames(vanil_50K)<-c("window","chrom","start","end","width")
vanil_50K %>% 
  select(-window) %>%
  mutate(chrom=as.numeric(chrom)) %>%
  mutate(chr=ifelse(chrom>=10,paste0("Ha412HOChr",chrom),paste0("Ha412HOChr0",chrom))) %>%
  mutate(window=paste0(chr,":",start,"-",end)) %>%
  select(window,chr,start,end) -> vanil_50K # all possible 5 kb windows in sunflower genome

fread("/data/users/mjahani/FST/Petiolaris_ID_list",header=F) %>%
  separate(V1,into = c("chr","pos"),sep = ":",remove = F) %>%
  mutate(pos=as.numeric(pos)) %>%
  select(ID=V1,chr,pos)-> ID_list # list od snps in petiolaris VCF (merged petpet and petfal) 


#which ID falls into which window
registerDoParallel(cores=45)
ID_list_window<-foreach(i=1:nrow(vanil_50K), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(ID_list,
                         chr == as.character(vanil_50K[i,2]) &
                         pos >= as.numeric(vanil_50K[i,3]) &
                         pos <= as.numeric(vanil_50K[i,4])),
                ID),
         window=as.character(vanil_50K[i,1]))
}
rm(ID_list)


fread("/data/users/mjahani/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>%
  group_by(chr) %>%
  arrange(pos,.by_group = TRUE) %>%
  mutate(delta_cM=cM-lag(cM, default = first(cM)),
         delta_bp=pos-lag(pos, default = first(pos)),
         start=lag(pos)) %>%
  ungroup() %>%
  mutate(recomb_rate=delta_cM/delta_bp) %>%
  select(chr,start,end=pos,recomb_rate) %>%
  filter(!is.na(start)) -> recombination #recombination rates for 1mb windows


registerDoParallel(cores=40)
window_recombination<-foreach(i=1:nrow(recombination), .combine='rbind', .errorhandling='stop') %dopar% {
  mutate(select(filter(vanil_50K,
                       chr == as.character(recombination[i,1]) &
                       start >= as.numeric(recombination[i,2]) & 
                         end <= as.numeric(recombination[i,3])),
                window),
         recomb_rate=as.character(recombination[i,4]))}
rm(recombination)
left_join(ID_list_window,window_recombination) %>%
  mutate(recomb_rate=as.numeric(recomb_rate)) -> ID_list_window_recombination #add recombination rate to 5kb windows

#add recombination quantile to 5kb windows
ID_list_window_recombination$quantile <- with(ID_list_window_recombination, cut(recomb_rate, 
                                                              breaks=quantile(recomb_rate, probs=seq(0,1, by=0.2), na.rm=TRUE),
                                                              labels=c("0-20","20-40","40-60","60-80","80-100"), 
                                                              include.lowest=TRUE))
  
#add inversion information to the SNP list (petfal)
registerDoParallel(cores=40)
petfal_inversion<-foreach(i=1:nrow(inversions.regions.petfal), .combine='rbind', .errorhandling='stop') %dopar% {
  mutate(select(filter(ID_list,
                       chr == as.character(inversions.regions.petfal[i,4]) &
                         pos >= as.numeric(inversions.regions.petfal[i,2]) & 
                         pos <= as.numeric(inversions.regions.petfal[i,3])),
                ID,chr),
         mds=as.character(inversions.regions.petfal[i,5]))}


#add inversion information to the SNP list (petpet)
registerDoParallel(cores=40)
petpet_inversion<-foreach(i=1:nrow(inversions.regions.petpet), .combine='rbind', .errorhandling='stop') %dopar% {
  mutate(select(filter(ID_list,
                       chr == as.character(inversions.regions.petpet[i,4]) &
                         pos >= as.numeric(inversions.regions.petpet[i,2]) & 
                         pos <= as.numeric(inversions.regions.petpet[i,3])),
                ID,chr),
         mds=as.character(inversions.regions.petpet[i,5]))}



ID_list_window_recombination %>%
  left_join(.,select(mutate(petfal_inversion,petfal_inv=paste0(chr,".",mds)),ID,petfal_inv)) %>%
  left_join(.,select(mutate(petpet_inversion,petpet_inv=paste0(chr,".",mds)),ID,petpet_inv)) %>%
  mutate(petfal_inv=replace_na(petfal_inv,"not_inversion")) %>%
  mutate(petpet_inv=replace_na(petpet_inv,"not_inversion")) %>%
  select(ID,window,recomb_bin=quantile,petfal_inv,petpet_inv)-> ID_list_window_recombination_inversion #data set of ID, window, recombination rate quantile, petfal inversion and petpet inversion


# system("vcftools --gzvcf /moonriseNFS/wild_gwas/Petiolaris.tranche90.snp.gwas.90.bi.remappedHa412HO.vcf.gz --freq --out /data/users/mjahani/FST/Petiolaris")
# 
# fread(paste0("/data/users/mjahani/FST/Petiolaris.frq"),header = F) %>%
#   mutate(ID=paste0(V1,":",V2)) %>%
#   separate(V5,into = c("allele1","freq1"),sep = ":",remove=F) %>%
#   separate(V6,into = c("allele2","freq2"),sep = ":",remove=F) %>%
#   mutate(freq1=as.numeric(freq1),freq2=as.numeric(freq2)) %>%
#   mutate(MAF=ifelse(freq1>=freq2,freq2,freq1)) %>%
#   select(ID,MAF) %>% 
#   right_join(.,ID_list_window_recombination_inversion) %>%
#   filter(MAF > 0.03) %>%
#   select(-MAF) -> ID_list_window_recombination_inversion_frq.3

#data set of ID, window, recombination rate quantile, the inversion name , taxa , inversion count
rbind(mutate(rename(select(ID_list_window_recombination_inversion,-petfal_inv),inversion=petpet_inv),taxa="petpet"),
      mutate(rename(select(ID_list_window_recombination_inversion,-petpet_inv),inversion=petfal_inv),taxa="petfal")
      ) %>% 
  left_join(.,inversion_count) %>%
  mutate(N_inversion=replace_na(N_inversion,0)) %>%
  fwrite("/data/users/mjahani/FST/ID_list_window_recombination_inversion_taxa", sep = "\t",col.names = T)
         
fread("/data/users/mjahani/FST/ID_list_window_recombination_inversion_taxa",header = T) ->  ID_list_window_recombination_inversion_frq.3_taxa
  
sample_inversion %>%
  distinct(inversion) %>%
  pull(inversion) -> INV # list of inversions

sample_inversion %>%
  distinct(taxa) %>%
  pull(taxa) -> TAXA # list of taxa

sample_inversion %>%
  filter(triangle_genotype!=1) %>%
  distinct(triangle_genotype) %>%
  pull(triangle_genotype) -> GENO #list of homozygote genotypes

#list of samples for each homozygote genotype and each allele
  registerDoParallel(cores=40)
foreach(i=1:length(INV), .combine='rbind', .errorhandling='stop') %dopar% { 
  foreach(j=1:length(GENO), .combine='rbind', .errorhandling='stop') %dopar% { 
  sample_inversion %>%
    filter(inversion==INV[i]) %>%
    filter(triangle_genotype==GENO[j]) %>%
    select(sample) %>%
    fwrite(paste0("/data/users/mjahani/FST/temporary_files/",INV[i],"_pop",j),sep = "\t",col.names = F)
  }
}
#list of samples for each homozygote genotype and each allele in a taxa
registerDoParallel(cores=40)
foreach(i=1:length(INV), .combine='rbind', .errorhandling='stop') %dopar% {
  foreach(j=1:length(GENO), .combine='rbind', .errorhandling='stop') %dopar% { 
    foreach(k=1:length(TAXA), .combine='rbind', .errorhandling='stop') %dopar% { 

  sample_inversion %>%
    filter(inversion==INV[i]) %>%
    filter(triangle_genotype==GENO[j]) %>%
    filter(taxa==TAXA[k]) %>%
    select(sample) %>%
    fwrite(paste0("/data/users/mjahani/FST/temporary_files/",INV[i],"_",TAXA[k],"_pop",j),sep = "\t",col.names = F)
    }
  }
}

#calculate FST between homozygote individuals for each inversion
foreach(i=1:length(INV), .combine='rbind', .errorhandling='stop') %dopar% {
  system(paste0("vcftools --gzvcf ",
         VCF,
         " --maf 0.03 --weir-fst-pop /data/users/mjahani/FST/temporary_files/",
         INV[i],"_pop1",
         " --weir-fst-pop /data/users/mjahani/FST/temporary_files/",
         INV[i],"_pop2",
         " --out  /data/users/mjahani/FST/temporary_files/",
         INV[i]))
}
#calculate FST between subspecies in a group of individuals with the same genotype for each inversion (only homozygote)
foreach(i=1:length(INV), .combine='rbind', .errorhandling='stop') %dopar% {
  foreach(j=1:length(GENO), .combine='rbind', .errorhandling='stop') %dopar% { 
  system(paste0("vcftools --gzvcf ",
                VCF,
                " --maf 0.03 --weir-fst-pop /data/users/mjahani/FST/temporary_files/",
                INV[i],"_",TAXA[1],"_pop",j,
                " --weir-fst-pop /data/users/mjahani/FST/temporary_files/",
                INV[i],"_",TAXA[2],"_pop",j,
                " --out  /data/users/mjahani/FST/temporary_files/" ,
                INV[i],"_pop",j))
  }
}

###################


Global_Inversion_FST <- NULL
for (i in 1:length(INV)) {
  fread(paste0("/data/users/mjahani/FST/temporary_files/",
               INV[i],
               ".weir.fst"),header = T) %>%
    mutate(ID=paste0(CHROM,":",POS)) %>%
             select(ID,Fst=WEIR_AND_COCKERHAM_FST) -> FST
           
    ID_list_window_recombination_inversion_frq.3_taxa %>%
      left_join(.,FST) %>%
      filter(!is.na(Fst)) -> window_FST
      #mutate(group=ifelse(inversion==INV[i],INV[i],"not_inversion"))-> window_FST
    rm(FST)
    
    # window_FST %>%
    #   group_by(window,group) %>%
    #   tally() %>%
    #   ungroup() %>%
    #   group_by(window,group) %>%
    #   filter(n==max(n)) %>%
    #   ungroup() %>%
    #   select(-n) %>%
    #   filter(group=="not_inversion") %>%
    #   distinct(window) %>%
    #   pull(window) -> noinversion_window #list of windows outside the inversion-region to FST
    
    window_FST %>%
      filter(inversion=="not_inversion") %>%
      distinct(window) %>%
      pull(window) -> noinversion_window #list of windows outside the inversion-region to FST
    
    
    # window_FST %>%
    #   group_by(window,group) %>%
    #   tally() %>%
    #   ungroup() %>%
    #   group_by(window,group) %>%
    #   filter(n==max(n)) %>%
    #   ungroup() %>%
    #   select(-n) %>%
    #   filter(group==INV[i]) %>%
    #   distinct(window) %>%
    #   pull(window) -> inversion_window  #list of windows within the inversion-region to FST
      
    window_FST %>%
      filter(inversion==as.character(INV[i])) %>%
      distinct(window) %>%
      pull(window) -> inversion_window   #list of windows within the inversion-region to FST
    
    
    window_FST %>%
      filter(window %in% noinversion_window) %>%
      distinct(ID,.keep_all = T) %>%
      group_by(window,recomb_bin) %>%
      summarise(mean_Fst=mean(as.numeric(Fst))) %>%
      ungroup() %>%
      group_by(recomb_bin) %>%
      sample_n(10000,replace = F) %>%
      ungroup() %>%
      select(recomb_bin,mean_Fst)-> null_dist
  
    rm(noinversion_window)
    
    window_FST %>%
      filter(window %in% inversion_window) %>% 
      distinct(recomb_bin) %>%
      pull(recomb_bin) -> RECOM_BIN
    
    for (j in 1:length(RECOM_BIN)) {
    window_FST %>%
      filter(window %in% inversion_window) %>%
      filter(recomb_bin == as.character(RECOM_BIN[j])) %>%  
      distinct(ID,.keep_all = T) %>%
      group_by(window,recomb_bin) %>%
      summarise(global_inversion_Fst=mean(as.numeric(Fst))) %>%
      ungroup() -> window_FST_inversion_recombination
      null_dist %>%
        filter(recomb_bin == as.character(RECOM_BIN[j])) -> null_dist_recom
      
      for (k in 1:nrow(window_FST_inversion_recombination)) {
        window_FST_inversion_recombination[k,4] <-  as.numeric((tally(filter(null_dist_recom,as.numeric(mean_Fst) < as.numeric(window_FST_inversion_recombination[k,3]))))/10000)
      }
        
      window_FST_inversion_recombination %>%
      mutate(inversion=INV[i]) %>%
      select(window,inversion,recomb_bin,global_inversion_Fst,P_global_inversion_Fst=V4) %>%
      rbind(.,Global_Inversion_FST) -> Global_Inversion_FST
    }
    rm(inversion_window,null_dist,window_FST,null_dist_recom)
} 

Species_Inversion_FST <- NULL
for (i in 1:length(INV)) {
  for (j in 1:length(GENO)) {

  fread(paste0("/data/users/mjahani/FST/temporary_files/",
               INV[i],"_pop",j,
               ".weir.fst"),header = T) %>%
    mutate(ID=paste0(CHROM,":",POS)) %>%
    mutate(Fst=ifelse(WEIR_AND_COCKERHAM_FST=="NaN",NA,WEIR_AND_COCKERHAM_FST)) %>%
    select(ID,Fst)  -> FST
  
  ID_list_window_recombination_inversion_frq.3_taxa %>%
    left_join(.,FST) %>%
    filter(!is.na(Fst)) -> window_FST
    #mutate(group=ifelse(inversion==INV[i],INV[i],"not_inversion"))-> window_FST
  rm(FST)
  # window_FST %>%
  #   group_by(window,group) %>%
  #   tally() %>%
  #   ungroup() %>%
  #   group_by(window,group) %>%
  #   filter(n==max(n)) %>%
  #   ungroup() %>%
  #   select(-n) %>%
  #   filter(group=="not_inversion") %>%
  #   distinct(window) %>%
  #   pull(window) -> noinversion_window
  window_FST %>%
    filter(inversion=="not_inversion") %>%
    distinct(window) %>%
    pull(window) -> noinversion_window

  # window_FST %>%
  #   group_by(window,group) %>%
  #   tally() %>%
  #   ungroup() %>%
  #   group_by(window,group) %>%
  #   filter(n==max(n)) %>%
  #   ungroup() %>%
  #   select(-n) %>%
  #   filter(group==INV[i]) %>%
  #   distinct(window) %>%
  #   pull(window) -> inversion_window
  
  window_FST %>%
    filter(inversion==as.character(INV[i])) %>%
    distinct(window) %>%
    pull(window) -> inversion_window 
  
  window_FST %>%
    filter(window %in% noinversion_window) %>%
    distinct(ID,.keep_all = T) %>%
    group_by(window,recomb_bin) %>%
    summarise(mean_Fst=mean(as.numeric(Fst))) %>%
    ungroup() %>%
    group_by(recomb_bin) %>%
    sample_n(10000,replace = F) %>%
    ungroup() %>%
    select(recomb_bin,mean_Fst)-> null_dist
  
  window_FST %>%
    filter(window %in% inversion_window) %>% 
    distinct(recomb_bin) %>%
    pull(recomb_bin) -> RECOM_BIN

  for (l in 1:length(RECOM_BIN)) {
    window_FST %>%
      filter(window %in% inversion_window) %>%
      filter(recomb_bin == as.character(RECOM_BIN[l])) %>%  
      distinct(ID,.keep_all = T) %>%
      group_by(window,recomb_bin) %>%
      summarise(species_inversion_Fst=mean(as.numeric(Fst))) %>%
      ungroup() -> window_FST_inversion_recombination
    null_dist %>%
      filter(recomb_bin == as.character(RECOM_BIN[l])) -> null_dist_recom
    
    for (k in 1:nrow(window_FST_inversion_recombination)) {
      window_FST_inversion_recombination[k,4] <-  as.numeric((tally(filter(null_dist_recom,as.numeric(mean_Fst) < as.numeric(window_FST_inversion_recombination[k,3]))))/10000)
    }
    
    window_FST_inversion_recombination %>%
      mutate(inversion=INV[i]) %>%
      mutate(genotype=GENO[j]) %>%
      select(window,inversion,recomb_bin,genotype,species_inversion_Fst,P_species_inversion_Fst=V4) %>%
      rbind(.,Species_Inversion_FST) -> Species_Inversion_FST
  }
  rm(inversion_window,null_dist,window_FST,null_dist_recom)
}}

sample_inversion %>% 
  filter(triangle_genotype!=1) %>% 
  group_by(inversion,taxa,triangle_genotype) %>% 
  tally() %>% 
  ungroup() %>%
  rename(allele=triangle_genotype,n_individuals=n) %>%
  arrange(inversion,allele) -> a
a %>% 
  group_by(inversion,allele) %>%
  summarise(total=sum(as.numeric(n_individuals))) %>%
  ungroup()-> b
left_join(a,b) %>%
  mutate(frequency=as.numeric(n_individuals)/as.numeric(total)) %>%
  select(inversion,taxa,allele,n_individuals,frequency) %>%
  fwrite(paste0("/data/users/mjahani/FST/inversion_numberofsamples"),sep = "\t",col.names = T,na = "NA")

fread("/data/users/mjahani/Null_W_Result/files/convergent_2binning",sep = "\t",header=T) %>%
  filter(comparison=="petfal_petpet") %>%
  distinct(window5) %>%
  rename(window=window5) %>%
  mutate(Convergent="TRUE")-> Convergent_window

fread("/data/users/mjahani/FST/inversion_numberofsamples",sep = "\t",header=T) %>% 
  filter(n_individuals < 5) %>%
  filter(allele == 0) %>%
  distinct(inversion) %>%
  pull(inversion) -> low_freq_all0
  
fread("/data/users/mjahani/FST/inversion_numberofsamples",sep = "\t",header=T) %>% 
  filter(n_individuals < 5) %>%
  filter(allele == 2) %>%
  distinct(inversion) %>%
  pull(inversion) -> low_freq_all2

Species_Inversion_FST %>% 
  filter(!is.na(window)) %>%
  filter(genotype==0) %>%
  select(-genotype) %>%
  rename(allele0_species_inversion_Fst=species_inversion_Fst,
         allele0_species_inversion_Fst_P=P_species_inversion_Fst) %>%
  filter(!inversion %in% low_freq_all0) %>%
  full_join(.,
            Species_Inversion_FST %>% 
              filter(!is.na(window)) %>%
              filter(genotype==2) %>%
              select(-genotype) %>%
              rename(allele2_species_inversion_Fst=species_inversion_Fst,
                     allele2_species_inversion_Fst_P=P_species_inversion_Fst) %>%
              filter(!inversion %in% low_freq_all2)) %>%
  full_join(.,
            Global_Inversion_FST %>%
              select(window,inversion,recomb_bin,global_inversion_Fst,global_inversion_Fst_P=P_global_inversion_Fst)) %>%
  left_join(.,
            ID_list_window_recombination_inversion_frq.3_taxa %>%
              distinct(window,inversion,N_inversion)) %>%
  left_join(.,
            Convergent_window) %>%
  mutate(is_convergent=replace_na(Convergent,"FALSE")) %>%
  mutate(is_inversion_present_in_both_species= ifelse(N_inversion==2,"TRUE","FALSE")) %>%
  select(window,
        inversion,
        is_inversion_present_in_both_species,
        is_window_convergent=is_convergent,
        global_inversion_Fst,
        global_inversion_Fst_P,
        allele0_species_inversion_Fst,
        allele0_species_inversion_Fst_P,
        allele2_species_inversion_Fst,
        allele2_species_inversion_Fst_P) %>% 
  fwrite(paste0("/data/users/mjahani/FST/inversion_window_Fst_summary_new"),sep = "\t",col.names = T,na = "NA")


#to correct is.convergent column based on corrected_GWAS and sprearman
fread("/data/users/mjahani/Null_W_Result/files/convergent_2binning",sep = "\t",header=T) %>% 
  filter(comparison=="petfal_petpet") %>%
  filter(analysis=="spearman" | analysis=="GWAS_corrected")  %>% 
  distinct(window5) %>%
  rename(window=window5) %>%
  mutate(Convergent="TRUE")-> Convergent_window

fread("/data/users/mjahani/FST/inversion_window_Fst_summary_new",sep = "\t",header=T) %>%
  select(-is_window_convergent) %>%
  left_join(.,Convergent_window) %>%
  mutate(is_window_convergent=replace_na(Convergent,"FALSE")) %>%
  select(window,
         inversion,
         is_inversion_present_in_both_species,
         is_window_convergent,
         global_inversion_Fst,
         global_inversion_Fst_P,
         allele0_species_inversion_Fst,
         allele0_species_inversion_Fst_P,
         allele2_species_inversion_Fst,
         allele2_species_inversion_Fst_P) %>% 
  fwrite(paste0("/data/users/mjahani/FST/inversion_window_Fst_summary_new_GWAcor_spearman"),sep = "\t",col.names = T,na = "NA")
  
  
########################

fread("/data/users/mjahani/FST/inversion_window_Fst_summary",sep = "\t",header=T) %>% nrow()
  mutate(test=ifelse(allele2_species_inversion_Fst<0,"negative","postive")) %>%
  group_by(inversion,test) %>%
  tally() %>%
  ungroup() %>%
  spread(test,n) %>%
  mutate(aa=ifelse(as.numeric(negative)>as.numeric(postive),"yes","no"),
         bb=as.numeric(negative)-as.numeric(postive)) %>%
  arrange(desc(aa),desc(bb))->a

global_inversion_Fst  
Ha412HOChr06.syn2 
Ha412HOChr17.syn4

allele0_species_inversion_Fst
Ha412HOChr01.syn1
Ha412HOChr13.syn1

allele2_species_inversion_Fst
Ha412HOChr07.syn1
Ha412HOChr08.syn1

list.files(path = "/data/users/mjahani/FST/temporary_files") -> file_list
data <- NULL
for (i in 1:length(file_list)) {
  fread(paste0("/data/users/mjahani/FST/temporary_files/",file_list[i]),sep="\t",header = F) %>%
    tally() %>%
    mutate(name=file_list[i]) %>%
    rbind(.,data) -> data }

data  %>% 
  filter(!grepl(".fst$", name)) %>%
  filter(!grepl(".log$", name)) %>%
  filter(grepl("*pet*", name)) %>%
  separate(name,into = c("inversion","taxa","allele"),sep="_",remove = T) %>%
  mutate(inv_alele=paste0(inversion,"_",allele)) %>%
  select(-inversion,-allele) %>% 
  spread(taxa,n) %>% 
  mutate(diff=abs(as.numeric(petfal)-as.numeric(petpet))) %>%
  filter(as.numeric(petfal)!=0) %>%
  filter(as.numeric(petpet)!=0) %>%
  rename(inversion=inv_alele,pop1=petfal,pop2=petpet)-> a

data  %>% 
  filter(!grepl(".fst$", name)) %>%
  filter(!grepl(".log$", name)) %>%
  filter(!grepl("*pet*", name)) %>%
  separate(name,into = c("inversion","allele"),sep="_",remove = T) %>%
  spread(allele,n) %>% 
  mutate(diff=abs(as.numeric(pop2)-as.numeric(pop1))) %>%
  filter(as.numeric(pop1)!=0) %>%
  filter(as.numeric(pop2)!=0)-> b

rbind(a,b) %>% 
  mutate(smaller_n=ifelse(as.numeric(pop1)>as.numeric(pop2),as.numeric(pop2),as.numeric(pop1))) %>%
           select(inversion,smaller_n,diff)-> n_diff

data.frame(file_list) %>% filter(grepl(".fst$", file_list)) -> file_list
data1 <- NULL
for (i in 1:nrow(file_list)) {
  fread(paste0("/data/users/mjahani/FST/temporary_files/",file_list[i,1]),sep="\t",header = T) %>%
    mutate(direction=ifelse(as.numeric(WEIR_AND_COCKERHAM_FST)<0,"negative","postive")) %>%
    select(direction) %>%
    rbind(.,data.frame(direction=c("negative","postive"))) %>%
    group_by(direction) %>%
    tally() %>%
    ungroup() %>%
    mutate(name=file_list[i,1]) %>%
    spread(direction,n) %>%
    mutate(total_n=(as.numeric(negative)-1)+(as.numeric(postive)-1)) %>%
    mutate(neg_ratio=(as.numeric(negative)-1)/as.numeric(total_n)) %>%
    select(name,neg_ratio,total_n) %>%
    rbind(.,data1) -> data1 }

data1  %>% 
  filter(total_n!=0) %>%
  mutate(inversion=gsub(".weir.fst","",name)) %>%
  select(inversion,neg_ratio,total_n) %>%
  full_join(.,n_diff) %>%
  fwrite(paste0("/data/users/mjahani/FST/inversion_negativeratio_popsize"),sep = "\t",col.names = T)
  
  fread("/Users/mojtabajahani/Downloads/inversion_negativeratio_popsize",sep = "\t",header = T) -> inversion_negativeratio_popsize

  ggplot(inversion_negativeratio_popsize,aes(x=smaller_n,y=neg_ratio*100)) + 
    geom_point() + 
    geom_smooth(method = lm) +
    theme_bw() +
    theme( 
      panel.border = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.y.right = element_text(color = "blue")
    )+
    xlab("Smaller Population Size")+
    ylab("Negative Loci Percentage")
   
  inversions.regions.petpet %>%
    mutate(size=floor(((as.numeric(end)-as.numeric(start))+1)/5000),
           inversion=paste0(chr,"_",mds)) %>% 
    select(inversion,size) %>%
    group_by(inversion) %>%
    summarise(N_win=sum(as.numeric(size)))-> a
  min(a$N_win)
  max(a$N_win)
    
  fread("/data/users/mjahani/Null_W_Result/files/convergent_2binning",sep = "\t",header=T) %>%
    filter(comparison=="petfal_petpet") %>%
    distinct(analysis,window5) %>%
    fwrite(paste0("/data/users/mjahani/converget&correction"),sep = "\t",col.names = T,na = "NA")
  
  