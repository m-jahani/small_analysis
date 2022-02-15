library(tidyverse)
library(data.table)


#wget https://raw.githubusercontent.com/owensgl/wild_gwas_2018/master/MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt
#wget https://raw.githubusercontent.com/owensgl/wild_gwas_2018/master/MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt

fread("/data/users/mjahani/Ha412HO_inv.v3.inversions.regions.v1.txt") -> inversions.regions
fread("/data/users/mjahani/Ha412HO_inv.v3.petinvfreq.txt") -> Petiolaris.frequency

inversions.regions %>%
  filter(species!="Niveus",species!="Petiolaris") -> inversions.regions.ann.arg

inversions.regions %>%
  filter(species=="Petiolaris") %>%
  select(-species) %>%
  mutate(species="PetFal") %>%
  left_join(.,Petiolaris.frequency)%>%
  filter(as.numeric(freq) > 0) %>%
  select(-species) %>%
  mutate(species="petfal") %>%
  select(group,start,end,chr,spe,species,dataset,mds) -> inversions.regions.petfal

inversions.regions %>%
  filter(species=="Petiolaris") %>%
  select(-species) %>%
  mutate(species="PetPet") %>%
  left_join(.,Petiolaris.frequency)%>%
  filter(as.numeric(freq) > 0) %>%
  select(-species) %>%
  mutate(species="petpet") %>%
  select(group,start,end,chr,spe,species,dataset,mds) -> inversions.regions.petpet


rbind(inversions.regions.petfal,inversions.regions.petpet) -> inversions.regions.Petiolaris
rbind(inversions.regions.Petiolaris,inversions.regions.ann.arg) -> inversions.regions.all


fread("/data/users/mjahani/Ha412HO_inv.v3.txt") %>%
  mutate(chr=ifelse(chr>=10,paste0("Ha412HOChr",chr),paste0("Ha412HOChr0",chr))) %>% mutate(mds=paste0(direction,mds)) %>%
  filter(species!="Niveus") %>% select(-species) %>% left_join(inversions.regions.all,.) %>% select(chr,start,end,species,sv_name) %>%
  fwrite("/data/users/mjahani/inversion.bed",sep = "\t",col.names= T)