library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly = TRUE)



GWA_result1 <- args[1] #a file with two columns,first column: sample name, second column: popname. No header
GWA_result2 <- args[2] #index file in the same directory
clump_result <- args[3] #directory to save the analysis. ends with /
SAVE_DIR <-  args[4] #directory to vcf2baypass.pl


fread(GWA_result1) %>% 
  select(SNP=V1,
         beta=V2,
         P=V4) %>% 
  filter((-1*log10(P))>4) %>% #7.310691
  rbind(.,
        fread(GWA_result2) %>% 
          select(SNP=V1,
                 beta=V2,
                 P=V4) %>% 
          filter((-1*log10(P))>4)) %>% #7.310691
  group_by(SNP) %>%
  filter(P == min(P)) %>%
  ungroup() %>%
  separate(SNP,
           into = c("CHR","POS"),
           sep = ":",
           remove = F) -> SNPS



fread(clump_result) %>% 
  select(SNP,
         POS,
         KB) %>%
  filter(KB > 0.001) %>%
  separate(SNP,
           into = c("CHR","POSITION"),
           sep = ":",
           remove = T) %>%
  separate(POS,
           into = c("CHR1","START_END"),
           sep = ":",
           remove = T) %>%
  separate(START_END,
           into = c("START","END"),
           remove = T) %>%
  mutate(clump_ID=paste0(CHR,":",START,"-",END)) %>%
  select(CHR,
         START,
         END,
         clump_ID)  -> CLUMPS


registerDoParallel(cores=48)
SNP_clump <- foreach(i=1:nrow(CLUMPS), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(SNPS,
                       CHR == as.character(CLUMPS[i,1]),
                       (POS >= as.numeric(CLUMPS[i,2]) &
                          POS <= as.numeric(CLUMPS[i,3]))
  ),
  CHR,POS),
  clump_ID=as.character(CLUMPS[i,4]))}

SNPS %>%
  left_join(.,SNP_clump) %>%
  filter(!is.na(clump_ID)) %>%
  select(SNP,
         beta,
         clump_ID) -> SNP_B_CLUMP


fread("/data/users/mjahani/febina/genotype_covar/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.tfam") %>%
  select(V1) %>%
  pull(V1) -> SAMPLES

SNP_B_CLUMP %>% 
  select(SNP) %>%
  fwrite(paste0(SAVE_DIR,"/temp_snp_list.txt"),
         col.names = F,
         quote = F,
         sep = "\t")

  registerDoParallel(cores=25)
foreach(i=1:length(SAMPLES), .combine='rbind', .errorhandling='stop') %dopar% { 
  
  data.frame(V1=SAMPLES[i],
             V2=SAMPLES[i]) %>%
    fwrite(paste0(SAVE_DIR,"/temp_sample_list_",
                  SAMPLES[i]),
           col.names = F,
           quote = F,
           sep = "\t")
  
  system(paste0("/home/mjahani/bin/plink --tfile /data/users/mjahani/febina/genotype_covar/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99 --keep ",
                SAVE_DIR,
                "/temp_sample_list_",
                SAMPLES[i],
                " --extract ",
                SAVE_DIR,
                "/temp_snp_list.txt --freq --out ",
                SAVE_DIR,
                "/temp_frequancy_",
                SAMPLES[i]),
         ignore.stdout = T)
  
  fread(paste0(SAVE_DIR,
               "/temp_frequancy_",
               SAMPLES[i],
               ".frq"
               ),
        header = T) %>%
    select(SNP,
           MAF) %>%
    left_join(.,
              SNP_B_CLUMP) %>%
    mutate(effect=(1-MAF)*beta) %>%
    group_by(clump_ID) %>%
    summarize(sum_effect=sum(effect)) %>%
    ungroup() %>%
    mutate(sample=SAMPLES[i]) %>%
    fwrite(paste0(SAVE_DIR,
                  "/",
                  gsub(".ps","",gsub(".*/","",GWA_result1)),
                  "_",
                  gsub(".ps","",gsub(".*/","",GWA_result2)),
                  "_clump_effect"),
           sep = "\t",
           col.names = F,
           quote = F,
           append = T,
           na = "NA")

    
  system(paste0("rm ",
                SAVE_DIR,
                "/temp_sample_list_",
                SAMPLES[i],
                " ",
                SAVE_DIR,
                "/temp_frequancy_",
                SAMPLES[i],
                "*"))
  
}

system(paste0("rm ",
       SAVE_DIR,
       "/temp_snp_list.txt"))
