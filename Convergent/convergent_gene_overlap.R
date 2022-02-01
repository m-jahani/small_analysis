library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)

gff3="/home/mjahani/Ha412HOv2.0-20181130_genes_v2.1.edited.gff3"
range=500
bed<-read.table(
  text = 
    gsub(";|=", "\t", # Sunbtitude different sepratots with tab
         readLines(gff3)),
  header =FALSE, 
  sep ='\t',
  fill=TRUE, 
  quote = "", 
  row.names = NULL, 
  stringsAsFactors = FALSE) %>% 
  filter(V3=="gene") %>% 
  filter(V1 %in% pull(mutate(data.frame(N=seq(1,17)),chr=ifelse(N<10,paste0("Ha412HOChr0",N),paste0("Ha412HOChr",N))),chr)) %>%
  select(chrom=V1,chromStart=V4,chromEnd=V5,name=V10) %>%
  mutate(chromStart=as.numeric(chromStart)-range,chromEnd=as.numeric(chromEnd)+range)



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
  select(window,chr,start,end) -> vanil_50K

registerDoParallel(cores=40)
window_gene<-foreach(i=1:nrow(bed), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(vanil_50K,
                       chr == as.character(bed[i,1]) &
                         ((start >= as.numeric(bed[i,2]) & start <= as.numeric(bed[i,3])) | 
                            (end >= as.numeric(bed[i,2]) & end <= as.numeric(bed[i,3])) | 
                            (start < as.numeric(bed[i,2]) & end > as.numeric(bed[i,3])))),
                window),
         Gene=as.character(bed[i,4]))}

 select(vanil_50K,window) %>% left_join(.,window_gene) -> vanil_50K_Gene



fread("/data/users/mjahani/Null_W_Result/files/convergent_2binning",sep = "\t",header=T) %>%
  left_join(.,rename(vanil_50K_Gene,window5=window)) -> convergent_window_gene
 
convergent_window_gene %>%
  fwrite("/data/users/mjahani/convergent_2binning_gene", sep = "\t",col.names = T,na = "NA") 

convergent_window_gene %>%
  group_by(data,type,analysis,comparison,variable) %>%
  summarise(N_convergent_window=n_distinct(window5)) %>%
  ungroup() -> a

convergent_window_gene %>%
  filter(!is.na(Gene)) %>%
  group_by(data,type,analysis,comparison,variable) %>%
  summarise(N_convergent_window_gene_overlap=n_distinct(window5)) %>%
  ungroup() -> b
  
full_join(a,b) %>%
  mutate(N_convergent_window_gene_overlap=replace_na(N_convergent_window_gene_overlap,0)) %>%
  mutate(overlap=N_convergent_window_gene_overlap/N_convergent_window) -> gene_overlap_summary
rm(a,b)

#distinct(vanil_50K_Gene,window,.keep_all = T) -> vanil_50K_Gene_distinct

vanil_50K_Gene %>%
  mutate(Gene_1=ifelse(is.na(Gene),0,1)) %>%
  group_by(window) %>%
  summarise(N_gene=sum(Gene_1)) -> vanil_50K_Gene_distinct

for (i in 1:nrow(gene_overlap_summary)) {
    registerDoParallel(cores=40)
    null_distrbution <-foreach(j=1:10000, .combine='rbind', .errorhandling='stop') %dopar% {
    vanil_50K_Gene_distinct %>%
      sample_n(as.numeric(gene_overlap_summary[i,6]),replace = F) %>%
      filter(as.numeric(N_gene) > 0) %>%
      tally() %>%
      mutate(round=j)}
      
    null_distrbution %>%
      filter(n < as.numeric(gene_overlap_summary[i,7])) %>%
      nrow() -> gene_overlap_summary[i,9]
    
    null_distrbution %>%
      filter(n > as.numeric(gene_overlap_summary[i,7])) %>%
      nrow() -> gene_overlap_summary[i,10]
    
    mean(as.numeric(null_distrbution$n)/as.numeric(gene_overlap_summary[i,6])) -> gene_overlap_summary[i,11]
    rm(null_distrbution)  
  }
  
gene_overlap_summary %>%
  mutate(P_value=(as.numeric(V9)/10000)) %>%
  mutate(Null_dist_mean_overlap=as.numeric(V11)/as.numeric(N_convergent_window)) %>%
  select(data,type,analysis,comparison,variable,
         N_convergent_window,N_convergent_window_gene_overlap,
         overlap,Null_dist_mean_overlap=V11,P_value) %>%
  fwrite(paste0("/data/users/mjahani/convergent_gene_overlap/convergent_gene_over"),sep = "\t",col.names = T,na = "NA")

fread("/Users/mojtabajahani/Downloads/convergent_gene_over",sep = "\t",header = T) %>%
  -> a

