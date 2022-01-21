library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(gtools)


tmp_add <- "/data/users/mjahani/baypass/"
rslt_add <- "/data/users/mjahani/Convergent/Baypass/"
Taxa_var<-read.table(paste0(tmp_add,"uniq_Taxa_var"),
                     header =F, 
                     sep ='\t',
                     row.names = NULL,
                     stringsAsFactors = FALSE) %>% filter(V1!="Taxa_var")  %>% 
  separate(V1,into=c("Taxa","Variable"),sep=":",extra="warn",remove = F) %>%
  filter(Variable %in% (read.table(paste0(tmp_add,"uniq_Taxa_var"),
                                   header =F, 
                                   sep ='\t',
                                   row.names = NULL,
                                   stringsAsFactors = FALSE) %>% filter(V1!="Taxa_var")  %>% 
                          separate(V1,into=c("Taxa","Variable"),sep=":",extra="warn",remove = T) %>%
                          group_by(Variable) %>% tally() %>% ungroup() %>% filter(n >= 2) %>% pull(Variable))) %>% pull(V1)

registerDoParallel(cores=25)

percentile <- foreach(i=1:length(Taxa_var), .combine='rbind', .errorhandling='stop') %dopar% {
  data.frame(Taxa_var=Taxa_var[i],Perc=10)
}


fread("/moonriseNFS/Null_W_Result/files/5kbwindow_recombination",sep="\t")-> window_recomb

registerDoParallel(cores=25)
SNP_count <- foreach(i=1:nrow(percentile), .combine='rbind', .errorhandling='stop') %dopar% {
  fread(paste0(tmp_add,percentile[i,1]),showProgress = F) %>% 
    select(Taxa_var=V6,Value=V4,window5=V5) %>%
    mutate(outlier=ifelse(Value >= percentile[i,2] , "outlier" , "non_outlier")) %>%
    group_by(Taxa_var,window5,outlier) %>% 
    summarise(snp_count = n()) %>%
    ungroup() %>%
    spread(outlier,snp_count) %>%
    mutate(outlier = replace_na(outlier, 0),non_outlier = replace_na(non_outlier, 0)) %>%
    mutate(snp_count =as.numeric(outlier)+as.numeric(non_outlier)) %>%
    select(-non_outlier) %>%
    rename(outlier_count=outlier) %>%
    left_join(.,window_recomb) %>%
    filter(!is.na(recomb_rate)) %>%
    mutate(quantile=cut(recomb_rate, 
                        breaks=quantile(recomb_rate, probs=seq(0,1, by=0.2), na.rm=TRUE),
                        labels=c("0-20","20-40","40-60","60-80","80-100"), 
                        include.lowest=TRUE)) %>%
    select(-recomb_rate)
}

snp_count2<- left_join(SNP_count,
                       SNP_count %>% 
                         group_by(Taxa_var,quantile) %>% 
                         summarise(totsnp=sum(snp_count),totout=sum(outlier_count)) %>% 
                         mutate(expected =totout/totsnp) %>% 
                         select(-totsnp,-totout)) %>% 
  mutate(P4=qbinom (0.9999,snp_count,expected)) %>%
  mutate(Candidates= ifelse(outlier_count > P4,"top","not"))

snp_count2 %>%
  select(Taxa_var,window5,outlier_count,snp_count,expected) %>%
  fwrite(paste0(rslt_add,"files/combination/Taxa_var_count_cutoff"),sep = "\t")
  

snp_count2 %>%
  select(Taxa_var,window5,Candidates) ->  snp_count2



fread("/moonriseNFS/mojtaba_tmp/LD/all_id_gene_window_number_maf.03",sep = "\t",header = T) %>%
  distinct(window5,window)->convert
convert2<-data.frame(species= c("H.annuus","H.argophyllus","H.petiolaris.fallax","H.petiolaris.petiolaris"),
                     Taxa= c("Annuus","Argophyllus","petfal","petpet"))
fread("/moonriseNFS/mojtaba_tmp/LD/all_id_gene_window_number_maf.03",sep = "\t",header = T) %>% distinct(window5,window)->convert

left_join(SNP_count,SNP_count %>%
            group_by(Taxa_var,quantile) %>%
            summarise(totsnp=sum(snp_count),totout=sum(outlier_count)) %>%
            mutate(expected =totout/totsnp) %>%
            select(-totsnp,-totout)) %>%
  mutate(P4=qbinom (0.9999,snp_count,expected)) %>%
  mutate(P8=qbinom (0.99999999,snp_count,expected)) %>%
  separate(Taxa_var,into=c("Taxa","Variable"),sep=":",extra="warn",remove = T) %>%
  left_join(.,convert2) %>%
  left_join(.,convert) %>%
  select(test_name=Variable,window_5K=window,snp_count,outlier_count,expected,p4=P4,p8=P8,species) %>%
  fwrite(paste0(rslt_add,"super_outlier/Superoutlier"),sep = "\t")




Taxas <- read.table(paste0(tmp_add,"uniq_Taxa_var"),
                    header =F, 
                    sep ='\t',
                    row.names = NULL,
                    stringsAsFactors = FALSE) %>% filter(V1!="Taxa_var")  %>% 
  separate(V1,into=c("Taxa","Variable"),sep=":",extra="warn",remove = F) %>%
  distinct(Taxa) %>% pull()

registerDoParallel(cores=25)
unique_window_taxa <- foreach(i= 1:length(Taxas), .combine='rbind' , .errorhandling='stop') %dopar% {
  snp_count2 %>%
    separate(Taxa_var,into=c("Taxa","Variable"),sep=":",extra="warn",remove = F) %>% 
    filter(Taxa==Taxas[i]) %>% 
    distinct(window5,.keep_all = T) %>%
    select(Taxa,window5) 
}

combination_Taxa <- permutations(n=4,r=2,v=Taxas,repeats.allowed=F)

orth_window_taxa <- foreach(i= 1:nrow(combination_Taxa), .combine='rbind') %dopar% {
  unique_window_taxa %>%
    filter(Taxa==combination_Taxa[i,1] | Taxa==combination_Taxa[i,2]) %>%
    group_by(window5) %>%
    tally() %>%
    filter(n == 2) %>%
    select(window5) %>%
    mutate(Taxa = combination_Taxa[i,1]) %>%
    mutate(Taxa2 = combination_Taxa[i,2]) 
}


top_compare<-snp_count2 %>%
  separate(Taxa_var,into=c("Taxa2","Variable"),sep=":",extra="warn",remove = T) %>% rename(Taxa2_candidate=Candidates)


registerDoParallel(cores=20)
foreach(i=1:nrow(percentile), .combine='rbind', .errorhandling='stop') %dopar% {
  fread(paste0(tmp_add,percentile[i,1]),showProgress = F) %>%
    select(Taxa=V1,Variable=V2,ID=V3,Value=V4,window5=V5,Taxa_var=V6) %>%
    left_join(.,SNP_count) %>%
    filter(!is.na(snp_count)) %>%
    select(-outlier_count,-quantile) %>%
    left_join(.,snp_count2) %>%
    left_join(.,orth_window_taxa) %>%
    filter(!is.na(Taxa2)) %>%
    left_join(.,top_compare) %>%
    filter(!is.na(Taxa2_candidate)) %>%
    fwrite(file = paste0(rslt_add,"files/alltaxa_annotate_outlier.ortho.uncorrected"), append = T,sep = "\t")}


rm(Taxa_var,percentile,SNP_count,snp_count2,Taxas,unique_window_taxa,top_compare)
#system(paste0("rm ",tmp_add,"alltaxa_annotate_freq.noinv"))

registerDoParallel(cores=25)
foreach(i=1:nrow(combination_Taxa), .combine='rbind', .errorhandling='stop') %dopar% {
  system(paste(paste("awk '{if(($1 == ",combination_Taxa[i,1],") && ($9 == ",combination_Taxa[i,2],")) {print $0}}'",sep='"'),paste0(rslt_add,"files/alltaxa_annotate_outlier.ortho.uncorrected > ",rslt_add,"files/combination/",combination_Taxa[i,1],"_",combination_Taxa[i,2])))
}

##############
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(gtools)  


#copy this directory in server and change the addredd below
"/data/users/mjahani/Convergent/Baypass/files/combination/"-> tmp_add
#copy these two files to server and change them below
"/moonriseNFS/Null_W_Result/files/ID_window_final" -> a
"/moonriseNFS/Null_W_Result/files/5kbwindow_recombination" -> b
#directory to save the results
"/moonriseNFS/Null_W_Result/binning/"-> result_directory
#######
Taxas <-c("Annuus","Argophyllus","petfal","petpet")

combination_Taxa <- permutations(n=4,r=2,v=Taxas,repeats.allowed=F)

Annotation<-fread(a,
                  header =F)  %>% select(window_name=V3,Window=V2) %>% distinct(Window, .keep_all = T) 

fread(b,sep="\t")-> window_recomb

for (i in 1:nrow(combination_Taxa)) {#loop through Taxas
  window_recom_ortho<-fread(paste0(tmp_add,combination_Taxa[i,1],"_",combination_Taxa[i,2]),showProgress = F) %>%
    select(window5=V5) %>%
    distinct(window5) %>%
    left_join(.,window_recomb) %>%
    filter(!is.na(recomb_rate)) 
  
  window_recom_ortho$quantile <- with(window_recom_ortho, cut(recomb_rate, 
                                                              breaks=quantile(recomb_rate, probs=seq(0,1, by=0.2), na.rm=TRUE),
                                                              labels=c("0-20","20-40","40-60","60-80","80-100"), 
                                                              include.lowest=TRUE))
  
  two_variable_orth<-fread(paste0(tmp_add,combination_Taxa[i,1],"_",combination_Taxa[i,2]),showProgress = F) %>% 
    select(Variable=V2,Value=V4,ID=V3,window5=V5,snp_count=V7,Candidates=V8,Taxa2_candidate=V10) %>% 
    left_join(.,window_recom_ortho) %>%
    filter(!is.na(recomb_rate)) %>% 
    select(-recomb_rate)
  
  rm(window_recom_ortho)
  
  #list of variables in each comparison
  variable_list  <- two_variable_orth %>%
    distinct(Variable) %>%
    pull(Variable)
  
  quantile_list<- two_variable_orth %>%
    distinct(quantile) %>%
    pull(quantile)
  
  for (z in 1:length(quantile_list)) {
    #list of windows which are top candidate atleast in one variable in each quantile
    top_list <- two_variable_orth %>% 
      filter(quantile == quantile_list[z]) %>%
      filter(Candidates=="top") %>% 
      distinct(window5) %>%
      pull(window5)
    #10000 random SNP_ID in windows that are not top windows in any variable in each quantile
    background_10_ID <- two_variable_orth %>%
      filter(quantile == quantile_list[z]) %>%
      filter(!window5 %in% top_list) %>%
      distinct(ID) %>%
      sample_n(10000) %>%
      pull(ID)
    
    for (j in 1:length(variable_list)) {#loop through variables
      background_10KB <- two_variable_orth %>%
        filter(Variable == variable_list[j]) %>%
        filter(quantile == quantile_list[z]) %>%
        filter(ID %in% background_10_ID) %>%
        pull(Value)
      #list windows for first variable that are not top and atleast have 5 snps
      window_list <- two_variable_orth %>%
        filter(Variable == variable_list[j]) %>%
        filter(quantile == quantile_list[z]) %>%
        filter(!window5 %in% top_list) %>%
        filter(!is.na(Value)) %>%
        filter(snp_count > 5) %>%
        distinct(window5) %>%
        pull(window5)
      
      Values<-two_variable_orth %>%
        filter(Variable==variable_list[j]) %>%
        filter(quantile == quantile_list[z]) %>%
        select(window5,Value)
      #list windows for first variable that are top for second species and atleast have 3 snps
      window_list1 <- two_variable_orth %>%
        filter(Variable == variable_list[j]) %>%
        filter(quantile == quantile_list[z]) %>%
        filter(Taxa2_candidate=="top") %>%
        filter(snp_count > 3) %>%
        distinct(window5) %>%
        pull(window5)
      # number of uinique windows in first variable which are are top for second species
      N_window <- two_variable_orth %>%
        filter(Variable == variable_list[j]) %>%
        filter(quantile == quantile_list[z]) %>%
        filter(Taxa2_candidate=="top") %>%
        distinct(window5) %>%
        nrow()
     
      Values1<-two_variable_orth %>%
        filter(quantile == quantile_list[z]) %>%
        filter(Variable==variable_list[j]) %>%
        select(window5,Value)
      
      if (length(window_list1)>0) {
        registerDoParallel(cores=40)
        nulldist<-foreach(k=1:length(window_list), .combine='rbind', .errorhandling='stop') %dopar% {#null distribution
          
          data.frame(N= sum(is.na(pull(filter(Values,window5==window_list[k])))==F),#how many snp for each not top window with  atleast  5 snps
                     W=as.vector(wilcox.test(background_10KB, pull(filter(Values,window5==window_list[k])))$statistic)) %>% #wilcox test, x=(10000 p_values in first variable that are not in top windows for first species) y=(P_values of first non_top window in first variable)
            filter(!is.na(W)) %>%                                                                                           
            mutate(Zscore= (2*W-10000*N) / sqrt(10000*N*(10000+N+1)/3)) %>%
            select(Zscore)
        }#null distribution
        rm(Values,window_list)
        
        Result<-foreach(l=1:length(window_list1), .combine='rbind', .errorhandling='stop') %dopar% {
          
          data.frame(Taxa1=combination_Taxa[i,1],
                     Taxa2=combination_Taxa[i,2],
                     Variable=variable_list[j],
                     quantile=as.character(quantile_list[z]),
                     Window=window_list1[l],
                     W=as.vector(wilcox.test (pull(filter(Values1,window5==window_list1[l])),background_10KB)$statistic),#wilcox test, x=(values of first not top window for the second species), y=(10000 p_values in first variable that are not in top windows for first species)
                     P=wilcox.test (pull(filter(Values1,window5==window_list1[l])),background_10KB)$p.value,
                     Mean_test=mean(pull(filter(Values1,window5==window_list1[l]))),
                     N_samples=length(pull(filter(Values1,window5==window_list1[l]))[!is.na(pull(filter(Values1,window5==window_list1[l])))]),
                     Mean_BG=mean(background_10KB),
                     N_BG=length(background_10KB[!is.na(background_10KB)]),
                     Count=N_window)
        }
        temporary_result<- Result %>% mutate(Pred_P=(2* W - N_BG * N_samples)/sqrt(N_BG * N_samples*(N_BG+N_samples+1)/3))
        
        rm(Values1,window_list1,background_10KB)
        
        Emp_result<-foreach(m=1:nrow(temporary_result), .combine='rbind', .errorhandling='stop') %dopar% {
          data.frame(Emp_p = 1-(sum(as.numeric(temporary_result[m,13]) > pull(nulldist,Zscore))/nrow(nulldist)),Window=temporary_result[m,5])}
        
        Emp_result%>% left_join(.,Annotation)-> Emp_result
        
        left_join(temporary_result,Emp_result) %>% filter(!is.na(Emp_p)) -> Final_result
        fwrite(Final_result, file = paste0(result_directory,"baypass_nullw_result_recom_bin"), append = T,sep = "\t")
        rm(nulldist,Result,temporary_result,Emp_result,Final_result,N_window)}
      else {
        data.frame(variable=variable_list[j],quantile=quantile_list[z],Taxa1=combination_Taxa[i,1],Taxa2=combination_Taxa[i,2]) %>%
          fwrite(file = paste0(result_directory,"zero_top_baypass"), append = T,sep = "\t")
        rm(background_10KB,window_list,Values,window_list1,N_window,Values1)}
    }#variables loop
  }#quantile loop 
}#Taxas loop

all_vars<-read.table(paste0(result_directory,"baypass_nullw_result_recom_bin"), header = TRUE,sep="\t")
var_types<-unique(all_vars[,3])
out_res <- NULL
q= 0.05
for (i in 1:length(var_types)){
  aa<-all_vars[all_vars[,3]==var_types[i],]
  L = length(aa[,14])
  aa_sorted<-aa[order(aa[,14]),]
  w = aa_sorted[aa_sorted[,14] < q * ((1:L) / L),]
  out_res <-rbind(out_res,filter(w,Variable==var_types[i]))}
write.table(out_res,paste0(result_directory,"baypass_nullw_result_recom_bin_FDR"),col.names= TRUE, row.names = FALSE,sep="\t")

########
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)


#copy this directory in server and change the addredd below
"/data/users/mjahani/Convergent/Baypass/files/combination/"-> tmp_add
"/moonriseNFS/Null_W_Result/files/ID_window_final" -> annot
"/moonriseNFS/Null_W_Result/dbinom_binning/Baypass/"-> result_directory
fread(paste0(tmp_add,"Taxa_var_count_cutoff"), header = TRUE,sep="\t")-> data

options(scipen = 9999)
a<-c(seq(1,nrow(data),100000),nrow(data)+1)
registerDoParallel(cores=200)
for (j in 1:(length(a)-1)) {
  dbinom_score <- foreach(i= a[j]:(a[j+1]-1), .combine='rbind' , .errorhandling='stop') %dopar% {
    data[i,] %>%
      mutate(score= sum(dbinom(as.numeric(data[i,3]):as.numeric(data[i,4]),as.numeric(data[i,4]),as.numeric(data[i,5]))))
  }
  fwrite(dbinom_score, file = paste0(result_directory,"temp"), append = T,sep = "\t")
  rm(dbinom_score)
}


data.frame(species = c("Annuus","Argophyllus","petfal","petpet"),
           Variable = c("fake","fake","fake","fake"),
           window5 = c("fake","fake","fake","fake"),
           chrom = c(0,0,0,0),
           start_window = c(0,0,0,0),
           end_window = c(0,0,0,0),
           log_score = c(0,0,0,0)) -> fake_taxa

fread(paste0(result_directory,"temp"),header = T, fill=T) %>%
  mutate(log_score= log10(score)*-1) %>%
  separate(Taxa_var,into=c("species","Variable"),sep=":",extra="warn",remove = T) %>%
  separate(window5,into=c("chrom","S_E"),sep=":",extra="warn",remove = F) %>%
  separate(S_E,into=c("start_window","end_window"),sep="-",extra="warn",remove = T) %>%
  select(-snp_count,-outlier_count,-score,-expected) -> dbinom

dbinom  %>% distinct(Variable) %>% pull(Variable)-> var

Annotation<-fread(annot,header =T)  %>%
  select(window_id=window,window5) %>%
  distinct(window5, .keep_all = T) 

foreach(j= 1:length(var), .combine='rbind' , .errorhandling='stop') %dopar% {
  dbinom %>%
    filter(Variable==var[j]) %>%
    distinct(window5,species,.keep_all = T) %>%
    rbind(.,fake_taxa) %>%
    spread(species,log_score) %>%
    filter(window5 != "fake") %>%
    left_join(.,Annotation) %>%
    select(Annuus,Argophyllus,petfal,
           petpet,chrom,start_window,end_window,window_id,variable=Variable) %>%
    arrange(chrom,start_window,end_window) %>%
    fwrite(paste0(result_directory,"dbinom_score_baypass","_",var[j]),sep = "\t",na = "NA",quote = F)
}


library(tidyverse)
library(data.table)
fread("/moonriseNFS/Null_W_Result/binning/baypass_nullw_result_recom_bin",sep="\t") %>%
  distinct(Taxa1,Taxa2,Variable,quantile) %>% nrow()/ 2430

library(tidyverse)
library(data.table)
fread("/moonriseNFS/Null_W_Result/binning/dbinom/tmep",sep="\t") %>% nrow() / 22765966



fread("/moonriseNFS/Null_W_Result/binning/baypass_nullw_result_recom_bin",sep="\t") %>%
  distinct(Taxa1,Taxa2,Variable,quantile) %>%
  group_by(Taxa1,Taxa2,quantile) %>% tally() -> a

