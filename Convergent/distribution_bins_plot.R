library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
files<-c("baypass_nullw_result","GWAScorrected_nullw_result","GWASuncorrected_nullw_result","spearman_nullw_result")
analy<-c("baypass","GWAS_corrected","GWAS_uncorrected","spearman")
datt<-c("env","gwa","gwa","env")
out_res <- NULL
for (j in 1:length(files)) {
  all_vars<-read.table(paste0("/moonriseNFS/Null_W_Result/null_w_binning/",files[j]), header = TRUE,sep="\t") %>%
    mutate(analysis=analy[j],data=datt[j])
  var_types<-unique(all_vars[,3])
  q= 0.05
  for (i in 1:length(var_types)){
    aa<-all_vars[all_vars[,3]==var_types[i],]
    L = length(aa[,14])
    aa_sorted<-aa[order(aa[,14]),]
    w = aa_sorted[aa_sorted[,14] < q * ((1:L) / L),]
    out_res <-rbind(out_res,filter(w,Variable==var_types[i]))} 
  rm(var_types,all_vars)
}


convert <-data.frame(Variable=c(c("OM","P1","P2","BICARB","K","MG","CA","Sodium",'PH',"CEC",
         "PERCENT_K","PERCENT_MG","PERCENT_CA","PERCENT_NA","SOL_SALTS"),
c("latitude","longitude","llevation","MAT","MWMT","MCMT",
         "TD","MAP","MSP","AHM","SHM","DD_0",
         "DD5","DD_18","DD18","NFFD","bFFP","eFFP",
         "FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH")),type=c(rep("soil",15),rep("climate",26)))

out_res %>%
  left_join(.,convert) %>%
  mutate(type=as.character(type)) %>%
  mutate(type=replace_na(type,"phenotype")) %>%
  mutate(binning=1) %>%
  select(binning,data,type,analysis,Taxa1,Taxa2,variable=Variable,window=window_name) -> one_bin

fread("/moonriseNFS/mojtaba_tmp/LD/convergent_windows_all",header = T) %>%
  mutate(binning=0) %>%
  select(binning,data,type,analysis,Taxa1,Taxa2,variable,window)-> no_bin


add<-c("Baypass","GWAS_corrected","GWAS_uncorrected","Spearman")
files<-c("baypass_nullw_result_recom_bin_FDR","GWAScorrected_nullw_result_recom_bin_FDR","GWASuncorrected_nullw_result_recom_bin_FDR","spearman_nullw_result_recom_bin_FDR")
analy<-c("baypass","GWAS_corrected","GWAS_uncorrected","spearman")
datt<-c("env","gwa","gwa","env")
two_bin <- NULL

for (j in 1:length(files)) {
  read.table(paste0("/moonriseNFS/Null_W_Result/null_w_topcandidate_binning/",add[j],"/",files[j]), header = TRUE,sep="\t") %>%
    mutate(analysis=analy[j],data=datt[j]) %>% 
    rbind(.,two_bin)-> two_bin
}

two_bin %>%
  left_join(.,convert) %>%
  mutate(type=as.character(type)) %>%
  mutate(type=replace_na(type,"phenotype")) %>%
  mutate(binning=2) %>%
  select(binning,data,type,analysis,Taxa1,Taxa2,variable=Variable,window=window_name) -> two_bin
  
rbind(no_bin,one_bin,two_bin) -> all_convergent


fread("/moonriseNFS/mojtaba_tmp/recombination/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>%
  group_by(chr) %>%
  arrange(pos,.by_group = TRUE) %>%
  mutate(delta_cM=cM-lag(cM, default = first(cM)),
         delta_bp=pos-lag(pos, default = first(pos)),
         start=lag(pos)) %>% 
  ungroup() %>%
  mutate(recomb_rate=delta_cM/delta_bp) %>%
  select(chr,start,end=pos,recomb_rate) %>%
  filter(!is.na(start)) -> recombination


fread("/moonriseNFS/mojtaba_tmp/LD/all_id_gene_window_number_maf.03",sep = "\t",header = T) %>%
  select(window) %>%
  distinct(window) %>%
  pull(window)-> windowww


vanil_50K<-read.table(file ="/moonriseNFS/soudi/for_Emily/bins_good_5000",header = FALSE)
vanil_50K$V6<-NULL
vanil_50K$V1<-paste("window",vanil_50K$V1,sep = "_")
vanil_50K$V2 <- gsub("chr", "", vanil_50K$V2)
colnames(vanil_50K)<-c("window","chrom","start","end","width")
vanil_50K %>% 
  mutate(chrom=as.numeric(chrom)) %>%
  mutate(chr=ifelse(chrom>=10,paste0("Ha412HOChr",chrom),paste0("Ha412HOChr0",chrom))) %>%
  select(window,chr,start,end) -> vanil_50K

registerDoParallel(cores=40)
window_recombination<-foreach(i=1:nrow(recombination), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(vanil_50K,
                       chr == as.character(recombination[i,1])  &
                         start >= as.numeric(recombination[i,2])&
                         end <= as.numeric(recombination[i,3])),
                window),
         recomb_rate=as.numeric(recombination[i,4]))}

window_recombination %>%
  filter(window %in% windowww) -> window_recombination


window_recombination$quartile <- with(window_recombination, cut(recomb_rate, 
                                                  breaks=quantile(recomb_rate, probs=seq(0,1, by=0.2), na.rm=TRUE), 
                                                  include.lowest=TRUE))

convert <- data.frame(quartile=c("[1.15e-10,4.76e-10]","(4.76e-10,2.67e-08]",
                                 "(2.67e-08,2.88e-07]","(2.88e-07,8.6e-07]","(8.6e-07,5.44e-06]"),
                      quantiles=c("0-20","20-40","40-60","60-80","80-100"))

left_join(window_recombination,convert) %>% 
  mutate(Q=as.character(quantiles)) %>%
  select(-quartile,-quantiles)  -> window_recombination
rm(convert)

left_join(all_convergent,window_recombination) -> all3_convergent_recomb_quant
fwrite(all3_convergent_recomb_quant,"/moonriseNFS/Null_W_Result/files/all3_convergent_recomb_quant",sep = "\t")
#####Macbook
library(tidyverse)
library(data.table)
 fread("/Users/mojtabajahani/Downloads/all_3_convergent_bin_plots/all3_convergent_recomb_quant",sep="\t") -> convergents

 convergents %>%
   group_by(binning,Q) %>%
   summarize(window_count=n_distinct(window)) %>%
   ungroup()-> tmp
 
 tmp %>%
   group_by(binning) %>%
   summarise(total=sum(window_count)) %>%
   ungroup() %>%
   full_join(.,tmp) %>%
   mutate(proportion=window_count/total) %>%
   select(binning,Q,proportion) ->binning_window_count
 rm(tmp)
   
  
 convergents %>%
  group_by(binning,type,analysis,Q) %>%
  summarize(window_count=n_distinct(window)) %>%
  ungroup() -> tmp
 
 tmp %>%
   group_by(binning,type,analysis) %>%
   summarise(total=sum(window_count)) %>%
   ungroup() %>%
   full_join(.,tmp) %>%
   mutate(proportion=window_count/total) %>%
   select(binning,type,analysis,Q,proportion) -> type_analysis_window_count
 rm(tmp)   
 
 convergents %>%
   group_by(binning,analysis,Q) %>%
   summarize(window_count=n_distinct(window)) %>%
   ungroup() -> tmp
 
 tmp %>%
   group_by(binning,analysis) %>%
   summarise(total=sum(window_count)) %>%
   ungroup() %>%
   full_join(.,tmp) %>%
   mutate(proportion=window_count/total) %>%
   select(binning,analysis,Q,proportion) -> analysis_window_count
 rm(tmp) 
 
 convergents %>%
   group_by(binning,type,Q) %>%
   summarize(window_count=n_distinct(window)) %>%
   ungroup() -> tmp
 
 tmp %>%
   group_by(binning,type) %>%
   summarise(total=sum(window_count)) %>%
   ungroup() %>%
   full_join(.,tmp) %>%
   mutate(proportion=window_count/total) %>%
   select(binning,type,Q,proportion) -> type_window_count
 rm(tmp) 
   
 ggplot(binning_window_count,aes(x=Q, y=proportion,fill=as.character(binning))) + 
   geom_bar(stat="identity", color="black", position=position_dodge())+
   xlab("Recombination Rate Quantiles")+
   ylab("Proportion of Convergent Windows")+
   scale_fill_discrete(labels=c("Not binned","Null-W binned","Top candidate and Null_W binned"))+
   guides(fill=guide_legend(title="Null-W Test"))+
   ggtitle("All the convergent windows")+
   theme_classic()
 
 ggplot(analysis_window_count,aes(x=Q, y=proportion,fill=as.character(binning))) + 
   geom_bar(stat="identity", color="black", position=position_dodge())+
   xlab("Recombination Rate Quantiles")+
   ylab("Proportion of Convergent Windows")+
   scale_fill_discrete(labels=c("Not binned","Null-W binned","Top candidate and Null_W binned"))+
   guides(fill=guide_legend(title="Null-W Test"))+
   facet_wrap(~analysis) +
   theme_classic() 
 
 ggplot(type_window_count,aes(x=Q, y=proportion,fill=as.character(binning))) + 
   geom_bar(stat="identity", color="black", position=position_dodge())+
   xlab("Recombination Rate Quantiles")+
   ylab("Proportion of Convergent Windows")+
   scale_fill_discrete(labels=c("Not binned","Null-W binned","Top candidate and Null_W binned"))+
   guides(fill=guide_legend(title="Null-W Test"))+
   facet_wrap(~type) +
   theme_classic()   
 
 ggplot(type_analysis_window_count,aes(x=Q, y=proportion,fill=as.character(binning))) + 
   geom_bar(stat="identity", color="black", position=position_dodge())+
   xlab("Recombination Rate Quantiles")+
   ylab("Proportion of Convergent Windows")+
   scale_fill_discrete(labels=c("Not binned","Null-W binned","Top candidate and Null_W binned"))+
   guides(fill=guide_legend(title="Null-W Test"))+
   facet_wrap(type~analysis) +
   theme_classic()   
 
 binning_window_count%>% group_by(binning) %>% summarise(a=sum(proportion)) %>% ungroup() %>% distinct(a)

 analysis_window_count%>% group_by(binning,analysis) %>% summarise(a=sum(proportion)) %>% ungroup() %>% distinct(a)

 type_window_count%>% group_by(binning,type) %>% summarise(a=sum(proportion)) %>% ungroup() %>% distinct(a)

 type_analysis_window_count%>% group_by(binning,analysis,type) %>% summarise(a=sum(proportion)) %>% ungroup() %>% distinct(a)




