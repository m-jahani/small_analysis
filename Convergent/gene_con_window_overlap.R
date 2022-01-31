library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)

data.frame(chrom=c("Ha412HOChr06","Ha412HOChr06","Ha412HOChr06","Ha412HOChr06","Ha412HOChr06","Ha412HOChr06","Ha412HOChr06","Ha412HOChr14","Ha412HOChr15"),
           chromStart=c(134215323,131658087,132943295,130707640,130724096,130727016,130729423,49106499,3315891),
           chromEnd=c(134218439,131660666,132944794,130710100,130726290,130729207,130731763,49108802,3319135),
           name=c("Ha412HOChr06g0281671","Ha412HOChr06g0281461","Ha412HOChr06g0281571","Ha412HOChr06g0281371_1","Ha412HOChr06g0281371_2","Ha412HOChr06g0281381_1","Ha412HOChr06g0281381_2","HaFT4","HaMYB111")) %>%
  mutate(chromStart=as.numeric(chromStart)-500,chromEnd=as.numeric(chromEnd)+500) -> marco_gene
# chrom chromStart chromEnd                  name
# Ha412HOChr06  134215323 134218439 Ha412HOChr06g0281671(HaFT1)
# Ha412HOChr06  131658087 131660666 Ha412HOChr06g0281461(HaFT2)
# Ha412HOChr06  132943295 132944794 Ha412HOChr06g0281571(HaFT3)
# Ha412HOChr06  130707640 130710100 Ha412HOChr06g0281371(HaFT5) 
# Ha412HOChr06  130724096 130726290 Ha412HOChr06g0281371(HaFT5)
# Ha412HOChr06  130727016 130729207 Ha412HOChr06g0281381(HaFT6)
# Ha412HOChr06  130729423 130731763 Ha412HOChr06g0281381(HaFT6)
# Ha412HOChr14  49106499  49108802 HaFT4
# Ha412HOChr15 3315891  3319135 HaMYB111


#gff3="/moonriseNFS/HA412/dovetail_assembly/genes/Ha412HOv2.0-20181130_genes_v2.1.gff3.gz"
gff3="/home/mjahani/Ha412HOv2.0-20181130_genes_v2.1.gff3" #uncompressd
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
  mutate(chromStart=as.numeric(chromStart)-range,chromEnd=as.numeric(chromEnd)+range) %>%
  rbind(.,marco_gene)


fread("/data/users/mjahani/Null_W_Result/files/convergent_2binning",sep = "\t",header=T) %>%
  separate(window5, into = c("Chr","start_end"), sep=":", remove = F)  %>%
  separate(start_end, into = c("start","end"), sep="-", remove = T) -> convergent

registerDoParallel(cores=48)
window_gene<-foreach(i=1:nrow(convergent), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(bed,
                       chrom == as.character(convergent[i,10]) &
                         ((chromStart >= as.numeric(convergent[i,11]) & chromStart <= as.numeric(convergent[i,12])) | 
                            (chromEnd >= as.numeric(convergent[i,11]) & chromEnd <= as.numeric(convergent[i,12])) | 
                            (chromStart < as.numeric(convergent[i,11]) & chromEnd > as.numeric(convergent[i,12])))),
                gene_start=chromStart,gene_end=chromEnd,gene_id=name),
         data=as.character(convergent[i,1]),
         type=as.character(convergent[i,2]),
         analysis=as.character(convergent[i,3]),
         Taxa1=as.character(convergent[i,4]),
         Taxa2=as.character(convergent[i,5]),
         comparison=as.character(convergent[i,6]),
         variable=as.character(convergent[i,7]),
         window_name=as.character(convergent[i,8]),
         window5=as.character(convergent[i,9]),
         Chr=as.character(convergent[i,10]),
         window_start=as.numeric(convergent[i,11]),
         window_end=as.numeric(convergent[i,12]))}

window_gene %>%
  select(data,
         type,
         analysis,
         Taxa1,
         Taxa2,
         comparison,
         variable,
         window_name,
         window5,
         Chr,
         window_start,
         window_end,
         gene_start,
         gene_end,
         gene_id) %>%
  fwrite("/data/users/mjahani/conv_window_gene_overlap/CONV_WIN_GENE_overlap",sep="\t",col.names = T)


#equivalant genes in arabidopsis
# Gene	Ha412	XRQv2_ID	
# HaFT1	Ha412HOChr06g0281671	HanXRQr2_Chr06g0274231	AT1G65480.1
# HaFT2	Ha412HOChr06g0281461	HanXRQr2_Chr06g0274071	AT1G65480.2
# HaFT3	Ha412HOChr06g0281571	HanXRQr2_Chr06g0274161	AT1G65480.1
# HaFT4	HaFT4	HanXRQr2_Chr14g0626821	AT1G65480.2
# HaFT5	Ha412HOChr06g0281371	HanXRQr2_Chr06g0274011	AT1G65480.1
# HaFT6	Ha412HOChr06g0281381	HanXRQr2_Chr06g0274021	AT1G65480.2
# HaMYB111	HaMYB111	HanXRQr2_Chr15g0672941	AT2G47460.1
#extra genes from marco that ha412V2.1 annotation does not include them
data.frame(Ha412HOv2.0_annot2.1=c("Ha412HOChr06g0281671",
            "Ha412HOChr06g0281461",
            "Ha412HOChr06g0281571",
            "HaFT4",
            "Ha412HOChr06g0281371",
            "Ha412HOChr06g0281381",
            "HaMYB111"),
  arabidopsis=c("AT1G65480.1",
                "AT1G65480.2",
                "AT1G65480.1",
                "AT1G65480.2",
                "AT1G65480.1",
                "AT1G65480.2",
                "AT2G47460.1"))-> MARCO

#genes that overlaps with convergent windows
fread("/data/users/mjahani/conv_window_gene_overlap/CONV_WIN_GENE_overlap",sep = "\t",header=T) %>%
  rename(Ha412HOv2.0_annot2.1=gene_id,
         Ha412HOv2.0_annot2.1_start=gene_start,
         Ha412HOv2.0_annot2.1_end=gene_end)-> window_gene

#arabifopsis against Ha412v2 BLAST hits
fread("/data/users/kaichi/GO/Ha412v2protein2_Ath_blastp.besthit",sep = "\t",header=F) %>%
  separate(V1,into = c("Ha412HOv2.0_annot2.1","remove"),sep = "-",remove = T) %>%
  select(Ha412HOv2.0_annot2.1,
         arabidopsis=V2) %>%
  rbind(.,
        MARCO)-> Ha412v2protein2_Ath_blastp.besthit

#other sunflower annotation and refrences genes
fread("/moonriseNFS/HA412/dovetail_assembly/genes/archive/Ha412HOv2.0-20181130_genes_id-map_v1.2.tsv",sep = "\t",header=T) %>%
  select(Ha412HOv2.0_annot2.1=Ha412HOv2.0,
         HanXRQr1.0,
         HanXRQv2) -> Sunflower_annotations



window_gene %>%
  rename(Ha412HOv2.0_annot2.1=gene_id) %>%
  left_join(.,Ha412v2protein2_Ath_blastp.besthit) %>%
  left_join(.,Sunflower_annotations) %>%
  select(data,
         type,
         analysis,
         Taxa1,
         Taxa2,
         comparison,
         variable,
         window_name,
         window5,
         Chr,
         window_start,
         window_end,
         Ha412HOv2.0_annot2.1_start=gene_start,
         Ha412HOv2.0_annot2.1_end=gene_end,
         Ha412HOv2.0_annot2.1,
         arabidopsis,
         HanXRQr1.0,
         HanXRQv2) %>%
  fwrite("/data/users/mjahani/conv_window_gene_overlap/CONV_WIN_GENE_overlap_Ath_HanXRQr1-2_Ha412HOv2_NA",sep="\t",col.names = T,na="NA")

#number of genes missing  
arabidopsis 30348
HanXRQv2  2832
HanXRQr1.0 3161


