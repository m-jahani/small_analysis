#run as (example): Rscript ~/input/NULL_LD.R small $SLURM_ARRAY_TASK_ID /scratch/mjahani/
#Mojtaba jahani  July 2020

suppressMessages(library(tidyverse))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(data.table))


args = commandArgs(trailingOnly = TRUE)

input_directory <- "~/input/"
result_directory <- "~/scratch/output/"
temporary_directory <- as.character(args[3])


list.files(path = input_directory)-> ped_file_list
gsub(".gz","",ped_file_list)-> ped_file_list
ped_file_list[grepl(".map",ped_file_list)]-> ped_file_list

fread(paste0(input_directory,"filesize"),
      sep = "\t",header = F) %>%
  filter(V1==as.character(args[1])) %>%
  arrange(V2) %>%
  slice(as.numeric(args[2]):as.numeric(args[2])) %>%
  select(V2)-> run_num
run_num[1,1]-> i
#suppressMessages(for (i in as.numeric(args[1]):as.numeric(args[1])) {

  Sys.time()-> start.time
  system(paste0("pigz -d ",
         input_directory,
         gsub(".map","",ped_file_list[i]),
         "*.gz"))

  data.frame(number= i,
             array=as.numeric(args[2]),
             tmp=gsub(".recode.vcf.map","",ped_file_list[i]),
             initial_snp_count= "Calculating",
             final_snp_count= "Calculating",
             window_pair_count= "Calculating",
             status="running",
             start_at=start.time,
             end_at="Calculating",
             time="Calculating") %>%
    separate(tmp,into=c("species","chr","data"),sep="_",remove=T) %>%
    select(number,array,data,species,chr,initial_snp_count,final_snp_count,window_pair_count,start_at,end_at,time,status) %>%
    fwrite(paste0(result_directory,
                  "report"),sep = "\t",col.names= F,append = T)
  
  # fread(paste0(result_directory,
  #              "report"),sep = "\t",header = F) %>%
  #   filter(V1 != "number") %>%
  #   rbind(.,c("number",
  #             "array",
  #             "data",
  #             "species",
  #             "chr",
  #             "initial_snp_count",
  #             "final_snp_count",
  #             "window_pair_count",
  #             "start_at",
  #             "end_at",
  #             "time",
  #             "status")) %>%
  #   arrange(desc(V1)) %>%
  #   fwrite(paste0(result_directory,
  #                 "report"),sep = "\t",col.names= F)
  
  system(paste0("plink --file ",
                input_directory,
                gsub(".map","",ped_file_list[i]),
                " --r2 --ld-window-kb 999999999 --ld-window 999999999 --ld-window-r2 0 --allow-extra-chr --out ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i])))
  
  system(paste0("pigz ",
                input_directory,
                gsub(".map","",ped_file_list[i]),
                "*"))

  system(paste0("wc -l ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i]),
                ".ld > ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i]),
                ".count &"))
  
  system(paste0("pigz ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i]),".ld"))
  
  system(paste0("pigz -dc ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i]),
                ".ld.gz | parallel --block 300M --pipe ' pigz > ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i]),
                "_part_{#}.gz'"))
  
  fread(paste0(input_directory,"all_id_gene_window_number_maf.03"),sep = "\t",header = T) %>%
    distinct(ID,window5) %>%
    separate(ID,into = c("chrom","pos"), sep=":", remove = F) %>%
    select(chrom,ID,window5) %>%
    filter(chrom==gsub("(_gwa.recode.vcf.map|_env.recode.vcf.map)","",gsub("[A-z]{6,11}_","",ped_file_list[i]))) %>%
    mutate(ID=gsub("Ha412HOChr","C",ID)) %>%
    select(ID,window5) -> ID_WIN_CHR
  
  
  list.files(path = temporary_directory) -> file_list
  file_list[grepl(gsub(".recode.vcf.map","",ped_file_list[i]),file_list)] -> file_list
  file_list[grepl("_part_",file_list)]-> file_list
  registerDoParallel(cores=28)
  foreach(z=1:length(file_list), .combine ='rbind', .errorhandling='stop') %dopar% {
    fread(paste0(temporary_directory,file_list[z]),
          header = F,showProgress=F) %>%
      select(SNP_A=V3,SNP_B=V6,R2=V7) %>%
      left_join(., rename(ID_WIN_CHR,SNP_A=ID,WIN_A=window5)) %>%
      filter(!is.na(WIN_A)) %>%
      left_join(., rename(ID_WIN_CHR,SNP_B=ID,WIN_B=window5)) %>%
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
      fwrite(paste0(temporary_directory,
                    gsub(".recode.vcf.map","",ped_file_list[i]),".LD_WINDOW"),sep = "\t",col.names= F,append = T)}
  rm(ID_WIN_CHR,file_list)
  system(paste0("rm ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i]),
                ".ld.gz"))
  system(paste0("rm ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i]),
                ".log"))
  system(paste0("rm ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i]),
                ".nosex"))
    system(paste0("rm ",
                temporary_directory,
                gsub(".recode.vcf.map","",ped_file_list[i]),
                "_part_*.gz"))
    
    fread(paste0(temporary_directory,
                 gsub(".recode.vcf.map","",ped_file_list[i]),".LD_WINDOW"),header = F) %>%
      rename(data=V1,Taxa=V2,window_A=V3,window_B=V4,LD=V5,snp_count=V6) %>%
      mutate(mult=(LD*snp_count)) %>%
      group_by(data,Taxa,window_A,window_B) %>%
      summarise(n_snp=sum(snp_count),multip=sum(mult)) %>%
      ungroup() %>% 
      mutate(mean_LD=(multip/n_snp)) %>%
      select(data,Taxa,window_A,window_B,LD=mean_LD,snp_count=n_snp) -> final_result
      fwrite(final_result,paste0(result_directory,
                    "window_pairwise_LD_result"),
             sep = "\t", col.names= F, append = T)

    read.table(paste0(temporary_directory,
                      gsub(".recode.vcf.map","",ped_file_list[i]),
                      ".count"),header = F)-> N
    Sys.time()-> end.time 
    end.time-start.time -> time_duration
    
    fread(paste0(result_directory,
                 "report"),sep = "\t",header = F) %>%
      filter(as.numeric(V1) != as.numeric(i)) %>%
      #rbind(.,report) %>%
      fwrite(paste0(result_directory,
                    "report"),sep = "\t",col.names= F)
    
    fread(paste0(result_directory,
                 "report"),sep = "\t",header = F) %>%
      mutate(sle=paste(V4,V5,V3,sep = "_")) %>%
      filter(sle != as.character(gsub(".recode.vcf.map","",ped_file_list[i]))) %>%
      select(-sle) %>%
      fwrite(paste0(result_directory,
                    "report"),sep = "\t",col.names= F)
    
    data.frame(number= i,
               array=as.numeric(args[2]),
               tmp=gsub(".recode.vcf.map","",ped_file_list[i]),
               initial_snp_count=N[1,1]-1,
               final_snp_count=(sum(final_result$snp_count)),
               start_at=start.time,
               end_at=end.time,
               window_pair_count=nrow(distinct(final_result,window_A,window_B))) %>%
      separate(tmp,into=c("species","chr","data"),sep="_",remove=T) %>%
      mutate(test=(initial_snp_count-final_snp_count)) %>%
      mutate(status=(ifelse(as.numeric(test)==0,"Done","Failed"))) %>%
      mutate(time=time_duration) %>%
      select(number,array,data,species,chr,initial_snp_count,final_snp_count,window_pair_count,start_at,end_at,time,status) %>%
      fwrite(paste0(result_directory,
                    "report"),sep = "\t",col.names= F,append = T)

    rm(final_result,report,N)
    system(paste0("rm ",
                  temporary_directory,
                  gsub(".recode.vcf.map","",ped_file_list[i]),
                  ".count"))
    system(paste0("rm ",
                  temporary_directory,
                  gsub(".recode.vcf.map","",ped_file_list[i]),
                  ".LD_WINDOW"))
  #  })



