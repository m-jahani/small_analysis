#Gene Ontology pipeline for Ha412v2.1
#GO terms extracted from TAIR homologs of Ha412v2.1 genes and InterPro terms from the Ha412v2.1 GFF
#Threshold for gene homologs: Blast e-value <= 1e-05

#Total genes in Ha412v2.1 = 49228

#Merged annotations (XRQ, TAIR, and InterPro with all ancestral annotations):
#Total genes with MF GO annotations = 32941 (66.92%)
#Total genes with BP GO annotations = 29171 (59.26%)
#Total genes with CC GO annotations = 30348 (61.65%)

#XRQ homolog blast2GO annotations with <3 parent node distances included only:
#Total genes with MF GO annotations = 28344 (57.58%)
#Total genes with BP GO annotations = 24203 (49.17%)
#Total genes with CC GO annotations = 17555 (35.66%)

#TAIR homolog annotations without added parent terms:
#Total genes with MF GO annotations = 23221 (47.17%)
#Total genes with BP GO annotations = 22750 (46.21%)
#Total genes with CC GO annotations = 25537 (51.87%)

######################### SET UP #########################

# (1) Load dependancies. Library dependancies: stats package, GO.db package.
library(stats)
library(GO.db)

# (2) Set Ha412v2.1_GO folder as working directory
setwd("/data/users/mjahani/GO_analysis/Ha412_v2.1_GO/")

# (3) Load necessary functions:
topGO_inverseList <- function(l){
  rId <- unlist(l, use.names = FALSE)
  lId <- rep(names(l), sapply(l, length))
  return(split(lId, rId))
}#Function borrowed from the topGO library (topGO::inverseList()), converts data frames into nested lists

stop_quietly <- function(){
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}#stop functions without error messages

run_helianthus_GO <- function(query_vector, library="Genome", ontology="BP", include.gene.names=T, p.adj="fdr", Source="Merged", subset_library=NULL, filter=0.05, statistic="Fisher", specific_terms_only=F, term_removal_distance=NULL){
  
  #First, input files based on ontology parameters:
  if(Source == "Merged"){
    input_file <- paste(ontology, "_v2.1_Tair_InterPro_XRQ.csv", sep="")
  }else if(Source == "XRQ"){
    input_file <- paste(ontology, "_v2.1_XRQ_GO.csv", sep="")
  }else if(Source == "TAIR"){
    input_file <- paste(ontology, "_v2.1_TAIR_GO.csv", sep="")
  }else{
    print("Error: 'Source' option is not recognized.")
    if(! class(Source) == "character"){
      print("'Source' option should be a string.")
      stop_quietly()
    }
  }
  
  #read in gene 2 GO term file:
  gene2GO <- read.csv(input_file, header = F, fill = T, stringsAsFactors = F, blank.lines.skip = TRUE, row.names = 1)
  #If the library is a subset of the genome, shrink the gen2GO DB to reflect this:
  if(library=="Subset"){
    gene2GO <- gene2GO[rownames(gene2GO)%in%as.vector(subset_library),]
    #print an error if the subset library does not contain any genes
    if(nrow(gene2GO) == 0){print("Subset library contains NO genes. Likely input error. Ensure that 'subset_library' is a vector of gene names in the format 'Ha412HOChrXXgXXXXXXXX.'")
      stop_quietly()}
  }
  #Reformat into GO2gene file:
  temp <- as.data.frame(t(gene2GO))
  GO2gene <-  topGO_inverseList(temp)
  rm(temp)
  
  #remove non-specific GO terms if specified to do so:
  if(specific_terms_only==T){
    if(term_removal_distance==1 | term_removal_distance==2 | term_removal_distance==3){
      #input specificity file
      spec_input <- paste(ontology, "_term_specificity.txt", sep ="")
      spec <- read.delim(spec_input)
      #extract GO terms to remove:
      remove_GO_terms <- spec[spec$distance <= term_removal_distance,]
      remove_GO_terms <- as.character(remove_GO_terms$Term)
      #remove extracted GO terms:
      GO2gene <- GO2gene[names(GO2gene) %in% remove_GO_terms == FALSE]
    } else {
      print("term_removal_distance must equal 1, 2, or 3.")
      stop_quietly()
    }
  }
  
  #Read query genes:
  if(! class(query_vector) == "character"){
    gene_vector <- as.character(query_vector)
    message <- paste("Converting query_vector from", class(query_vector), "to character.")
    print(message)
  }
  #Do not allow repeated genes to be counted multiple times
  unique_genes <- unique(as.vector(query_vector))
  #Subset the gene2GO genomic df into a gene2GO df for only the query genes:
  GO_Annots_Temp <- gene2GO[rownames(gene2GO)%in%unique_genes,]
  
  #extract all the GO terms in the query genes into a vector called "output_vect:
  for(row in 1:nrow(GO_Annots_Temp)){
    k <- unique(as.character(as.vector(GO_Annots_Temp[row,])))
    k <- k[grepl("GO", k)]  
    if(exists("output_vect") == F){
      output_vect <- k
    } else {
      output_vect <- c(output_vect, k)}
  }#for loop end
  
  #Make a frequency dataframe of all GO terms in the query genes called "w":
  GO.ID <- factor(output_vect)
  rm(output_vect, k)
  w = as.data.frame(table(GO.ID))
  
  #remove excluded GO terms from w if any are specified:
  if(specific_terms_only==T){
    kept <- names(GO2gene)
    w <- w[w$GO.ID %in% kept,]
  }
  
  #if no GO terms exist in the query genes (i.e., only unannotated genes), end the function here.
  if(nrow(w) == 0){
    message <- paste("No GO terms are annotated to any of the query genes for the", ontology, "ontology")
    print(message)
    stop_quietly()
  } 
  
  #Make columns for Actual # genes with GO term, total # genes with GO term in library, expected # genes within GO term, and p-value of enrichment
  w$GO.ID<- as.character(w$GO.ID)
  colnames(w)[2] <- "Actual"
  w <- subset(w, w$Actual > 1)
  #if no GO terms exist in the query genes (i.e., only unannotated genes), end the function here.
  if(nrow(w) == 0){
    message <- paste("No GO terms are annotated to any of the query genes for the", ontology, "ontology")
    print(message)
    stop_quietly()
  } 
  
  #total genes in Query with annotations:
  signif_genes = nrow(GO_Annots_Temp)
  #total genes in library:
  total_genes = nrow(gene2GO)
  
  #Find the expected number of genes to have that query number given the freqency of the GO term in the library:
  #Need some sub-functions to do this:
  #function for finding total # of genes with a GO ID
  total <- function(GO.ID){
    t <- GO2gene[GO.ID]
    t <- as.character(unlist(t))
    t <- t[t%in%as.vector(rownames(gene2GO))]
    Total <- length(t)
    return(Total)
  }
  #function for finding total expected genes in a subset of genes
  expected <- function(Total){
    exp <- (Total/total_genes)*signif_genes
    return(exp)
  }
  #fisher test:
  get_fisher <- function(w){
    mat <- matrix(as.numeric(unlist(w[c(2,5:7)])), ncol=2)
    f <- fisher.test(as.table(mat), alt="greater")
    return(c(w[1], f$p.value))
  }
  #Chi squared:
  get_Chi <- function(w){
    mat <- matrix(as.numeric(unlist(w[c(2,5:7)])), ncol=2)
    c <- chisq.test(as.table(mat))
    return(c(w[1], c$p.value))
  }
  
  #total = genes in genome with annotation
  w$Total <- sapply(w$GO.ID, total)
  #exp = (total/total_genes)*signif_genes
  w$Expected <- sapply(w$Total, expected)
  #only include samples where the expected value is smaller than the actual value
  w <- subset(w, (w$Actual > w$Expected))
  #temporary columns for the contigency table:
  w$allGO <- w$Total-w$Actual
  w$perGO <- signif_genes- w$Actual
  w$pernonSig <- (total_genes - (w$Total - w$Actual))
  #function for your statistical test:
  if(statistic == "Fisher"){
    p_vect <- t(apply(w, 1, get_fisher))} 
  else if(statistic == "Chi"){
    p_vect <- t(apply(w, 1, get_Chi))}
  else{print("Statistical test not recognized.")
    stop_quietly()}
  
  #Reformat output file:
  w <- merge(w, p_vect, by = "GO.ID")
  colnames(w)[8] <- "p_val"
  w$p_val <- as.numeric(as.character(w$p_val))
  w <- w[c(-5,-6,-7)]
  #Add adjusted p-values:
  w$p.adj <- p.adjust(w$p_val,method=p.adj)
  #Remove insignificant terms:
  w <- subset(w, w$p.adj <= filter)
  #Add GO term anrow=max(lengths(t))nnotations from GO.db:
  w$Term <- Term(w$GO.ID)
  
  if(include.gene.names == T & nrow(w) > 0){
    w$Gene_Names <- NA
    for(row in 1:nrow(w)){
      temp_go_term <- w$GO.ID[row]
      temp <- GO2gene[names(GO2gene)%in%temp_go_term]
      temp <- as.vector(unlist(temp))
      temp <- unique_genes[unique_genes%in%temp]
      temp <- unique(temp)
      temp <- paste(temp,collapse=",")
      w$Gene_Names[row] <- temp 
    }#for loop end
  }
  
  #Return the final dataframe:
  return(w)
}#Ha412 Gene ontology function. See below for parameters and functions


################ run_helianthus_GO parameters: ################

# > query_vector: A character vector of genes in which to search for GO enrichment

# > library: the genes the query_vector will be search against. Default = "Genome". Options = ("Genome", "Subset"). "Subset" used if only a selection of genes are to be used as the library. Use of "Subset" requires a vector of gene names to be input in the "subset_library" option. 

# > ontology: Options = "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Compartment). Default = "BP"

# > include.gene.names: Boolean. Default = T. Whether or not the names of genes with each GO term should be included in the output file.

# > p.adj: Options: see p.adjust.methods. Default = "fdr"

# > Source: Source of GO terms. Options = "XRQ" (GO terms from XRQ homolog blast2GO only), "TAIR" (GO terms from TAIR homologs), or "Merged" (GO terms from TAIR homologs, XRQ homologs, and InterPro terms). Default = "Merged". The "XRQ" library features GO with a maximum of 2 parent node steps away from directly annotated GO terms. The "TAIR" library contains all GO terms directly annotated to Tair homologs of Ha412 genes, without ancestral nodes added. Due to the variation in annotation methods of Tair contributers, annotations may follow varying paradigms and thus may not be entirely representitive. In the "Merged" annotations, ALL ancestral terms are used as the mixed annotation paradigms for TAIR, XRQ, and InterPro cause a lot of variability and including all terms is the best way to minimize resultant noise. The "Merged" GO library thus contains far more GO terms than the others, and will be more likely to give false positives. When using this library, use of a more stringent filtering threshold, statistical test, and p.adjustment method would be wise. 

# > subset_library: If entire genome is not being used as the background library, input a character vector of gene names to be used as the background library. Subset library must contain ALL genes in the query, otherwise results may be problematic.

# > filter: p-adjust threshold above which GO terms should be discarded from the output file. Default = 0.05.

# > statistic: Which test to use to generate a p-value. Options = "Fisher" (Fisher's Exact Test), "Chi" (Chi-squared). Default = "Fisher", Fisher's test will give more conservative p-values than Chi-squared. 

# > specific_terms_only: Boolean. Whether or not to remove non-specific GO terms from the analysis. Requires and input for term_removal_distance. Deafult = FALSE.

# > term_removal_distance: numeric, one of (1,2,3). When using the specific_terms_only = TRUE option, term_removal_distance allows the user to determine how specific output terms will be. 1, 2, and 3 refer to term distances from the overall ontology annotation. For example, term_removal_distance=1 in a "BP" ontology would remove all GO terms that are direct children of the "biological_process" annotation. term_removal_distance=2 removes direct children of the "biological_process" annotation, and the children of direct children of of the "biological_process" annotation. 1 removes the least specific terms, while 3 removes more specific terms. Specifically useful when using the "merged" library which is prone to false positives. Default = NULL.


############### Examples: ###############

#query some random genes against the genome for biological process:
#query <- dplyr::sample_n(all_genes, 1000)
#query <- query$Gene
#GO_output <- run_helianthus_GO(query_vector=query, include.gene.names=F)
#GO.ID        Actual Total Expected   fisher_p   p.adj      Annotation
#GO:0009416      3   166 0.3322745 0.004449566 0.04152929 response to light stimulus
#GO:0016192      2   122 0.2442017 0.024757417 0.09041839 vesicle-mediated transport

#query some random genes against a subset of the genome:
#lib <- dplyr::sample_n(all_genes, 10000)
#lib <- lib$Gene
#lib <- unique(c(query, lib))
#GO_output <- run_helianthus_GO(query_vector=query, library="Subset", subset_library=lib, include.gene.names=F)
#GO.ID        Actual Total Expected    fisher_p   p.adj        Annotation
#GO:0009416      3     8 0.7419355 0.023704985 0.1738817 response to light stimulus
#GO:0016192      2     2 0.1854839 0.007059498 0.1738817 vesicle-mediated transport

#include gene names in output file:
#GO_output <- run_helianthus_GO(query_vector=query)
#GO.ID        Actual Total Expected   fisher_p   p.adj         Annotation                 Gene_Names
#GO:0009416      3   166 0.3322745 0.004449566 0.04152929 response to light stimulus Ha412HOChr03g00007336,Ha412HOChr10g00024597,Ha412HOChr12g00051837
#GO:0016192      2   122 0.2442017 0.024757417 0.09041839 vesicle-mediated transport                       Ha412HOChr12g00031889,Ha412HOChr04g00012162

#Chi-squared test:
#GO_output <- run_helianthus_GO(query_vector=query, include.gene.names=F, statistic = "Chi")

#p-value cut-off of 0.001:
#GO_output <- run_helianthus_GO(query_vector=query, filter = 0.001, include.gene.names=F)

#Molecular function ontology:
#GO_output <- run_helianthus_GO(query_vector=query, ontology = "MF", include.gene.names=F)

#only look at more specific GO terms:
#GO_output <- run_helianthus_GO(query_vector=query, include.gene.names=F, specific_terms_only=T, term_removal_distance = 3)

#XRQ library:
#GO_output <- run_helianthus_GO(query_vector=query, Source = "XRQ", include.gene.names=F)


library(tidyverse)
library(data.table)

fread("/data/users/mjahani/convergent_2binning_gene", sep = "\t",header = T) %>%
  filter(!is.na(Gene)) %>% 
  mutate(dtacv=paste(data,type,analysis,comparison,variable,sep = "-")) %>%
  select(dtacv,Gene) -> convergent_genes

convergent_genes %>%
  # group_by(dtacv) %>%
  # tally() %>%
  # ungroup() %>%
  # filter(n>1)%>%
  distinct(dtacv) %>%
  pull(dtacv) -> DTACY

# GO_result <- NULL
# for (i in 66:68) {
#   convergent_genes %>%
#     filter(dtacv==DTACY[i]) %>%
#     distinct(Gene) %>%
#     pull(Gene) -> GENES
#   run_helianthus_GO(query_vector=GENES) %>% 
#     mutate(dtacv=DTACY[i]) %>%
#     rbind(.,GO_result) -> GO_result
# }

GO_result <- NULL
for (i in 1:length(DTACY)) {
  convergent_genes %>%
    filter(dtacv==DTACY[i]) %>%
    distinct(Gene) %>%
    pull(Gene) -> GENES
  tryCatch(run_helianthus_GO(query_vector=GENES), error=function(e) NULL) -> data
  if (is.null(data)) {
    data.frame(GO.ID="No_GO_terms",Actual=NA,Total=NA,Expected=NA,p_val=NA,p.adj=NA,Term=NA,Gene_Names=NA,dtacv=DTACY[i]) %>%
  rbind(.,GO_result) -> GO_result
  } else if (as.numeric(nrow(data))==0) {
    data.frame(GO.ID="zero_line",Actual=NA,Total=NA,Expected=NA,p_val=NA,p.adj=NA,Term=NA,Gene_Names=NA,dtacv=DTACY[i]) %>%
      rbind(.,GO_result) -> GO_result
  } else {
    data %>%
    mutate(dtacv=DTACY[i]) %>%
      rbind(.,GO_result) -> GO_result
    }
  rm(GENES,data)
}

  
GO_result %>%
  filter(is.na(Gene_Names)) %>%
  mutate(dtacv=gsub("Leaf_height_mid-width","Leaf_height_mid_width",dtacv)) %>%
  mutate(dtacv=gsub("Leaf_width_mid-height","Leaf_width_mid_height",dtacv)) %>%
  separate(dtacv,into = c("data","type","analysis","comparison","variable"),sep = "-",remove = T) %>% 
  mutate(variable=gsub("Leaf_height_mid_width","Leaf_height_mid-width",variable)) %>%
  mutate(variable=gsub("Leaf_width_mid_height","Leaf_width_mid-height",variable)) %>%
  select(data,
         type,
         analysis,
         comparison,
         variable,
         GO.ID,
         Actual,
         Total,
         Expected,
         p_val,
         p.adj,
         Term,
         Gene_Names) %>%

