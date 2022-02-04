setwd("~/Dropbox/desktop/convergence_nullW")
options(stringsAsFactors = F)

#dd <- read.table ("~/Dropbox/desktop/convergence_nullW/sunflower_LD_clusters.csv",T,sep = ",")
#dd <- read.table ("~/Dropbox/desktop/convergence_nullW/LD_clustering_merged_0.95_5cM.csv",T,sep = ",")
dd <- read.table ("~/Dropbox/desktop/convergence_nullW/LD_clustering_0.9_5cM.csv",T,sep = ",")
#dd <- read.table ("~/Dropbox/desktop/convergence_nullW/LD_clustering_merged_0.9_5cM.csv",T,sep = ",")


inv <- read.table ("Ha412HO_inv.v3.inversions.regions.v1.txt",T)

inv <- inv[inv$species != "Niveus",]

invname <- inv$spe

invname[invname == "annuus"] <- "ann"
invname[invname == "argophyllus"] <- "arg"
invname[invname == "petiolaris"] <- "pet"

inv$invname <- invname

#ff <- dd[dd$analysis == "baypass",]
ff <- dd[dd$analysis != "baypass",]


int1 <- strsplit (as.character (ff$range), split = ":")
ff$low <- as.numeric (sapply (int1,"[[",1))
ff$high <- as.numeric (sapply (int1,"[[",2))

the_vars <- unique (ff$variable)

#barplot (ff$variable)

the_dirs2 <- ff$direction
the_dirs2[the_dirs2 == "Argophyllus_Annuus"] <- "Annuus_Argophyllus"
the_dirs2[the_dirs2 == "petfal_Annuus"] <- "Annuus_petfal"
the_dirs2[the_dirs2 == "petpet_Annuus"] <- "Annuus_petpet"
the_dirs2[the_dirs2 == "petfal_Argophyllus"] <- "Argophyllus_petfal"
the_dirs2[the_dirs2 == "petpet_Argophyllus"] <- "Argophyllus_petpet"
the_dirs2[the_dirs2 == "petpet_petfal"] <- "petfal_petpet"

the_dirs2[the_dirs2 == "Annuus_Argophyllus"] <- "ann_arg"
the_dirs2[the_dirs2 == "Annuus_petfal"] <- "ann_fal"
the_dirs2[the_dirs2 == "Annuus_petpet"] <- "ann_pet"
the_dirs2[the_dirs2 == "Argophyllus_petfal"] <- "arg_fal"
the_dirs2[the_dirs2 == "Argophyllus_petpet"] <- "arg_pet"
the_dirs2[the_dirs2 == "petfal_petpet"] <- "fal_pet"


#colours
the_dirs2_col <- ff$direction
the_dirs2_col[the_dirs2_col == "Argophyllus_Annuus"] <- "black"
the_dirs2_col[the_dirs2_col == "petfal_Annuus"] <- "black"
the_dirs2_col[the_dirs2_col == "petpet_Annuus"] <- "black"
the_dirs2_col[the_dirs2_col == "petfal_Argophyllus"] <- "black"
the_dirs2_col[the_dirs2_col == "petpet_Argophyllus"] <- "black"
the_dirs2_col[the_dirs2_col == "petpet_petfal"] <- "black"
the_dirs2_col[the_dirs2_col != "black"] <- "grey50"

ff$the_dirs2 <- the_dirs2
ff$the_dirs2_col <- the_dirs2_col

the_dirs <- unique (ff$direction)

uniq_the_dirs2 <- unique(the_dirs2)

#use this to match with inversion naming (no petfal/petpet)
temp1 <- strsplit (the_dirs2,split = "_")
the_dirs_inv1 <- sapply(temp1,"[[",1)
the_dirs_inv2 <- sapply(temp1,"[[",2)
the_dirs_inv1[the_dirs_inv1 == "fal"] <- "pet"
the_dirs_inv2[the_dirs_inv2 == "fal"] <- "pet"

ff$the_dirs_inv1 <- the_dirs_inv1
ff$the_dirs_inv2 <- the_dirs_inv2


the_res <- array (0, c (length (the_vars),length (the_dirs),length(the_dirs)))
the_res2 <- array (0, c (length (the_vars),length (the_dirs),length(the_dirs)))

chroms <- unique (ff$chromosome)

the_max <- tapply (ff$high,ff$chromosome,max)

#pdf("sunflower_LD_clusters_0.9cutoff.pdf")
#pdf("sunflower_LD_clusters_merged_0.95cutoff.pdf")
pdf("sunflower_LD_clusters_merged_0.9cutoff.pdf")
#pdf("sunflower_LD_clusters_baypass.pdf")

col_inv <- c("red","blue")

for (i in 1:length (the_vars)){

	sub1 <- ff[ff$variable == the_vars[i],]	

	#loop through the chromosomes:
	for (j in 1:length (chroms)){
	
		if (j == 1 | j == 10){
			par (mfcol = c (3,3))
			j_axis <- 1
		}
		
		sub2 <- sub1[sub1$chromosome == chroms[j],]
		
		this_max <- the_max[names(the_max) == chroms[j]]
		
		plot (1,1,xlim = c(0, this_max),ylim = c (0,1.2),col = "white", main = paste(the_vars[i],chroms[j],sep = "  "),yaxt = "n",ylab = "", xlab = "Position (Mbp)", cex.main = 0.5,xaxt = "n")

		axis(1,at = c(0,50000000,100000000,150000000,200000000,250000000), labels = c(0,50,100,150,200,250))

		if (j_axis < 4){
			par (las = 1)
			axis(2,labels = uniq_the_dirs2,at=((1:6)/6))
		} 
		j_axis <- j_axis + 1
		
		if (nrow (sub2) > 0){
		
		
		for (k in 1:nrow(sub2)){
			
			ypos <- which(uniq_the_dirs2 == sub2$the_dirs2[k])/length (uniq_the_dirs2)
			
			#check inversion chromsome, position, and species:
			inv2 <- inv[inv$chr == chroms[j],]
			check1 <- sub2$low[k] > inv$start & sub2$low[k] < inv$end & chroms[j] == inv$chr & (sub2$the_dirs_inv1[k] == inv$invname | sub2$the_dirs_inv2[k] == inv$invname)
			check2 <- sub2$high[k] > inv$start & sub2$high[k] < inv$end & chroms[j] == inv$chr & (sub2$the_dirs_inv1[k] == inv$invname | sub2$the_dirs_inv2[k] == inv$invname)
			check3 <- inv$start > sub2$low[k]  & inv$start < sub2$high[k] & chroms[j] == inv$chr & (sub2$the_dirs_inv1[k] == inv$invname | sub2$the_dirs_inv2[k] == inv$invname)
			check4 <- inv$end > sub2$low[k]  & inv$end < sub2$high[k] & chroms[j] == inv$chr & (sub2$the_dirs_inv1[k] == inv$invname | sub2$the_dirs_inv2[k] == inv$invname)
			
			checksum1 <- sum(check1 + check2 + check3 + check4)
			
			check_ind <- which ((check1 + check2 + check3 + check4) > 0)
			
			sub_inv <- inv[check_ind,]
			
			if (checksum1 != 0){
				
				for (bb in 1:nrow(sub_inv)){
					arrows(sub_inv$start[bb],(ypos + (bb/36)),sub_inv$end[bb],(ypos + (bb/36)),length = 0.015,lwd = 1, angle = 90, code = 3, col = col_inv[bb])
				}
			} 
				
			arrows(sub2$low[k],ypos,sub2$high[k],ypos,length = 0.03,lwd = 2, angle = 90, code = 3, col = sub2$the_dirs2_col[k])
						
							
		} #end for
		}	#end if
	}
}


dev.off()