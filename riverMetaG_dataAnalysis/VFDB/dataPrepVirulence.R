setwd("/Users/JohnChod2130/Downloads/riverMetaG/VFDB/")

#Read in file
Files <- read.csv(file="Output_vfdb_megaHit.txt",sep="\t",header=TRUE)

#Identity and coverage cutoff as shown in resFinder article:
Files.filt <- Files[Files$X.COVERAGE>59.9 & Files$X.IDENTITY>79.9,]

#Need to add one filter pass for all samples so that we don't lose samples. 
#In other words, we need to make sure we retain samples with 0s across the board
unique(Files.filt$X.FILE)
#Sample 5 needs to be added back to table

toRemove <- read.csv(file="toRemove.txt",sep="\t",header=TRUE)

Files.filt <- rbind(Files.filt,toRemove)

#load vfdb class file
vfdb.class <- read.csv("VFDBs_mapFile.csv",header=TRUE,sep=",")
#Match data to mapping file so that we can obtain overall class of virulence
Files.filt$class <- vfdb.class$Virulence_factor[match(Files.filt$GENE,vfdb.class$Gene)]

library(dplyr)

#Group by class
vrl.melt <- Files.filt %>%
  group_by(X.FILE) %>%
  count(class)
vrl.melt.df <- as.data.frame(vrl.melt)


library(reshape2)
#recast by class
cast_data = as.data.frame(dcast(vrl.melt, X.FILE ~ class, sum))
#reorder data frame
cast_data.o <- cast_data[order(as.character(cast_data$X.FILE)),]

#Convert File names to rownames and remove from dataframe and add new names
cast_data.f <- cast_data.o[,-1]

#load Mapping
mapping <- read.csv(file="map.csv",sep=",",header=TRUE)

#Re-order mapping file based on name
mapping_o <- mapping[order(mapping$name),]

#Make sure to re-confirm that order in mapping file is sequential, like the order in the dataframe 
#We are replacing names. Double check

mapping_o$name
cast_data.o$X.FILE

#Replace names
rownames(cast_data.f) <- mapping_o[,2]
cast_data.f

#Normalize by genome equivalents
vec <- mapping_o$Genome_equivalents
cast_data.f.norm <- sweep(cast_data.f,1,vec,FUN="/")

#Remove "remove" column
final.df <- cast_data.f.norm[ , !names(cast_data.f.norm) %in% 
      c("remove")]

#Save dataframe for downstream PCoA and ordiR2step analysis
write.csv(final.df,"Virulence_classNormalized.csv", row.names = TRUE)
