setwd("/Users/JohnChod2130/Downloads/riverMetaG/ARGs_(ResFinder)/")

Files <- read.csv(file="Output_resfinder_megaHit.txt",sep="\t",header=TRUE)

#Identity and coverage cutoff as shown in resFinder article:
Files.filt <- Files[Files$X.COVERAGE>59.9 & Files$X.IDENTITY>79.9,]

#Need to add one filter pass for all samples so that we don't lose samples. 
#In other words, we need to make sure we retain samples with 0s across the board
unique(Files.filt$X.FILE)
#Note: samples 4 needs to be added back to table

#Upload temp sample
toRemove <- read.csv(file="toRemove.txt",sep="\t",header=TRUE)

#Bind dataframes
Files.filt <- rbind(Files.filt,toRemove)

#load antimicrobial class file
amc <- read.csv("antiMicrobialClass.csv",header=TRUE,sep=",")
#Match genes to antimicrobial class and add as a new column to Files
Files.filt$class <- amc$Antimicrobial.class[match(Files.filt$PRODUCT,amc$Gene)]

library(dplyr)

#Group genes by antibiotic class 
ABR.melt <- Files.filt %>%
  group_by(X.FILE) %>%
  count(class)
ABR.melt.df <- as.data.frame(ABR.melt)

#recast by class
cast_data = as.data.frame(dcast(ABR.melt, X.FILE ~ class, sum))
#reorder data frame
cast_data.o <- cast_data[order(as.character(cast_data$X.FILE)),]

#Convert File names to rownames and remove from dataframe and add new names
cast_data.f <- cast_data.o[,-1]

#load Mapping File
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
write.csv(final.df,"AMR_classNormalized.csv", row.names = TRUE)
