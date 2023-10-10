setwd("/Users/JohnChod2130/Downloads/riverMetaG/EffluxPumpEncodingGenes_CARD/")

Files <- read.csv(file="Output_card_megaHit.txt",sep="\t",header=TRUE)

#Load list of efflux genes
efflux <- read.csv("Efflux_pumps.csv",header=TRUE,sep=",")

#Filter for efflux genes only since the CARD output detects more than just efflux genes
Files.e <- Files[which(Files$GENE %in% efflux$Gene),]

#Identity and coverage cutoff as shown in resFinder article:
Files.filt <- Files.e[Files.e$X.COVERAGE>59.9 & Files.e$X.IDENTITY>79.9,]

#Need to add one filter pass for all samples so that we don't lose samples. 
#In other words, we need to make sure we retain samples with 0s across the board
unique(Files.filt$X.FILE)

#Upload temp sample
toRemove <- read.csv(file="toRemove.txt",sep="\t",header=TRUE)

#Bind dataframes
Files.filt <- rbind(Files.filt,toRemove)

#Match data to mapping file so that we can obtain overall class of efflux
Files.filt$System <- efflux$System[match(Files.filt$GENE,efflux$Gene)]
nrow(Files.filt)


library(dplyr)

#Group by efflux system

efflux.melt <- Files.filt %>%
  group_by(X.FILE) %>%
  count(System)
efflux.melt.df <- as.data.frame(efflux.melt)

#Recast by efflux system
cast_data = as.data.frame(dcast(efflux.melt, X.FILE ~ System, sum))

#Convert File names to rownames and remove from dataframe and add new names
cast_data.f <- cast_data[,-1]

#load Mapping
mapping <- read.csv(file="map.csv",sep=",",header=TRUE)

#Make sure to re-confirm that order in mapping file is sequential, like the order in the dataframe 
#We are replacing names. Double check

mapping$name
cast_data$X.FILE

#Replace names
rownames(cast_data.f) <- mapping[,2]
cast_data.f

#Normalize by genome equivalents
vec <- mapping_o$Genome_equivalents
cast_data.f.norm <- sweep(cast_data.f,1,vec,FUN="/")

#Remove "remove" column
final.df <- cast_data.f.norm[ , !names(cast_data.f.norm) %in% 
                           c(" remove")]

#Save dataframe for downstream PCoA and ordiR2step analysis
write.csv(final.df,"Efflux_classNormalized.csv", row.names = TRUE)
