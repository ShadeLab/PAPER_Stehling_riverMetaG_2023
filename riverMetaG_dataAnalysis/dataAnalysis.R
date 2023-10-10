setwd("/Users/JohnChod2130/Downloads/riverMetaG")

amrN <- read.csv(file="AMR_classNormalized.csv",sep=",",header=TRUE, row.names =1, check.names = FALSE)
efxN <- read.csv(file="Efflux_classNormalized.csv",sep=",",header=TRUE, row.names =1, check.names = FALSE)
virN <- read.csv(file="Virulence_classNormalized.csv",sep=",",header=TRUE, row.names =1, check.names = FALSE)

#Combine datasets
genesC <- cbind(amrN,efxN,virN)

#load Mapping
mapping <- read.csv(file="map.csv",sep=",",header=TRUE)

#double check order in df matches order in mapping file
match(row.names(genesC),mapping$name)

#Load vegan
library(vegan)
dist.genesC <- vegdist(genesC, method="bray")
#dist.ABR <- vegdist(cast_data.f.norm, method="euclidean")

#Create groups
groups <- c(mapping$ID)
#Calculate betadisper
mod <- betadisper(dist.genesC, groups)
#Calculates spatial median, not center of mass

#Calculate variance explained on PCoA axes 1 & 2
PC1var <- mod$eig[1]/sum(mod$eig) #Checks out
PC2var <- mod$eig[2]/sum(mod$eig) #Checks out

#Extract scores
modScores <- scores(mod)

#Extract centroid
centroids <- as.data.frame(modScores$centroids)

#Add ID
centroids$Label <- mapping$ID
#Add water quality
centroids$WaterQ <- mapping$waterQ
#Add river type
centroids$type <- mapping$type
#Add time
centroids$sampling <- mapping$sampling

library(ggplot2)
#Change time from factor to numeric
#centroids$Time <- rep(c(12.5,25,30,35,40,45),4)

library(viridis)

PCA_NonPolarNeg <- ggplot(centroids, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(color = type, fill =type, shape = WaterQ, size=as.factor(sampling)))+
  #geom_point(aes(shape=factor(WaterQ),fill=factor(type)),colour="black",pch=21)+
  xlab(label = paste("PC1"," (", format(round(PC1var,digits = 3)*100,nsmall=1), "% var. explained)", sep = ""))+
  ylab(label = paste("PC2"," (", format(round(PC2var,digits = 3)*100,nsmall=1), "% var. explained)", sep = ""))+
  scale_x_continuous(breaks = pretty(centroids$PCoA1, n = 5)) +
  scale_y_continuous(breaks = pretty(centroids$PCoA2, n = 7)) +
  scale_shape_discrete(breaks=c('good', 'regular', 'poor')) +
  scale_shape_manual(values = c(22, 19, 24)) +   #19 is circle (poor), 22 is square (Good), 24 is triangle (regular)
  scale_fill_manual(values=c("#CC79A7", "#F0E442","#56B4E9")) +
  scale_color_manual(values=c("#CC79A7", "#F0E442","#56B4E9")) +
  scale_size_manual(values=c(2,3.5)) +
  theme_bw()+
  theme(legend.position="right",axis.title = element_text(size = 10),axis.text = element_text(size = 8),
        legend.title=element_text(size=8))



#Stil need to figure out a smart way to get actual paired samples for geom_seg. These are just in
#random order and do not represent pairs 
segA <- mapping[mapping$segments=="A",]
segA <- segA[order(segA$order),]

segB <- mapping[mapping$segments=="B",]
segB <- segB[order(segB$order),]


segment_data = data.frame(
  x = c(centroids$PCoA1[match(segA$ID,mapping$ID)]),
  xend = c(centroids$PCoA1[match(segB$ID,mapping$ID)]), 
  y = c(centroids$PCoA2[match(segA$ID,mapping$ID)]),
  yend = c(centroids$PCoA2[match(segB$ID,mapping$ID)])
)

#Generate PCoA
PCA_NonPolarNeg+
  geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend)) +
  labs(size = "Sampling") + 
  labs(color = "River") + 
  labs(fill = "River") +
  labs(shape = "Water quality") 

#ggsave("CA_NPN_Plot.eps",plot=PCA_NonPolarNegFinal,device="eps",width=4,height=4, units = "in", dpi=300)

#Run variation partitioning.
NPNvarPart <- varpart(dist.genesC,~ID,~waterQ,data=mapping)

adonis(dist.ABR ~ waterQ * type, data=mapping,permutations=999)

fit <- envfit(modScores, mapping[,13:30], perm=999)
fit 

AMR.dca <- decorana(cast_data.f.norm)
summary(AMR.dca)
plot(AMR.dca)

#create mapping file with env variables
envVar <- mapping[,13:30]
testRDA <- rda(cast_data.f.norm ~ ., envVar)
anova(testRDA)


capscaleTest <- capscale(genesC ~ E_coli + pH + BOD + Phosphor + Temperature
                         + Turbidity + Total_solids + Dissolved_Oxygen, mapping, dist="bray",sqrt.dist=TRUE)

capscale.empty <- capscale(genesC ~ 1, mapping, dist="bray",sqrt.dist=TRUE)


sel.osR2 <- ordiR2step(capscale.empty, scope = formula(capscaleTest), direction="forward", permutations=10000)
sel.osR2$anova


      