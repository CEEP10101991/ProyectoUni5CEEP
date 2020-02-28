### load the necesary packages to perform the analysis

library(SNPRelate)
library(ape)
library(ggplot2)


##### For using the SNPRelate pachage its necesary to generate a ´gds´ file
## Generate gds file from plink format files
snpgdsBED2GDS("../data/plinkWolves.bed", 
              "../data/plinkWolves.fam", 
              "../data/plinkWolves.bim", 
              out.gdsfn="../data/plinkWolves.gds")
# Creat an object that reads the previously generated gds file
genofile <- snpgdsOpen("../data/plinkWolves.gds")

# Check snp.ids
head(read.gdsn(index.gdsn(genofile, "snp.id")))

# Check sample.ids
head(read.gdsn(index.gdsn(genofile, "sample.id")))

# obtain a lis of the gdsn samples
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample.id

### Run the PCA 

pca <- snpgdsPCA(genofile, num.thread=2)

# Calculate the percentage variation of the first compnents
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

x<-round(pc.percent, 2)
sum(x[1:4])
sum(x[1:10])
sum(x[1:30])


# Set the results at datafarame format 
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

# Plot
ggplot(data = tab, aes(x=EV2, y=EV1)) + geom_point() +
  ylab(paste0("eigenvector 1 explaining ", round(pc.percent, 2)[1], "%")) +
  xlab(paste0("eigenvector 2 explaining ", round(pc.percent, 2)[2], "%"))

####Creating a Map, that allows at first instance to visualice the ubication of the samples and the Ecotype or Status for each One  

#Package requiered
library(maptools)
#maptools contain some basic maps that aer usefull for first visualizations
data(wrld_simpl)
plot(wrld_simpl)

#Lobos is an object obtained from the meta data of wolves, that contains the georeferenciation of the samples
lobos<-read.csv("../data/AllSamples_n117_wLatLong.csv",header=TRUE,sep=",")
puntos<-SpatialPoints(lobos[, c("Longitude","Latitude")])#To create an object that contains spatial points from the Longitude an Latitude columns at "lobos"
#to generate the map
plot(puntos, pch=17, col=lobos$Ecotype.or.Status,add=TRUE)
legend("right","botom",unique(lobos$Ecotype.or.Status),col=1:length(lobos$Ecotype.or.Status),pch=1,cex = 0.4) #the legend performed shows the  Ecotype or Status of each ploted sample
