library(ggplot2)
library(vegan)

setwd("~/Estatistica_R/mel_teste")
mat <- read.csv("pcoa.csv", header=TRUE, sep=";",dec=",")
head(mat)
colnames(mat)
rownames(mat) <- mat$Nome
mat <- mat[, -1]
mat2 <- t(mat)
head(mat2)
colnames(mat2)
rownames(mat2)

PCOA <- capscale(mat2~1, distance="bray") ## Chose distance
summary(PCOA)
plot(PCOA)

#----------------------------------------------#
#----------- SCORES PCOA-------------#
scores <- scores(PCOA)


scores <- list(species=scores(PCOA,display=c("species")),
               sites=scores(PCOA,display=c("sites")))

scoresSites <- rownames(scores$sites)
scoresSitesGrp <- as.character()
for (i in 1:length(scoresSites)) 
  scoresSitesGrp[i] <- substr(as.character(scoresSites[i]),1,nchar(scoresSites[i])-1)
scoresSitesGrp <- as.factor(scoresSitesGrp)
#scoresSiteColors <- sample(colors(),length(levels(scoresSitesGrp)),replace=F)
scoresSiteColors <- sample_colors<-c("red","red","green","green","blue","blue","brown","brown")

scoresSiteSymbols<-c(15,22,18,23,17,24,16,21)

#------------PLOT-----------------------#

ordiplot(scores, type="none", main="PCoA")

orditorp(scores,display="species", air=0.01)


points(scores$sites,
       col=scoresSiteColors[scoresSitesGrp], label=scoresSiteSymbols,
       pch=scoresSiteSymbols[scoresSitesGrp], cex=1)


legend("topleft",legend=levels(scoresSitesGrp),pch=scoresSiteSymbols,col=scoresSiteColors, cex=0.8)



#------------------------------------------------------------------------#
#----------SCORES NMDS-------------------------#
nmds <- metaMDS(mat2, k=2, trymax=100) 
stressplot(nmds)
plot(nmds)
nmds$stress   #k=2, stress=0.04
nmds$species
nmds$points

scores <- list(species=scores(nmds,display=c("species")),
               sites=scores(nmds,display=c("sites")))

scoresSites <- rownames(scores$sites)
scoresSitesGrp <- as.character()
for (i in 1:length(scoresSites)) 
  scoresSitesGrp[i] <- substr(as.character(scoresSites[i]),1,nchar(scoresSites[i])-1)
scoresSitesGrp <- as.factor(scoresSitesGrp)
#scoresSiteColors <- sample(colors(),length(levels(scoresSitesGrp)),replace=F)
scoresSiteColors <- sample_colors<-c("blue","red","blue","red","blue","red","blue","red")

scoresSiteSymbols<-c(15,15,18,18,17,17,16,16)

#------------PLOT-----------------------#
pdf("nmds2.pdf") 

ordiplot(scores,type="n", main="nMDS")

orditorp(scores,display="species", air=0.01)
#orditorp(scores,display="species", pch=3, labels = "F",air=0.01)

points(scores$sites,
       col=scoresSiteColors[scoresSitesGrp], 
       pch=scoresSiteSymbols[scoresSitesGrp], cex=1)


#Teste ordicluster
ordicluster(nmds,hclust(vegdist(mat2,"bray")))

legend("bottomright",legend=levels(scoresSitesGrp),pch=scoresSiteSymbols,col=scoresSiteColors, cex=0.8)

dev.off()


# Define a group variable (first 12 samples belong to group 1, last 12 samples to group 2)
group = c(rep("Group1", 4), rep("Group2", 4))

# Create a vector of color values with same length as the vector of group values
colors = c(rep("red", 4), rep("blue", 4))

# Plot convex hulls with colors based on the group identity
ordiplot(nmds, type = "n")
for(i in unique(group)) {
  ordihull(nmds$points[grep(i, group),], draw="polygon",
           groups = group[group == i],col = colors[grep(i,group)],label=F) } 

orditorp(nmds, display = "species", col = "red", air = 0.01)
orditorp(nmds, display = "sites", col = c(rep("red",4), rep("blue", 4)), air = 0.01, cex = 1.25)


data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- scoresSitesGrp  #  add the grp variable created earlier
head(data.scores) 
species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)

ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.3) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=site, shape=site), size=5) + scale_shape_manual(values=c(15,15,18,18,17,17,16,16)) +# add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=3,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("A" = "red", "B" = "blue", "C"="green", "D"= "yellow", "E"="black", "F"="navy", "G"="orange", "H"="maroon")) +
  coord_equal() +
  theme_bw()