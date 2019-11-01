library(poppr);library(adegenet);library(pegas);library(hierfstat);library(ggplot2);library("magrittr");library("vegan");library(PopGenReport);library(HWxtest);library(genepop)
setwd("~/Dropbox/Marco_Crotti/Evolutionary genomics of whitefish/Microsatellites")
setwd("C:/Users/mario/Dropbox/Marco_Crotti/Evolutionary genomics of whitefish/Microsatellites")

whitefish <- read.genalex("uk_whitefish_micro.csv")
splitStrata(whitefish) <- ~Lake/Country
setPop(whitefish) <- ~Lake

whitefish_genind <- genclone2genind(whitefish)
whitefish.hfstat <- genind2hierfstat(whitefish)
whitefish_genepop <- genind2genpop(whitefish_genind)
whitefish.loci <- genind2loci(whitefish_genind)

whitefish

### data investigation
(whitefishlt <- locus_table(whitefish))
info_table(whitefish, type = "missing", plot = TRUE)

poppr(whitefish)


### monomorphic loci

locus_table(whitefish)


locus_table(whitefish, pop = "lom")
locus_table(whitefish, pop = "eck")
locus_table(whitefish, pop = "lte")
locus_table(whitefish, pop = "rta")
locus_table(whitefish, pop = "haw")
locus_table(whitefish, pop = "bwa")
locus_table(whitefish, pop = "uwa")

pop <- c("lom","eck","lte","rta","hwa","bwa","uwa")
nulls <- c(3,4,1,1,0,0,0)
names(nulls) <- pop
barplot(nulls, cex.axis = 2, cex.names = 2)


### HWE 

hwe_test <- hwx.test(whitefish_genind, B=1000000, histobins=T, statName = "LLR", showCurve = TRUE)

LLRwhitefish <- hwdf(hwe_test)
p.adjust(LLRwhitefish$`P-val(LLR)`,method = "bonferroni")



test_HW("whitefish.genepop.txt",outputFile = "hwe_test.txt",dememorization = 100000, batches = 1000, iterations = 10000)


### Linkage disequilibrium

test_LD("whitefish.genepop.txt",dememorization = 100000, batches = 1000, iterations = 10000)


#################
#################
#################  Population structure

## AMOVA

whitefish.amova <- poppr.amova(whitefish, ~Country/Lake)
whitefish.amova
set.seed(1999)
amova.test <- randtest(whitefish.amova, nrepet = 999)
amova.test
plot(amova.test)

## Fst
write.csv(whitefish.hfstat, "whitefish.hfstat.csv")
whitefish.hfstat <- read.csv("whitefish.hfstat.csv", header=TRUE)
fstwhitefish <- pairwise.WCfst(whitefish.hfstat[,-1],diploid = TRUE )
fstwhitefish
wc(whitefish.hfstat[,-1])$FST
varcomp.glob(data.frame(whitefish.hfstat$pop), whitefish.hfstat[,-c(1:2)])$F
boot.vc(whitefish.hfstat[,2],whitefish.hfstat[,-c(1:2)], nboot=10000)$ci
boot.ppfst(whitefish.hfstat[,-1])$ul
boot.ppfst(whitefish.hfstat[,-1])$ll

test.within(whitefish.hfstat[,-c(1:2)], within = whitefish.hfstat[,1], test.lev = whitefish.hfstat[,2])

### DAPC
whitefish_genind <- missingno(whitefish_genind, type = "genotype", cutoff = 0)

# no priors
(groups1 <- find.clusters(whitefish_genind, max.n.clust=10, n.pca = 300,
                         choose.n.clust = FALSE, criterion = "min"))
groups1
xval1 <- xvalDapc(tab(whitefish_genind, NA.method="mean"), groups1$grp, n.pca.max = 300,
                  result = "groupMean", center = TRUE, scale = FALSE, parallel = "multicore",
                  n.pca = NULL, n.rep = 100, xval.plot = TRUE)
xval1
dapc1 <- dapc(whitefish_genind, pop = groups1$grp, n.pca=20, n.da = 5)
scatter(dapc1, 1,2, posi.da = NULL, scree.pca = TRUE, posi.pca = NULL, solid = .9, bg = "white",col=cols1)
summary(dapc1)
assignplot(dapc1)
dapca1.eig.perc <- 100*dapc1$eig/sum(dapc1$eig)
dapca1.eig.perc[1:5]

compoplot(dapc1, posi="bottomright",
          txt.leg=paste("Cluster", 1:6), lab="",
          ncol=1, xlab="individuals", col=funky(6))

posterior <- data.frame(dapc1$posterior)
colnames(posterior) <- c("Group1","Group2","Group3","Group4","Group5","Group6")
cols1 <- c("Group1"="firebrick1","Group2"="#e3655b","Group3"="gold1","Group4"="#74BEE9","Group5"="#DC8C91","Group6"="#2372CD")
barplot(t(dapc1$posterior), col=cols,
        xlab="Lake", ylab="Ancestry", border=NA,cex.names = 0.1)
write.csv(dapc1$ind.coord,"dapc1_DFs.csv")


# priors
xval2 <- xvalDapc(tab(whitefish_genind, NA.method="mean"), whitefish_genind$pop, n.pca.max = 300,
                  result = "groupMean", center = TRUE, scale = FALSE, parallel = "multicore",
                  n.pca = NULL, n.rep = 100, xval.plot = TRUE)
xval2
dapc2 <- dapc(whitefish_genind, pop = whitefish_genind$pop, n.pca=40, n.da = 6)
scatter(dapc2, 1,2, posi.da = NULL, scree.pca = TRUE, posi.pca = NULL, solid = .9, bg = "white",col=cols2, cellipse = 0)
summary(dapc2)
assignplot(dapc2)
dapca2.eig.perc <- 100*dapc2$eig/sum(dapc2$eig)
dapca2.eig.perc[1:5]

contrib1 <- loadingplot(dapc1$var.contr, axis=1,
                       thres=.07, lab.jitter=1)
contrib2 <- loadingplot(dapc1$var.contr, axis=2,
                        thres=.07, lab.jitter=1)

col <- funky(length(levels(dapc1$grp)))

compoplot(dapc2, posi="bottomright",
          txt.leg=paste("Cluster", 1:6), lab="",
          ncol=1, xlab="individuals", col=funky(6))

posterior2 <- data.frame(dapc2$posterior)
colnames(posterior2) <- c("LLA","ECK","LTE","RTA","HWA","BWA","UWA")
cols2 <- c("LLA"="#2372CD","ECK"="#74BEE9","LTE"="gold1","RTA"="firebrick1","HWA"="#DC8C91","BWA"="#e3655b","UWA"="firebrick4")
barplot(t(posterior2), col = cols2,
        xlab="Lake", ylab="Ancestry", border=NA,cex.names = 0.1)

## PCA
X <- scaleGen(whitefish_genind, NA.method="mean")
class(X)
pca1 <- dudi.pca(X,cent=TRUE,scale=TRUE,scannf=TRUE,nf=3)
barplot(pca1$eig[1:5],main="PCA eigenvalues", col=heat.colors(50))
s.label(pca1$li,xax=1, yax=2)
s.class(pca1$li, fac=pop(whitefish_genind), col=cols)
eig.perc <- 100*pca1$eig/sum(pca1$eig)
eig.perc[1:4]

contrib1 <- loadingplot(pca1$co, axis=1,
                        thres=.07, lab.jitter=1)
contrib2 <- loadingplot(pca1$co, axis=2,
                        thres=.07, lab.jitter=1)

write.csv(pca1$li,"whitefish.pca2.csv")

whitefish.pca <- read.csv("whitefish.pca1.csv",header=TRUE)

whitefish.pca$Pop <- factor(whitefish.pca$Pop, levels = c("LLA", "ECK", "LTE","RTA","HAW","BWA","UWA"))
cols <- c("HAW"="#DC8C91","BWA"="#e3655b","ECK"="#74BEE9","LLA"="#2372CD","RTA"="firebrick1","LTE"="gold1","UWA"="firebrick4")
ggplot(whitefish.pca, aes(x=PC1, y=PC2, color=Pop)) + geom_point(size=3) + scale_color_manual(values=cols) + theme_bw() +
  theme(axis.title.y = element_text(size = 30),axis.title.x = element_text(size = 30),axis.text.x = element_text(size=21),axis.text.y = element_text(size=21)) +
  theme(legend.title = element_text(size=10)) + theme(legend.text = element_text(size=18)) + labs(x="PC1 12.4%",y="PC2 6.6%")
  
ggplot(whitefish.pca, aes(x=PC3, y=PC4, color=Pop)) + geom_point(size=3) + scale_color_manual(values=cols) + theme_bw() +
  theme(axis.title.y = element_text(size = 30),axis.title.x = element_text(size = 30),axis.text.x = element_text(size=21),axis.text.y = element_text(size=21)) +
  theme(legend.title = element_text(size=10)) + theme(legend.text = element_text(size=18)) + labs(x="PC1 12.4%",y="PC3 2.9%")


## NJ tree
library(phytools)
dist <- genet.dist(whitefish.hfstat) 
njtr <- nj(dist)
plot(njtr, "cladogram")
matrix1 <- as.matrix(whitefish.hfstat)

tree_out <- aboot(out_genind,strata = ~ Lake, distance = "edwards.dist",tree = "nj", sample = 1000, cex.names =2)

write.tree(tree_out,"whitefish_out.tre")


########### Isolation by distance

coords <- read.csv("lakes_coords.csv", header=TRUE)
lake_dist <- dist(cbind(coords$Long,coords$Lat))
gen_dist <- dist.genpop(whitefish_genepop, method = 2)

ibd <- mantel.randtest(gen_dist,lake_dist, nrepet = 10000)
ibd
plot(ibd)

plot(lake_dist, gen_dist); abline(lm(gen_dist~lake_dist), col="red",lty=2)


### genetic diversity ~ latitude correlation

lat.data <- read.csv("genetics_latitude.csv",header=TRUE)

HO <- cor.test(lat.data$Het,lat.data$Lat,method = "spearman")
plot(lat.data$Het ~ lat.data$Lat)
abline(lm(lat.data$Het ~ lat.data$Lat))

AR <- cor.test(lat.data$AR, lat.data$Lat, method = "spearman")
plot(lat.data$AR ~ lat.data$Lat)

par(mfrow=c(2,1))
plot(lat.data$Het ~ lat.data$Lat, ylab="Observed heterozygosity",
     xlab="Latitude",pch=19,cex=2,ylim=c(0.2,0.55))
abline(lm(lat.data$Het ~ lat.data$Lat))
plot(lat.data$AR ~ lat.data$Lat, ylab=" Allelic richness",
     xlab="Latitude",pch=19,cex=2,ylim=c(2,5.5))
abline(lm(lat.data$AR ~ lat.data$Lat))
