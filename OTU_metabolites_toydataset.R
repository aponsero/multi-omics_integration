library(mixOmics)
library(igraph)
library(tidyverse)
library(vegan)
library(factoextra)

#### function to perform pre-filtering
low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent=0.01 # cutoff chosen
){
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# Function to perform TSS
TSS.divide = function(x){
  x/sum(x)
}

#### loading data
data_tax_raw<- read.csv(file="/Users/aponsero/Documents/UA_POSTDOC/projects/Malak_test_dataset/data_test/taxon_data_OTU_fakeraw.csv", header=TRUE, sep=",")
data_metab<- read.csv(file="/Users/aponsero/Documents/UA_POSTDOC/projects/Malak_test_dataset/data_test/metabolite_data.csv", header=TRUE, sep=",")
### factors
depth <- factor(c("0-20","0-20","0-20","0-20","90-110","90-110","90-110","90-110"))
time <-factor(c("TI","TF","TI","TF","TI","TF","TI","TF"))
treatment <- factor(c("Control","Control","Glucose","Glucose","Control","Control","Glucose","Glucose"))
names <- c("0-20-TI","0-20-TF","0-20-U-TI","0-20-U-TF","90-110-TI","90-110-TF","90-110-U-TI","90-110-U-TF")
samples <- c("s1","s1","s2","s2","s4","s4","s5","s5")

#### OTU count preparation
data_tax_raw = data_tax_raw+1
head(data_tax_raw)
sum(which(data_tax_raw == 0))
# remove low OTU
result.filter <- low.count.removal(data_tax_raw, percent=0.01)
data.filter <- result.filter$data.filter
length(result.filter$keep.otu)
head(data.filter)
# TSS
data.TSS <- t(apply(data.filter, 1, TSS.divide))
sum(which(data.TSS == 0))
# logratio
summary(data.TSS)
data_tax_log_CLR <- logratio.transfo(data.TSS, logratio = "CLR")
data_tax_log_ILR <- logratio.transfo(data.TSS, logratio = "ILR")

summary(data_tax_log_CLR)
summary(data_tax_log_ILR)
sum(which(data_tax_log_ILR == "-Inf"))
boxplot(log(t(data.TSS))+1)

#### PCA analysis OTU count
tune.pca(data_tax_log_CLR, ncomp = 8, center = FALSE, scale = FALSE)

pca.taxon<- pca(data_tax_log_CLR, ncomp = 3, center = FALSE, scale = FALSE)
pca.taxon
plotVar(pca.taxon,var.names = NULL)
plotIndiv(pca.taxon, comp = c(1, 2), ind.names=treatment, group=depth, legend=TRUE, title = '16s OTU, PCA comp 1 - 2')
plotIndiv(pca.taxon, comp = c(1, 2), ind.names=names, group=depth, legend=TRUE, title = '16s OTU, PCA comp 1 - 2')
plotIndiv(pca.taxon, comp = c(1, 2), ind.names=names, group=time, legend=TRUE, title = '16s OTU, PCA comp 1 - 2')

#sparse PCA
spca.taxon<- spca(data_tax_log_CLR, ncomp = 3, center = FALSE, scale = FALSE)
spca.taxon
plotVar(spca.taxon,var.names = NULL)
plotIndiv(spca.taxon, comp = c(1, 2), ind.names=treatment, group=depth, legend=TRUE, title = '16s OTU, PCA comp 1 - 2')
plotIndiv(spca.taxon, comp = c(1, 2), ind.names=names, group=depth, legend=TRUE, title = '16s OTU, PCA comp 1 - 2')
plotIndiv(spca.taxon, comp = c(1, 2), ind.names=names, group=time, legend=TRUE, title = '16s OTU, PCA comp 1 - 2')


#### PCA analysis metabolite
tune.pca(data_metab, ncomp = 8, center = FALSE, scale = FALSE)
tune.pca(data_metab, ncomp = 8, center = TRUE, scale = TRUE)

pca.metab<- pca(data_metab, ncomp = 3, center = TRUE, scale = TRUE)
pca.metab
plotVar(pca.metab,var.names = NULL)
plotIndiv(pca.metab, comp = c(1, 2), ind.names=treatment, group=depth, legend=TRUE, title = 'metabolites, PCA comp 1 - 2')
plotIndiv(pca.metab, comp = c(1, 2), ind.names=names, group=depth, legend=TRUE, title = 'metabolites, PCA comp 1 - 2')
plotIndiv(pca.metab, comp = c(1, 2), ind.names=names, group=time, legend=TRUE, title = 'metabolites, PCA comp 1 - 2')

#############procrustes analysis#################
X <- data_tax_log_CLR
Y <- data_metab
row.names(X) <- c("0-20-TI","0-20-TF","0-20-U -TI","0-20-U-TF","90-110-TI","90-110-TF","90-110-U-TI","90-110-U-TF")
row.names(Y) <- c("0-20-TI","0-20-TF","0-20-U -TI","0-20-U-TF","90-110-TI","90-110-TF","90-110-U-TI","90-110-U-TF")
head(cbind(rownames(X), rownames(Y))) 

#run PCA
X.pca <- prcomp(X, scale=FALSE, center=FALSE)
fviz_eig(X.pca)
fviz_pca_ind(X.pca,repel = TRUE, center=TRUE)

Y.pca <- prcomp(Y, scale=TRUE, center=TRUE)
fviz_eig(Y.pca)
fviz_pca_ind(Y.pca,repel = TRUE)

#get the Procrustes analysis done
analysis.pro <- procrustes(X.pca, Y.pca, symmetric = TRUE)
analysis.pro
summary(analysis.pro)
plot(analysis.pro, kind = "1")
plot(analysis.pro, kind = "2")

analysis2.pro <- protest(X.pca, Y.pca, symmetric = TRUE)
analysis2.pro
summary(analysis2.pro)
plot(analysis2.pro, kind = "1")
plot(analysis2.pro, kind = "2")


##### PLS analysis regression ######
toy_data.pls <- pls(data_tax_log_CLR, data_metab, ncomp = 6, mode = "regression")
tune.pls <- perf(toy_data.pls, validation = "loo", progressBar = FALSE, nrepeat = 50)
plot(tune.pls$Q2.total)
abline(h = 0.0975)
tune.pls$Q2.total
plotIndiv(toy_data.pls, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'Y-variate',
         title = 'Data_toy, PLS comp 1 - 2, Y-space')

plotIndiv(toy_data.pls, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'X-variate',
          title = 'Data_toy, PLS comp 1 - 2, X-space')

plotIndiv(toy_data.pls, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'XY-variate', 
          title = 'Data_toy, PLS comp 1 - 2, XY-space')

plotVar(toy_data.pls, comp = 1:2, title = 'Data_toy PLS_reg')

##### PLS analysis canonical ######
toy_data.pls_c <- pls(data_tax_log_CLR, data_metab, ncomp = 6, mode = "canonical")

plotIndiv(toy_data.pls_c, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'Y-variate',
          title = 'Data_toy, PLS comp 1 - 2, Y-space')

plotIndiv(toy_data.pls_c, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'X-variate',
          title = 'Data_toy, PLS comp 1 - 2, X-space')

plotIndiv(toy_data.pls_c, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'XY-variate', 
          title = 'Data_toy, PLS comp 1 - 2, XY-space')

plotVar(toy_data.pls_c, comp = 1:2, title = 'Data_toy PLS_reg')

#### sPLS analysis canonical ####
toy_data.spls_c <- spls(data_tax_log_CLR, data_metab,keepX=c(50,50,50),
                        keepY=c(10,10,10),
                        ncomp = 3, mode = "canonical")

plotIndiv(toy_data.spls_c, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'Y-variate',
          title = 'Data_toy, sPLS comp 1 - 2, Y-space')

plotIndiv(toy_data.spls_c, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'X-variate',
          title = 'Data_toy, sPLS comp 1 - 2, X-space')

plotIndiv(toy_data.spls_c, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'XY-variate', 
          title = 'Data_toy, sPLS comp 1 - 2, XY-space')

plotVar(toy_data.spls_c, comp = 1:2, title = 'Data_toy sPLS_can')

#### impact of the glucose treatment ?
### PLS-DA on OTU
toy_data.plsda <- plsda(data_tax_log_CLR, treatment)
toy_data.plsda
plotIndiv(toy_data.plsda, group = treatment, legend = TRUE,
          ind.names = names)
plotVar(toy_data.plsda, title = 'Data_toy PLS-DA_reg')

### PLS-DA on metabolites
toy_data_metab.plsda <- plsda(data_metab, treatment)
toy_data_metab.plsda
plotIndiv(toy_data_metab.plsda, group = treatment, legend = TRUE,
          ind.names = names)
plotVar(toy_data_metab.plsda, title = 'Data_toy PLS-DA_reg')
plotLoadings(toy_data_metab.plsda, comp = 1, contrib = 'max')
# depth_max -> X1.Butanol, Butyrate, IsoButyrate, Isopropanol
plotLoadings(toy_data_metab.plsda, comp = 2, contrib = 'max')
# glucose -> Methanol, Propyla.alcohol


### Run multiblock spls-da
# set number of variables to select, per component and per data set
toyData <- list(taxons = data_tax_log_CLR, metabolites = data_metab)

list.keepX <- list(taxons = rep(20, 2), metabolites = rep(10,2))

toyData.block.splsda <- block.splsda(X = toyData, Y = treatment,
                                  ncomp = 2, keepX = list.keepX)
toyData.block.splsda

# Individuals plot
plotIndiv(toyData.block.splsda, ind.names = depth, legend=TRUE)

# Variables plot
# discrimination treatments on comp2
# depth discriminated on comp1
plotLoadings(toyData.block.splsda, comp = 1, contrib = 'max')
plotLoadings(toyData.block.splsda, comp = 2, contrib = 'max')
plotVar(toyData.block.splsda, legend = TRUE)

circosPlot(toyData.block.splsda, cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

network(toyData.block.splsda, blocks = c(1,2),
        color.node = c('darkorchid', 'brown1'), cutoff = 0.7)

cim(toyData.block.splsda)

###########################################
##### Multilevel analysis
design <- data.frame(sample = samples)

##PCA multilevel OTU
pca.multilevel_OTU <- pca(data_tax_log_CLR, ncomp = 3, scale = FALSE, center = FALSE, 
                      multilevel = design)
pca.multilevel_OTU
plotIndiv(pca.multilevel_OTU, ind.names = names,legend=TRUE, 
          group = depth,
          title = 'PCA multilevel OTU, comp 1 - 2')

##PCA multilevel metabolome
pca.multilevel_meta <- pca(data_metab, ncomp = 3, scale = FALSE, center = FALSE, 
                      multilevel = design)
pca.multilevel_meta
plotIndiv(pca.multilevel_meta, ind.names = names,legend=TRUE, 
          group = depth,
          title = 'PCA multilevel, comp 1 - 2')
?plotIndiv
plotIndiv(pca.multilevel_meta, comp=c(1:3) ,ind.names = names,legend=TRUE, 
          group = depth,
          title = 'PCA multilevel, comp 1 - 3', style='3d')

# multilevel PLS
pls.multilevel <- pls(data_tax_log_CLR, data_metab, multilevel = design, ncomp = 3)
plotIndiv(pls.multilevel, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'Y-variate',
          title = 'Data_toy, PLS multilevel comp 1 - 2, Y-space')

plotIndiv(pls.multilevel, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'X-variate',
          title = 'Data_toy, PLS multilevel comp 1 - 2, X-space')

plotIndiv(pls.multilevel, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'XY-variate', 
          title = 'Data_toy, PLS multilevel comp 1 - 2, XY-space')

plotVar(pls.multilevel, comp = 1:2, title = 'Data_toy PLS multilevel')

# multilevel sPLS
spls.multilevel <- spls(data_tax_log_CLR, data_metab, multilevel = design,
                        keepX = c(20,20,20), keepY= c(10,10,10), ncomp = 3)
plotIndiv(spls.multilevel, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'Y-variate',
          title = 'Data_toy, sPLS multilevel comp 1 - 2, Y-space')

plotIndiv(spls.multilevel, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'X-variate',
          title = 'Data_toy, sPLS multilevel comp 1 - 2, X-space')

plotIndiv(spls.multilevel, comp = 1:2, ind.names=names, group=depth, legend=TRUE, rep.space= 'XY-variate', 
          title = 'Data_toy, sPLS multilevel comp 1 - 2, XY-space')

plotVar(spls.multilevel, comp = 1:2, title = 'Data_toy sPLS multilevel')


#multilevel sPLS-DA on metabolites
multi_data_metab.plsda <- plsda(data_metab, treatment, ncomp=3)
multi_data_metab.plsda
plotIndiv(multi_data_metab.plsda, group = treatment, legend = TRUE,
          ind.names = names, comp=c(1:3), style="3d")
plotVar(multi_data_metab.plsda, title = 'Data_toy PLS-DA_reg')
plotLoadings(multi_data_metab.plsda, comp = 1, contrib = 'max')
plotLoadings(multi_data_metab.plsda, comp = 2, contrib = 'max')
