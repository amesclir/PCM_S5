library(geiger)
## read data matrix
sqData<-read.csv("squamate-data.csv",row.names=1)
## print dimensions of our data frame
dim(sqData)
## read phylogenetic tree
sqTree<-read.nexus("squamate.tre")
print(sqTree,printlen=2)
## plot our tree
plotTree(sqTree,type="fan",lwd=1,fsize=0.3,ftype="i")
## check name matching
chk<-name.check(sqTree,sqData)
summary(chk)
## drop tips of tree that are missing from data matrix
sqTree.pruned<-drop.tip(sqTree,chk$tree_not_data)
## drop rows of matrix that are missing from tree
sqData.pruned<-sqData[!(rownames(sqData)%in%chk$data_not_tree),,drop=FALSE]
## extract discrete trait
toes<-setNames(as.factor(sqData.pruned[,"rear.toes"]),rownames(sqData.pruned))
head(toes)

## create design matrix for bi-directional
## ordered model 2 parameters

ordered.model2<-matrix(c(0,1,0,0,0,0,2,0,1,0,0,0,0,2,0,1,0,0,0,0,2,0,1,0,0,0,0,2,0,1,0,0,0,0,2,0),6,6,byrow=TRUE,dimnames=list(0:5,0:5))
ordered.model2

## create design matrix for directional ordered
## model one parameter
directional.model1<-matrix(c(0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0),6,6,byrow=TRUE,dimnames=list(0:5,0:5))
directional.model1

## fit bi-directional ordered model
fitOrdered2<-fitDiscrete(sqTree.pruned,toes,model=ordered.model2,surpressWarnings=TRUE)
print(fitOrdered2,digits=3)

## fit directional (loss only) ordered model
fitDirectional1<-fitDiscrete(sqTree.pruned,toes,model=directional.model1,surpressWarnings=TRUE)
print(fitDirectional1,digits=3)

## split plot area into two panels
par(mfrow=c(1,2))
## plot ordered and directional models
plot(fitOrdered2,show.zeros=FALSE,signif=5,mar=c(0.1,1.1,0.1,0.1))
mtext("(a)",line=-2,adj=0,cex=1.5)
plot(fitDirectional1,show.zeros=FALSE,signif=5,mar=c(0.1,1.1,0.1,0.1))
mtext("(b)",line=-2,adj=0,cex=1.5)


## accumulate AIC scores of all five models into
## a vector
aic<-setNames(c(AIC(fitER),AIC(fitDirectional),AIC(fitOrdered),AIC(fitSYM),AIC(fitARD),AIC(fitDirectional1),AIC(fitOrdered2)),c("ER","Directional","Ordered","SYM","ARD", "Directional1", "Ordered2"))
aic
aic.w(aic)
