
###### External shape data analysis ##########
##Library packages
library(geomorph)
library(shapes)
library(svd)
library(scatterplot3d)
library(rgl)
library(MASS)
library(ape)
library(vegan)
library(ggplot2)
library(ggfortify)
library(abind)
library(car)


###################################################################
################# External whole body ############################
###Reading data and general data preparation

dat<-readland.tps("Whole_f1.tps", negNA = TRUE, specID = c("imageID"))

#### Estimate missing landmarks
dat<-estimate.missing(dat,method = "TPS")


#### Flipping of individuals where it is needed.
PC4200<-dat[,,40]
PC4200r<-rotate.coords(PC4200, "flipX")
dat[,,40]<-PC4200r

PC42261<-dat[,,100]
PC42261r<-rotate.coords(PC42261, "flipX")
dat[,,100]<-PC42261r


#### Specifying which landmarks are sliding. LM =6,18,20,24,25,26,27
sliding<-matrix(c(5,6,7,17,18,19,18,19,20,26,24,25,24,25,26,25,26,27,26,27,24),ncol=3,byrow=TRUE)

####Performing Generlized Procrustes analysis
dat.lands<-gpagen(dat, curves = sliding)


### Adding character information
datCh<-read.csv("external_ch.csv", stringsAsFactors = TRUE, header = TRUE, sep = ",")

gdf.dat<-geomorph.data.frame(dat.lands, morph = datCh$morph, 
                             csize = dat.lands$Csize,logcsize = log(dat.lands$Csize),
                             ind = dimnames(dat.lands$coords)[[3]], length = datCh$length_cm, 
                             weight = datCh$weight_g, id = datCh$id_me, sex = datCh$sex, 
                             age = datCh$age)


#General examination of data.
datCh$logCsize<-gdf.dat$logcsize
ggplot(datCh, aes(x=log(length_cm), y=logCsize)) +geom_point()+ geom_smooth(method="lm", se=T) + xlab("Length")

plotAllSpecimens(dat.lands$coords, mean = TRUE)

plotOutliers(dat.lands$coords,inspect.outliers = TRUE)
outliers <- plotOutliers(gdf.dat$coords, groups = gdf.dat$morph, inspect.outliers = TRUE)



## Testing for Morph effect.
fit.dat0<-procD.lm(coords ~ logcsize, data = gdf.dat,iter = 999, RRPP=TRUE)
fit.dat1<-procD.lm(coords ~ logcsize + morph, data = gdf.dat,iter = 999, RRPP=TRUE)
fit.dat2<-procD.lm(coords ~ logcsize * morph, data = gdf.dat,iter = 999, RRPP=TRUE)

anova(fit.dat0, fit.dat1)
anova(fit.dat1, fit.dat2)

summary(fit.dat1)
summary(fit.dat2)

###Pairwise analysis
gp.dat<-interaction(gdf.dat$morph)

PW2<-pairwise(fit.dat1, fit.dat0,groups = gp.dat, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

## Testing allometry further
#Homogeneity of Slopes test, parallel (fit.den1) or non-parallel *includes covariate
gp.dat<-interaction(gdf.dat$morph)

PW<-pairwise(fit.dat2, fit.dat1,groups = gp.dat, covariate = gdf.dat$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
#Significant differences in vector angles (where the shape trajectories are going, convergence or divergence) *included covariate
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

### Plotting
palettething=c("#00BA38", #LB
               "#8B008B", #PI
               "#F8766D", #PL
               "#00BFC4") #SB
test<-plotAllometry(fit.dat2,gdf.dat$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.dat$morph)],xlab="Predictor")


spreds<-shape.predictor(gdf.dat$coords, test$PredLine, Intercept = TRUE, predmin = min(test$PredLine),predmax = max(test$PredLine))
M <- mshape(gdf.dat$coords)
plotRefToTarget(M, spreds$predmin, mag=1) #left
plotRefToTarget(M, spreds$predmax, mag=1) #right


## Plotting by morph (with size effects.)
palettec<-c("#00BA38","magenta4","#F8766D","#00BFC4")
gp.dat<-gdf.dat$morph
col.gp<-palettec
names(col.gp)<-levels(gp.dat)
col.gp.dat<-col.gp[match(gp.dat,names(col.gp))]
TS.dat<-gm.prcomp(gdf.dat$coords)
plot(TS.dat,pch=21, bg=col.gp.dat)

scatterplot(TS.dat$x[,1],TS.dat$x[,2],groups = gp.dat, grid = FALSE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.1,0.09),ylim=c(-0.05,0.05),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TS.dat$x[,2],TS.dat$x[,3],groups = gp.dat, grid = FALSE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.05,0.06),ylim=c(-0.05,0.05),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(gdf.dat$coords)


### Plot Warps
##### What are shape differences for each PC.
ref.dat<-mshape(gdf.dat$coords)
PC <- TS.dat$x[,1]
preds <- shape.predictor(gdf.dat$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.dat, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.dat, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
#### Run also for other PC.


## Allometry free shape
## Plotting by morph (without size effects.) 
### Can only be used if there is no difference in allometry between the morphs.
dat.Anova<-procD.lm(coords~logcsize, data = gdf.dat,iter=499,RRPP=TRUE)
shape.resid<-arrayspecs(dat.Anova$residuals,p=dim(dat.lands$coords)[1],k=dim(dat.lands$coords)[2])
adj.shape.dat<-shape.resid + array(dat.lands$consensus,dim(shape.resid))
TSadj.dat<-gm.prcomp(adj.shape.dat)
plot(TSadj.dat,pch=21,bg=col.gp.dat)

scatterplot(TSadj.dat$x[,1],TSadj.dat$x[,2],groups = gp.dat, grid = FALSE, smooth = FALSE, regLine = FALSE, col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.1,0.1),ylim=c(-0.05,0.05),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TSadj.dat$x[,2],TSadj.dat$x[,3],groups = gp.dat, grid = FALSE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.05,0.05),ylim=c(-0.04,0.04),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(adj.shape.dat)


### Plot Warps (for size adjusted)
##### What are shape differences for each PC.
ref.dat<-mshape(adj.shape.dat)
PC <- TSadj.dat$x[,1]
preds <- shape.predictor(adj.shape.dat, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.dat, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.dat, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
#### Run also for other PC.


## Allometry Regression
fit.size.null<-procD.lm(coords~logcsize,data=gdf.dat)
fit.size<-procD.lm(coords~logcsize*morph,groups=datCh$morph,data=gdf.dat)
A<-plotAllometry(fit.size,size=dat.lands$Csize,logsz=TRUE, method="RegScore",pch=19,col=col.gp.dat)
B<-plotAllometry(fit.size,size=dat.lands$Csize,logsz=TRUE, method="PredLine",pch=19,col=col.gp.dat)
C<-plotAllometry(fit.size,size=dat.lands$Csize,logsz=TRUE, method="size.shape",pch=19,col=col.gp.dat)
D<-plotAllometry(fit.size,size=dat.lands$Csize,logsz=TRUE, method="CAC",pch=19,col=col.gp.dat)


####################################################################
#####################################################################
############### Head morphology ################################

###Reading data and general data preparation

hau<-readland.tps("Head.tps", negNA = TRUE, specID = c("imageID"))

hau<-estimate.missing(hau,method = "TPS")


PC4200<-hau[,,40]
PC4200r<-rotate.coords(PC4200, "flipX")
hau[,,40]<-PC4200r

PC42261<-hau[,,100]
PC42261r<-rotate.coords(PC42261, "flipX")
hau[,,100]<-PC42261r

##### LM: 10, 11, 12, 13.
slidingH<-matrix(c(13,10,11,10,11,12,11,12,13,12,13,10),ncol=3,byrow=TRUE)

hau.lands<-gpagen(hau, curves = slidingH)


hauCh<-read.csv("external_ch.csv", stringsAsFactors = TRUE, header = TRUE, sep = ",")

gdf.hau<-geomorph.data.frame(hau.lands, morph = hauCh$morph, 
                             csize = hau.lands$Csize,
                             logcsize = log(hau.lands$Csize),
                             ind = dimnames(hau.lands$coords)[[3]],
                             length = hauCh$length_cm,weight = hauCh$weight_g, 
                             id = hauCh$id_me, sex = hauCh$sex, age = hauCh$age)


hauCh$logCsize<-gdf.hau$logcsize
ggplot(hauCh, aes(x=length_cm, y=logCsize)) +geom_point()+ geom_smooth(method="lm", se=T) + xlab("Length")

plotAllSpecimens(hau.lands$coords, mean = TRUE)
plotOutliers(hau.lands$coords,inspect.outliers = TRUE)
outliers <- plotOutliers(gdf.hau$coords, groups = gdf.hau$morph, inspect.outliers = TRUE)

########## Testing morph effects ###########################
fit.hau0<-procD.lm(coords ~ logcsize, data = gdf.hau,iter = 999, RRPP=TRUE)
fit.hau1<-procD.lm(coords ~ logcsize + morph, data = gdf.hau,iter = 999, RRPP=TRUE)
fit.hau2<-procD.lm(coords ~ logcsize * morph, data = gdf.hau,iter = 999, RRPP=TRUE)

anova(fit.hau0, fit.hau1)
anova(fit.hau1, fit.hau2)

summary(fit.hau1)
summary(fit.hau2)


gp.hau<-interaction(gdf.hau$morph)

PW2<-pairwise(fit.hau1, fit.hau0,groups = gp.hau, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means

## Testing allometry further

gp.hau<-interaction(gdf.hau$morph)

#Homogeneity of Slopes test, parallel (fit.den1) or non-parallel *includes covariate
PW<-pairwise(fit.hau2, fit.hau1,groups = gp.hau, covariate = gdf.hau$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

#Significant differences in vector angles (where the shape trajectories are going, convergence or divergence) *included covariate
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

palettething=c("#00BA38", #LB
               "#8B008B", #PI
               "#F8766D", #PL
               "#00BFC4") #SB
plotAllometry(fit.hau2,gdf.hau$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.hau$morph)], xlab="Predictor") 

spreds<-shape.predictor(gdf.hau$coords, test$PredLine, Intercept = TRUE, predmin = min(test$PredLine),predmax = max(test$PredLine))
M <- mshape(gdf.hau$coords)
plotRefToTarget(M, spreds$predmin, mag=1) #left
plotRefToTarget(M, spreds$predmax, mag=1) #right


## Plotting by morph (with size effects.)
palettec<-c("#00BA38","magenta4","#F8766D","#00BFC4")
gp.hau<-gdf.hau$morph
col.gp<-palettec
names(col.gp)<-levels(gp.hau)
col.gp.hau<-col.gp[match(gp.hau,names(col.gp))]
TS.hau<-gm.prcomp(gdf.hau$coords)
plot(TS.hau,pch=21, bg=col.gp.hau)

scatterplot(TS.hau$x[,1],TS.hau$x[,2],groups = gp.hau, grid = FALSE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.12,0.12),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TS.hau$x[,2],TS.hau$x[,3],groups = gp.hau, grid = FALSE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.1,0.09),ylim=c(-0.1,0.1),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(gdf.hau$coords)

### Plot Warps
ref.dat<-mshape(gdf.hau$coords)
PC <- TS.hau$x[,1]
preds <- shape.predictor(gdf.hau$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.dat, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.dat, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry free shape
## Plotting by morph (without size effects.) 
hau.Anova<-procD.lm(coords~logcsize, data = gdf.hau,iter=499,RRPP=TRUE)
shape.resid<-arrayspecs(hau.Anova$residuals,p=dim(hau.lands$coords)[1],k=dim(hau.lands$coords)[2])
adj.shape.hau<-shape.resid + array(hau.lands$consensus,dim(shape.resid))
TSadj.hau<-gm.prcomp(adj.shape.hau)
plot(TSadj.hau,pch=21,bg=col.gp.dat)

scatterplot(TSadj.hau$x[,1],TSadj.hau$x[,2],groups = gp.hau, grid = FALSE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.1,0.1),ylim=c(-0.1,0.1),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TSadj.hau$x[,2],TSadj.hau$x[,3],groups = gp.hau, grid = FALSE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.1,0.1),ylim=c(-0.1,0.1),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(adj.shape.hau)

### Plot Warps (for size adjusted)
ref.hau<-mshape(adj.shape.hau)
PC <- TSadj.hau$x[,1]
preds <- shape.predictor(adj.shape.hau, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.hau, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.hau, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry Regression
fit.size.null<-procD.lm(coords~logcsize,data=gdf.hau)
fit.size<-procD.lm(coords~logcsize*morph,groups=hauCh$morph,data=gdf.hau)
A<-plotAllometry(fit.size,size=hau.lands$Csize,logsz=TRUE, method="RegScore",pch=19,col=col.gp.hau)
B<-plotAllometry(fit.size,size=hau.lands$Csize,logsz=TRUE, method="PredLine",pch=19,col=col.gp.hau)
C<-plotAllometry(fit.size,size=hau.lands$Csize,logsz=TRUE, method="size.shape",pch=19,col=col.gp.hau)
D<-plotAllometry(fit.size,size=hau.lands$Csize,logsz=TRUE, method="CAC",pch=19,col=col.gp.hau)

