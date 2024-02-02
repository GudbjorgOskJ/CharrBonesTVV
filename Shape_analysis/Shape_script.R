###### Internal shape data analysis ##########
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
library(patchwork)
library(emmeans)
library(lsmeans)


###################################################################
################# Dentary shape ############################
###Reading data and general data preparation

##Reading data
den<-readland.tps("Dentary_LM.tps", negNA = TRUE, specID = c("imageID"))

#### Estimate missing landmarks
den1<-estimate.missing(den,method = "TPS")

#### Flipping of individuals where it is needed.
Need_flipping_den<-den1[,,40:41]
Flipped_den<-rotate.coords(Need_flipping_den, "flipX")
den1[,,40:41]<-Flipped_den

den1<-rotate.coords(den1, type = c("rotateC"))
den1<-rotate.coords(den1, type = c("flipX"))

#### Specifying which landmarks are sliding. LM = 10, 14 and 15 for Dentary
sliding <- matrix(c(9, 10, 11, 13, 14, 15, 14, 15, 1), ncol = 3, byrow = TRUE)

####Performing Generlized Procrustes analysis
gpa.lands<-gpagen(den1, curves = sliding)


### Adding character information
denCh<- read.csv("Dentary_Ch.csv", stringsAsFactors = TRUE, header=TRUE)


gdf.den<-geomorph.data.frame(gpa.lands, morph = denCh$morph, 
                             csize = gpa.lands$Csize,logcsize = log(gpa.lands$Csize),
                             ind = dimnames(gpa.lands$coords)[[3]], length = denCh$length_cm,
                             weight = denCh$weight_g, id = denCh$id_me, sex = denCh$sex, 
                             age = denCh$age, loglength = log(denCh$length_cm))

gdf.den<-geomorph.data.frame(gpa.lands, morph = denCh$morph, 
                             csize = gpa.lands$Csize,logcsize = log(gpa.lands$Csize),
                             ind = dimnames(gpa.lands$coords)[[3]], length = denCh$length_cm,
                             weight = denCh$weight_g, id = denCh$id_me, sex = denCh$sex, 
                             age = denCh$age, loglength = log(denCh$length_cm), teethd = denCh$teeth_den,
                             teethm =denCh$teeth_max, teethpr = denCh$teeth_pre, teethpa = denCh$teeth_pal,
                             teethv = denCh$teeth_vo, teethg = denCh$teeth_glos)


#General examination of data.
denCh$logCsize<-gdf.den$logcsize
ggplot(denCh, aes(x=length_cm, y=logCsize)) +geom_point()+ geom_smooth(method="lm", se=T) + xlab("Length")

plotAllSpecimens(gpa.lands$coords, mean = TRUE)

plotOutliers(gpa.lands$coords,inspect.outliers = TRUE)
outliers <- plotOutliers(gdf.den$coords, groups = gdf.den$morph, inspect.outliers = TRUE) # shows outliers by group.

## Testing for Morph effect.
fit.den0<-procD.lm(coords ~ logcsize, data = gdf.den,iter = 999, RRPP=TRUE)
fit.den1<-procD.lm(coords ~ logcsize + morph, data = gdf.den,iter = 999, RRPP=TRUE)
fit.den2<-procD.lm(coords ~ logcsize * morph, data = gdf.den,iter = 999, RRPP=TRUE)

anova(fit.den0, fit.den1)
anova(fit.den1, fit.den2)

summary(fit.den1)
summary(fit.den2)

###Pairwise analysis
gp.den<-interaction(gdf.den$morph)

PW2<-pairwise(fit.den2, fit.den0, groups = gp.den, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

PW<-pairwise(fit.den2, fit.den1,groups = gp.den, covariate = gdf.den$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)


## Testing allometry further
#Homogeneity of Slopes test, parallel (fit.den1) or non-parallel *includes covariate
gp.den<-interaction(gdf.den$morph)

PW<-pairwise(fit.den2, fit.den1,groups = gp.den, covariate = gdf.den$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

#Significant differences in vector angles (where the shape trajectories are going, convergence or divergence) *included covariate
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

### Plotting
palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB

test<-plotAllometry(fit.den2,gdf.den$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.den$morph)], xlab="Predictor")

spreds<-shape.predictor(gdf.den$coords, test$PredLine, Intercept = TRUE, predmin = min(test$PredLine),predmax = max(test$PredLine))
M <- mshape(gdf.den$coords)
plotRefToTarget(M, spreds$predmin, mag=2) #left
plotRefToTarget(M, spreds$predmax, mag=2) #right



## Plotting by morph (with size effects.)
palettec<-c("#74add1","#f46d43","#a50026","#313695")
gp.den<-gdf.den$morph
col.gp<-palettec
names(col.gp)<-levels(gp.den)
col.gp.den<-col.gp[match(gp.den,names(col.gp))]
TS.den<-gm.prcomp(gdf.den$coords)
plot(TS.den,pch=21, bg=col.gp.den)

scatterplot(TS.den$x[,1],TS.den$x[,2],groups = gp.den, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.1,0.1),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TS.den$x[,2],TS.den$x[,3],groups = gp.den, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.10,0.12),ylim=c(-0.10,0.12),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(gdf.den$coords)

### Plot Warps
##### What are shape differences for each PC.
ref.den<-mshape(gdf.den$coords)
PC <- TS.den$x[,1] 
preds <- shape.predictor(gdf.den$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.den, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.den, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
#### Run also for other PC.


## Allometry free shape
## Plotting by morph (without size effects.) 
### Can only be used if there is no difference in allometry between the morphs.
den.Anova<-procD.lm(coords~logcsize, data = gdf.den,iter=499,RRPP=TRUE)
shape.resid<-arrayspecs(den.Anova$residuals,p=dim(gpa.lands$coords)[1],k=dim(gpa.lands$coords)[2])
adj.shape.den<-shape.resid + array(gpa.lands$consensus,dim(shape.resid))
TSadj.den<-gm.prcomp(adj.shape.den)
plot(TSadj.den,pch=21,bg=col.gp.den)

scatterplot(TSadj.den$x[,1],TSadj.den$x[,2],groups = gp.den, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.1,0.1),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TSadj.den$x[,2],TSadj.den$x[,3],groups = gp.den, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.1,0.1),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(adj.shape.den)

### Plot Warps (for size adjusted)
##### What are shape differences for each PC.
ref.den<-mshape(adj.shape.den)
PC <- TSadj.den$x[,1]
preds <- shape.predictor(adj.shape.den, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.den, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.den, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
#### Run also for other PC.

## Allometry Regression
fit.size.null<-procD.lm(coords~logcsize,data=gdf.den)
fit.size<-procD.lm(coords~logcsize*morph,groups=denCh$morph,data=gdf.den)
A<-plotAllometry(fit.size,size=gpa.lands$Csize,logsz=TRUE, method="RegScore",pch=19,col=col.gp.den)
B<-plotAllometry(fit.size,size=gpa.lands$Csize,logsz=TRUE, method="PredLine",pch=19,col=col.gp.den)
C<-plotAllometry(fit.size,size=gpa.lands$Csize,logsz=TRUE, method="size.shape",pch=19,col=col.gp.den)
D<-plotAllometry(fit.size,size=gpa.lands$Csize,logsz=TRUE, method="CAC",pch=19,col=col.gp.den)


####################################################################
#####################################################################
############### Maxilla morphology ################################

###Reading data and general data preparation

max<-readland.tps("Maxilla_LM.tps", negNA = TRUE, specID = c("imageID"))

max1<-estimate.missing(max,method = "TPS")

Need_flipping_max<-max1[,,131]
Flipped_max<-rotate.coords(Need_flipping_max, "flipY")
max1[,,131]<-Flipped_max

max1<-rotate.coords(max1, type = c("flipY"))

##### LM: 15, 16, 17
sliding1 <-matrix(c(14, 15, 16, 15, 16, 7, 7, 17, 5),ncol = 3,byrow = TRUE) # Landmarks 15, 16, 17
max.land1<-gpagen(max1, curves = sliding1)

maxCh<- read.csv("Maxilla_Ch.csv", stringsAsFactors = TRUE, header=TRUE)

gdf.max<-geomorph.data.frame(max.land1, morph = maxCh$morph, 
                             csize = max.land1$Csize,logcsize = log(max.land1$Csize),
                             ind = dimnames(max.land1$coords)[[3]], length = maxCh$length_cm,
                             weight = maxCh$weight_g, id = maxCh$id_me, sex = maxCh$sex, 
                             age = maxCh$age, loglength = log(maxCh$length_cm))

gdf.max<-geomorph.data.frame(max.land1, morph = maxCh$morph, 
                             csize = max.land1$Csize,logcsize = log(max.land1$Csize),
                             ind = dimnames(max.land1$coords)[[3]], length = maxCh$length_cm,
                             weight = maxCh$weight_g, id = maxCh$id_me, sex = maxCh$sex, 
                             age = maxCh$age, loglength = log(maxCh$length_cm), 
                             teethd = maxCh$teeth_den, teethm =maxCh$teeth_max, 
                             teethpr = maxCh$teeth_pre, teethpa = maxCh$teeth_pal, 
                             teethv = maxCh$teeth_vo, teethg = maxCh$teeth_glos, 
                             teethang = maxCh$teeth_ang)


maxCh$logCsize<-gdf.max$logcsize
ggplot(maxCh, aes(x=length_cm, y=logCsize)) +geom_point()+ geom_smooth(method="lm", se=T) + xlab("Length")

plotAllSpecimens(max.land1$coords, mean = TRUE)
plotOutliers(max.land1$coords,inspect.outliers = TRUE)

plotOutliers(max.land1$coords, groups = gdf.max$morph, inspect.outliers = TRUE)

########## Testing morph effects ###########################
fit.max0<-procD.lm(coords ~ logcsize, data = gdf.max,iter = 999, RRPP=TRUE)
fit.max1<-procD.lm(coords ~ logcsize + morph, data = gdf.max,iter = 999, RRPP=TRUE)
fit.max2<-procD.lm(coords ~ logcsize * morph, data = gdf.max,iter = 999, RRPP=TRUE)

anova(fit.max0, fit.max1)
anova(fit.max1, fit.max2)

summary(fit.max1)
summary(fit.max2)

gp.max<-interaction(gdf.max$morph)
PW2<-pairwise(fit.max1, fit.max0,groups = gp.max, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means


## Testing allometry further
gp.max<-interaction(gdf.max$morph)

#Homogeneity of Slopes test, parallel (fit.den1) or non-parallel *includes covariate
PW<-pairwise(fit.max2, fit.max1,groups = gp.max, covariate = gdf.max$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

#Significant differences in vector angles (where the shape trajectories are going, convergence or divergence) *included covariate
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB

test<-plotAllometry(fit.max2,gdf.max$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.max$morph)], xlab="Predictor")
spreds<-shape.predictor(gdf.max$coords, test$PredLine, Intercept = TRUE, predmin = min(test$PredLine),predmax = max(test$PredLine))

M <- mshape(gdf.max$coords)

plotRefToTarget(M, spreds$predmin, mag=1) #left
plotRefToTarget(M, spreds$predmax, mag=1) #right

## Plotting by morph (with size effects)
palettec<-c("#74add1","#f46d43","#a50026","#313695")
gp.max<-gdf.max$morph
col.gp<-palettec
names(col.gp)<-levels(gp.max)
col.gp.max<-col.gp[match(gp.max,names(col.gp))]
TS.max<-gm.prcomp(gdf.max$coords)
plot(TS.max,pch=21, bg=col.gp.max)

scatterplot(TS.max$x[,1],TS.max$x[,2],groups = gp.max, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.1,0.1),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TS.max$x[,2],TS.max$x[,3],groups = gp.max, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.1,0.1),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(gdf.max$coords)

### Plot Warps
ref.max<-mshape(gdf.max$coords)
PC <- TS.max$x[,1]
preds <- shape.predictor(gdf.max$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.max, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.max, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)

## Allometry free shape
## Plotting by morph (without size effects.) 
max.Anova<-procD.lm(coords~logcsize, data = gdf.max,iter=499,RRPP=TRUE)
shape.resid<-arrayspecs(max.Anova$residuals,p=dim(max.land1$coords)[1],k=dim(max.land1$coords)[2])
adj.shape.max<-shape.resid + array(max.land1$consensus,dim(shape.resid))
TSadj.max<-gm.prcomp(adj.shape.max)
plot(TSadj.max,pch=21,bg=col.gp.max)

scatterplot(TSadj.max$x[,1],TSadj.max$x[,2],groups = gp.max, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.12,0.12),ylim=c(-0.1,0.1),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TSadj.max$x[,2],TSadj.max$x[,3],groups = gp.max, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.12,0.12),ylim=c(-0.1,0.1),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(adj.shape.max)

### Plot Warps (for size adjusted)
ref.max<-mshape(adj.shape.max)
PC <- TSadj.max$x[,1]
preds <- shape.predictor(adj.shape.max, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.max, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.max, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)

## Allometry Regression
fit.size.null<-procD.lm(coords~logcsize,data=gdf.max)
fit.size<-procD.lm(coords~logcsize * morph, groups=maxCh$morph,data=gdf.max)
A<-plotAllometry(fit.size,size=max.land1$Csize,logsz=TRUE, method="RegScore",pch=19,col=col.gp.max)
B<-plotAllometry(fit.size,size=max.land1$Csize,logsz=TRUE, method="PredLine",pch=19,col=col.gp.max)
C<-plotAllometry(fit.size,size=max.land1$Csize,logsz=TRUE, method="size.shape",pch=19,col=col.gp.max)
D<-plotAllometry(fit.size,size=max.land1$Csize,logsz=TRUE, method="CAC",pch=19,col=col.gp.max)


####################################################################
#####################################################################
############### Articular-angular morphology ################################

###Reading data and general data preparation

art<-readland.tps("Articular_Angular_LM.tps", negNA = TRUE, specID = c("imageID"))

art1<-estimate.missing(art,method = "TPS")

#Flipping
x<-c(26, 28, 31, 63, 116)
Need_flipping_art<-art1[,,x]
Flipped_art<-rotate.coords(Need_flipping_art, "flipX")
art1[,,x]<-Flipped_art


#### LM: 9, 13, 14
sliding2<-matrix(c(12, 13, 14, 13, 14, 1), ncol = 3, byrow = TRUE)

art.lands<-gpagen(art1, curves = sliding2)


artCh<- read.csv("Articular_Angular_Ch.csv", stringsAsFactors = TRUE, header=TRUE)


gdf.art<-geomorph.data.frame(art.lands, morph = artCh$morph, csize = art.lands$Csize,
                             logcsize = log(art.lands$Csize),ind = dimnames(art.lands$coords)[[3]],
                             length = artCh$length_cm,weight = artCh$weight_g,
                             id = artCh$id_me, sex = artCh$sex,
                             age = artCh$age, loglength = log(artCh$length_cm))


artCh$logCsize<-gdf.art$logcsize
ggplot(artCh, aes(x=length_cm, y=logCsize)) +geom_point()+ geom_smooth(method="lm", se=T) + xlab("Length")

plotAllSpecimens(art.lands$coords, mean = TRUE)
plotOutliers(art.lands$coords,inspect.outliers = TRUE)
outliers <- plotOutliers(gdf.art$coords, groups = gdf.art$morph, inspect.outliers = TRUE)


########## Testing morph effects ###########################
fit.art0<-procD.lm(coords ~ logcsize, data = gdf.art,iter = 999, RRPP=TRUE)
fit.art1<-procD.lm(coords ~ logcsize + morph, data = gdf.art,iter = 999, RRPP=TRUE)
fit.art2<-procD.lm(coords ~ logcsize * morph, data = gdf.art,iter = 999, RRPP=TRUE)

anova(fit.art0, fit.art1)
anova(fit.art1, fit.art2)

summary(fit.art1)
summary(fit.art2)

gp.art<-interaction(gdf.art$morph)
PW2<-pairwise(fit.art1, fit.art0,groups = gp.art, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means


## Testing allometry further
gp.art<-interaction(gdf.art$morph)
#Homogeneity of Slopes test, parallel (fit.den1) or non-parallel *includes covariate
PW<-pairwise(fit.art2, fit.art1,groups = gp.art, covariate = gdf.art$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
#Significant differences in vector angles (where the shape trajectories are going, convergence or divergence) *included covariate
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
test<-plotAllometry(fit.art2,gdf.art$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.art$morph)], xlab="Predictor")

spreds<-shape.predictor(gdf.art$coords, test$PredLine, Intercept = TRUE, predmin = min(test$PredLine),predmax = max(test$PredLine))

M <- mshape(gdf.art$coords)

plotRefToTarget(M, spreds$predmin, mag=1) #left
plotRefToTarget(M, spreds$predmax, mag=1) #right

## Plotting by morph (with size effects.)
palettec<-c("#74add1","#f46d43","#a50026","#313695")
gp.art<-gdf.art$morph
col.gp<-palettec
names(col.gp)<-levels(gp.art)
col.gp.art<-col.gp[match(gp.art,names(col.gp))]
TS.art<-gm.prcomp(gdf.art$coords)
plot(TS.art,pch=21, bg=col.gp.art)

scatterplot(TS.art$x[,1],TS.art$x[,2],groups = gp.art, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.15,0.15),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TS.art$x[,2],TS.art$x[,3],groups = gp.art, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.15,0.15),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(gdf.art$coords) 

### Plot Warps
ref.art<-mshape(gdf.art$coords)
PC <- TS.art$x[,1]
preds <- shape.predictor(gdf.art$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.art, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.art, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry free shape
## Plotting by morph (without size effects.) 
art.Anova<-procD.lm(coords~logcsize, data = gdf.art,iter=499,RRPP=TRUE)
shape.resid<-arrayspecs(art.Anova$residuals,p=dim(art.lands$coords)[1],k=dim(art.lands$coords)[2])
adj.shape.art<-shape.resid + array(art.lands$consensus,dim(shape.resid))
TSadj.art<-gm.prcomp(adj.shape.art)
plot(TSadj.art,pch=21,bg=col.gp.art)

scatterplot(TSadj.art$x[,1],TSadj.art$x[,2],groups =gp.art, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.16,0.12),ylim=c(-0.1,0.1),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TSadj.art$x[,2],TSadj.art$x[,3],groups =gp.art, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.16,0.12),ylim=c(-0.1,0.1),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(adj.shape.art)

### Plot Warps (for size adjusted)
ref.art<-mshape(adj.shape.art)
PC <- TSadj.art$x[,1] # Run for other PC
preds <- shape.predictor(adj.shape.art, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.art, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.art, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry Regression
fit.size.null<-procD.lm(coords~logcsize,data=gdf.art)
fit.size<-procD.lm(coords~logcsize*morph,groups=artCh$morph,data=gdf.art)
A<-plotAllometry(fit.size,size=art.lands$Csize,logsz=TRUE, method="RegScore",pch=19,col=col.gp.art)
B<-plotAllometry(fit.size,size=art.lands$Csize,logsz=TRUE, method="PredLine",pch=19,col=col.gp.art)
C<-plotAllometry(fit.size,size=art.lands$Csize,logsz=TRUE, method="size.shape",pch=19,col=col.gp.art)
D<-plotAllometry(fit.size,size=art.lands$Csize,logsz=TRUE, method="CAC",pch=19,col=col.gp.art)


####################################################################
#####################################################################
############### Quadrate morphology ################################

###Reading data and general data preparation

qu<-readland.tps("Quadrate_LM.tps", negNA = TRUE, specID = c("imageID"))

qu1<-estimate.missing(qu,method = "TPS")

x<-c(26,44,49:50,52,58:60,66,68,70,228)
Need_flipping_qu<-qu1[,,x]
Flipped_qu<-rotate.coords(Need_flipping_qu, "flipX")
qu1[,,x]<-Flipped_qu


##### LM: 4
sliding<-matrix(c(3,4,6),ncol=3,byrow=TRUE)
qu.land<-gpagen(qu1, curves = sliding)

quCh<- read.csv("Quadrate_Ch.csv", stringsAsFactors = TRUE, header=TRUE) # Yes use the same data set as Articular angular


gdf.qu<-geomorph.data.frame(qu.land, morph = quCh$morph, csize = qu.land$Csize,
                            logcsize = log(qu.land$Csize),ind = dimnames(qu.land$coords)[[3]],
                            length = quCh$length_cm,weight = quCh$weight_g, id = quCh$id_me,
                            sex = quCh$sex, age = quCh$age, loglength = log(quCh$length_cm))

quCh$logCsize<-gdf.qu$logcsize
ggplot(quCh, aes(x=length_cm, y=logCsize)) +geom_point()+ geom_smooth(method="lm", se=T) + xlab("Length")

plotAllSpecimens(qu.land$coords, mean=TRUE)
plotOutliers(qu.land$coords,inspect.outliers = TRUE)
outliers <- plotOutliers(qu.land$coords, groups = gdf.qu$morph, inspect.outliers = TRUE)

########## Testing morph effects ###########################
fit.qu0<-procD.lm(coords ~ logcsize, data = gdf.qu,iter = 999, RRPP=TRUE)
fit.qu1<-procD.lm(coords ~ logcsize + morph, data = gdf.qu,iter = 999, RRPP=TRUE)
fit.qu2<-procD.lm(coords ~ logcsize * morph, data = gdf.qu,iter = 999, RRPP=TRUE)

anova(fit.qu0, fit.qu1)
anova(fit.qu1, fit.qu2)

summary(fit.qu1)
summary(fit.qu2)

gp.qu<-interaction(gdf.qu$morph)
PW2<-pairwise(fit.qu1, fit.qu0,groups = gp.qu, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means


## Testing allometry further
gp.qu<-interaction(gdf.qu$morph)

#Homogeneity of Slopes test, parallel (fit.den1) or non-parallel *includes covariate
PW<-pairwise(fit.qu2, fit.qu1,groups = gp.qu, covariate = gdf.qu$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

#Significant differences in vector angles (where the shape trajectories are going, convergence or divergence) *included covariate
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB

test<-plotAllometry(fit.qu2,gdf.qu$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.qu$morph)], xlab="Predictor")

spreds<-shape.predictor(gdf.qu$coords, test$PredLine, Intercept = TRUE, predmin = min(test$PredLine),predmax = max(test$PredLine))
M <- mshape(gdf.qu$coords)
plotRefToTarget(M, spreds$predmin, mag=1) #left
plotRefToTarget(M, spreds$predmax, mag=1) #right


## Plotting by morph (with size effects.)
palettec<-c("#74add1","#f46d43","#a50026","#313695")
gp.qu<-gdf.qu$morph
col.gp<-palettec
names(col.gp)<-levels(gp.qu)
col.gp.qu<-col.gp[match(gp.qu,names(col.gp))]
TS.qu<-gm.prcomp(gdf.qu$coords)
plot(TS.qu,pch=21, bg=col.gp.qu)

scatterplot(TS.qu$x[,1],TS.qu$x[,2],groups = gp.qu, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.21,0.15),ylim=c(-0.17,0.15),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TS.qu$x[,2],TS.qu$x[,3],groups = gp.qu, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.21,0.15),ylim=c(-0.17,0.15),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(gdf.qu$coords)

### Plot Warps
ref.qu<-mshape(gdf.qu$coords)
PC <- TS.qu$x[,1]
preds <- shape.predictor(gdf.qu$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.qu, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.qu, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry free shape
## Plotting by morph (without size effects.) 
qu.Anova<-procD.lm(coords~logcsize, data = gdf.qu,iter=499,RRPP=TRUE)
shape.resid<-arrayspecs(qu.Anova$residuals,p=dim(qu.land$coords)[1],k=dim(qu.land$coords)[2])
adj.shape.qu<-shape.resid + array(qu.land$consensus,dim(shape.resid))
TSadj.qu<-gm.prcomp(adj.shape.qu)
plot(TSadj.qu,pch=21,bg=col.gp.qu)

scatterplot(TSadj.qu$x[,1],TSadj.qu$x[,2],groups = gp.qu, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.2,0.2),ylim=c(-0.18,0.15),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TSadj.qu$x[,2],TSadj.qu$x[,3],groups = gp.qu, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.2,0.2),ylim=c(-0.18,0.15),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(adj.shape.qu)

### Plot Warps (for size adjusted)
ref.qu<-mshape(adj.shape.qu)
PC <- TSadj.qu$x[,1] # Run for other PC
preds <- shape.predictor(adj.shape.qu, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.qu, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.qu, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry Regression
fit.size.null<-procD.lm(coords~logcsize,data=gdf.qu)
fit.size<-procD.lm(coords~logcsize*morph,groups=quCh$morph,data=gdf.qu)
A<-plotAllometry(fit.size,size=qu.land$Csize,logsz=TRUE, method="RegScore",pch=19,col=col.gp.qu)
B<-plotAllometry(fit.size,size=qu.land$Csize,logsz=TRUE, method="PredLine",pch=19,col=col.gp.qu)
C<-plotAllometry(fit.size,size=qu.land$Csize,logsz=TRUE, method="size.shape",pch=19,col=col.gp.qu)
D<-plotAllometry(fit.size,size=qu.land$Csize,logsz=TRUE, method="CAC",pch=19,col=col.gp.qu)


####################################################################
#####################################################################
############### Premaxilla morphology ################################

###Reading data and general data preparation

prem<-readland.tps("Premaxilla_LM.tps", negNA = TRUE, specID = c("imageID"))

prem1<-estimate.missing(prem,method = "TPS")


x<-c(1,9,28,64,79,82,87,97,99,198,206,216,219)
Need_flipping_prem<-prem1[,,x]
Flipped_prem<-rotate.coords(Need_flipping_prem, "flipY")
prem1[,,x]<-Flipped_prem


prem1<-rotate.coords(prem1, type = c("rotateC"))

#### No sliding LM
pre.land<-gpagen(prem1)


preCh<- read.csv("Premaxilla_Ch.csv", stringsAsFactors = TRUE, header=TRUE) # Yes same data set as Articular angular


gdf.pre<-geomorph.data.frame(pre.land, morph = preCh$morph, csize = pre.land$Csize,
                             logcsize = log(pre.land$Csize),ind = dimnames(pre.land$coords)[[3]],
                             length = preCh$length_cm,weight = preCh$weight_g, id = preCh$id_me,
                             sex = preCh$sex, age = preCh$age, loglength = log(preCh$length_cm))

gdf.pre<-geomorph.data.frame(pre.land, morph = preCh$morph, csize = pre.land$Csize,
                             logcsize = log(pre.land$Csize),ind = dimnames(pre.land$coords)[[3]],
                             length = preCh$length_cm,weight = preCh$weight_g, id = preCh$id_me,
                             sex = preCh$sex, age = preCh$age, loglength = log(preCh$length_cm),
                             teethd = preCh$teeth_den, teethm =preCh$teeth_max,
                             teethpr = preCh$teeth_pre, teethpa = preCh$teeth_pal,
                             teethv = preCh$teeth_vo, teethg = preCh$teeth_glos)


preCh$logCsize<-gdf.pre$logcsize
ggplot(preCh, aes(x=length_cm, y=logCsize)) +geom_point()+ geom_smooth(method="lm", se=T) + xlab("Length")


plotAllSpecimens(pre.land$coords, mean = TRUE)
plotOutliers(pre.land$coords,inspect.outliers = TRUE)
outliers <- plotOutliers(pre.land$coords, groups = gdf.pre$morph, inspect.outliers = TRUE)


########## Testing morph effects ###########################
fit.pre0<-procD.lm(coords ~ logcsize, data = gdf.pre,iter = 999, RRPP=TRUE)
fit.pre1<-procD.lm(coords ~ logcsize + morph, data = gdf.pre,iter = 999, RRPP=TRUE)
fit.pre2<-procD.lm(coords ~ logcsize * morph, data = gdf.pre,iter = 999, RRPP=TRUE)

anova(fit.pre0, fit.pre1)
anova(fit.pre1, fit.pre2)

summary(fit.pre1)
summary(fit.pre2)

gp.pre<-interaction(gdf.pre$morph)
PW2<-pairwise(fit.pre1, fit.pre0,groups = gp.pre, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means

PW<-pairwise(fit.pre2, fit.pre1,groups = gp.pre, covariate = gdf.pre$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means


## Testing allometry further
gp.pre<-interaction(gdf.pre$morph)

#Homogeneity of Slopes test, parallel (fit.den1) or non-parallel *includes covariate
PW<-pairwise(fit.pre2, fit.pre1,groups = gp.pre, covariate = gdf.pre$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

#Significant differences in vector angles (where the shape trajectories are going, convergence or divergence) *included covariate
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB

test<-plotAllometry(fit.pre2,gdf.pre$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.pre$morph)], xlab="Predictor")

plot(gdf.pre$logcsize,(test$PredLine)*-1,pch = 19, col = palettething[as.numeric(gdf.pre$morph)],xlab = c("Predictor"),ylab = c("PC 1 for fitted values"))

Predladj<-(test$PredLine)*-1
spreds<-shape.predictor(gdf.pre$coords, Predladj, Intercept = TRUE, predmin = min(Predladj),predmax = max(Predladj))

M <- mshape(gdf.pre$coords)

plotRefToTarget(M, spreds$predmin, mag=1) #left
plotRefToTarget(M, spreds$predmax, mag=1) #right


## Plotting by morph (with size effects.)
palettec<-c("#74add1","#f46d43","#a50026","#313695")
gp.pre<-gdf.pre$morph
col.gp<-palettec
names(col.gp)<-levels(gp.pre)
col.gp.pre<-col.gp[match(gp.pre,names(col.gp))]
TS.pre<-gm.prcomp(gdf.pre$coords)
plot(TS.pre,pch=21, bg=col.gp.pre)

scatterplot(TS.pre$x[,1],TS.pre$x[,2],groups = gp.pre, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.22),ylim=c(-0.17,0.2),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TS.pre$x[,2],TS.pre$x[,3],groups = gp.pre, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.22),ylim=c(-0.17,0.2),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(gdf.pre$coords)


### Plot Warps
ref.pre<-mshape(gdf.pre$coords)
PC <- TS.pre$x[,1] 
preds <- shape.predictor(gdf.pre$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.pre, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.pre, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry free shape
## Plotting by morph (without size effects.) 
pre.Anova<-procD.lm(coords~logcsize, data = gdf.pre,iter=499,RRPP=TRUE)
shape.resid<-arrayspecs(pre.Anova$residuals,p=dim(pre.land$coords)[1],k=dim(pre.land$coords)[2])
adj.shape.pre<-shape.resid + array(pre.land$consensus,dim(shape.resid))
TSadj.pre<-gm.prcomp(adj.shape.pre)
plot(TSadj.pre,pch=21,bg=col.gp.pre)

scatterplot(TSadj.pre$x[,1],TSadj.pre$x[,2],groups = gp.pre, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.15,0.15),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TSadj.pre$x[,2],TSadj.pre$x[,3],groups = gp.pre, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.15,0.15),ylim=c(-0.15,0.15),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(adj.shape.pre)

### Plot Warps (for size adjusted)
ref.pre<-mshape(adj.shape.pre)
PC <- TSadj.pre$x[,1]
preds <- shape.predictor(adj.shape.pre, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.pre, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.pre, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry Regression
fit.size.null<-procD.lm(coords~logcsize,data=gdf.pre)
fit.size<-procD.lm(coords~logcsize*morph,groups=preCh$morph,data=gdf.pre)
A<-plotAllometry(fit.size,size=pre.land$Csize,logsz=TRUE, method="RegScore",pch=19,col=col.gp.pre)
B<-plotAllometry(fit.size,size=pre.land$Csize,logsz=TRUE, method="PredLine",pch=19,col=col.gp.pre)
C<-plotAllometry(fit.size,size=pre.land$Csize,logsz=TRUE, method="size.shape",pch=19,col=col.gp.pre)
D<-plotAllometry(fit.size,size=pre.land$Csize,logsz=TRUE, method="CAC",pch=19,col=col.gp.pre)


####################################################################
#####################################################################
############### Supramaxilla morphology ################################

###Reading data and general data preparation

sup<-readland.tps("Supramaxilla_LM.tps", negNA = TRUE, specID = c("imageID"))

sup1<-estimate.missing(sup,method = "TPS")

x<-c(2,23,27,30,37,43,45:47,51,53,55,57:59,64,67,69,75:76,81,90,94,100,150,
     152,156,169,179,207,213)
Need_flipping_sup<-sup1[,,x]
Flipped_sup<-rotate.coords(Need_flipping_sup, "flipY")
sup1[,,x]<-Flipped_sup


### No sliding LM
sup.land<-gpagen(sup1)

supCh<- read.csv("Supramaxilla_Ch.csv", stringsAsFactors = TRUE, header=TRUE)


gdf.sup<-geomorph.data.frame(sup.land, morph = supCh$morph, csize = sup.land$Csize,
                             logcsize = log(sup.land$Csize),ind = dimnames(sup.land$coords)[[3]],
                             length = supCh$length_cm, weight = supCh$weight_g,
                             id = supCh$id_me, sex = supCh$sex, age = supCh$age,
                             loglength = log(supCh$length_cm))


supCh$logCsize<-gdf.sup$logcsize
ggplot(supCh, aes(x=length_cm, y=logCsize)) +geom_point()+ geom_smooth(method="lm", se=T) + xlab("Length")

plotAllSpecimens(sup.land$coords, mean = TRUE)
plotOutliers(sup.land$coords,inspect.outliers = TRUE)
outliers <- plotOutliers(sup.land$coords, groups = gdf.sup$morph, inspect.outliers = TRUE)

########## Testing morph effects ###########################
fit.sup0<-procD.lm(coords ~ logcsize, data = gdf.sup,iter = 999, RRPP=TRUE)
fit.sup1<-procD.lm(coords ~ logcsize + morph, data = gdf.sup,iter = 999, RRPP=TRUE)
fit.sup2<-procD.lm(coords ~ logcsize * morph, data = gdf.sup,iter = 999, RRPP=TRUE)

anova(fit.sup0, fit.sup1)
anova(fit.sup1, fit.sup2)

summary(fit.sup1)
summary(fit.sup2)

gp.sup<-interaction(gdf.sup$morph)
PW2<-pairwise(fit.sup1, fit.sup0,groups = gp.sup, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means


## Testing allometry further
gp.sup<-interaction(gdf.sup$morph)
#Homogeneity of Slopes test, parallel (fit.den1) or non-parallel *includes covariate
PW<-pairwise(fit.sup2, fit.sup1,groups = gp.sup, covariate = gdf.sup$logcsize)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

#Significant differences in vector angles (where the shape trajectories are going, convergence or divergence) *included covariate
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB

test<-plotAllometry(fit.sup2,gdf.sup$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.sup$morph)], xlab="Predictor")

spreds<-shape.predictor(gdf.sup$coords, test$PredLine, Intercept = TRUE, predmin = min(test$PredLine),predmax = max(test$PredLine))
M <- mshape(gdf.sup$coords)
plotRefToTarget(M, spreds$predmin, mag=1) #left
plotRefToTarget(M, spreds$predmax, mag=1) #right

## Plotting by morph (with size effects.)
palettec<-c("#74add1","#f46d43","#a50026","#313695")
gp.sup<-gdf.sup$morph
col.gp<-palettec
names(col.gp)<-levels(gp.sup)
col.gp.sup<-col.gp[match(gp.sup,names(col.gp))]
TS.sup<-gm.prcomp(gdf.sup$coords)
plot(TS.sup,pch=21, bg=col.gp.sup)

scatterplot(TS.sup$x[,1],TS.sup$x[,2],groups = gp.sup, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.25,0.2),ylim=c(-0.2,0.2),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TS.sup$x[,2],TS.sup$x[,3],groups = gp.sup, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.25,0.2),ylim=c(-0.2,0.2),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(gdf.sup$coords)

### Plot Warps
ref.sup<-mshape(gdf.sup$coords)
PC <- TS.sup$x[,1] 
preds <- shape.predictor(gdf.sup$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.sup, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.sup, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry free shape
## Plotting by morph (without size effects.) 
sup.Anova<-procD.lm(coords~logcsize, data = gdf.sup,iter=499,RRPP=TRUE)
shape.resid<-arrayspecs(sup.Anova$residuals,p=dim(sup.land$coords)[1],k=dim(sup.land$coords)[2])
adj.shape.sup<-shape.resid + array(sup.land$consensus,dim(shape.resid))
TSadj.sup<-gm.prcomp(adj.shape.sup)
plot(TSadj.sup,pch=21,bg=col.gp.sup)

scatterplot(TSadj.sup$x[,1],TSadj.sup$x[,2],groups = gp.sup, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.2,0.2),ylim=c(-0.18,0.18),asp=1, xlab = "PC 1", ylab = "PC 2")
scatterplot(TSadj.sup$x[,2],TSadj.sup$x[,3],groups = gp.sup, grid = TRUE, smooth = FALSE, regLine = FALSE,col = palettec, pch=c(19,1,19,1,19,1), ellipse =list(levels=c(.95), robust=TRUE, fill=TRUE, fill.alpha=0.2),xlim=c(-0.2,0.2),ylim=c(-0.18,0.18),asp=1, xlab = "PC 2", ylab = "PC 3")

### Information on each PC
gm.prcomp(adj.shape.sup)

### Plot Warps (for size adjusted)
ref.sup<-mshape(adj.shape.sup)
PC <- TSadj.sup$x[,1]
preds <- shape.predictor(adj.shape.sup, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(ref.sup, preds$pred1, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)
plotRefToTarget(ref.sup, preds$pred2, method = "TPS", gridPars=gridPar(pt.size=1,pt.bg="black"),mag=1)


## Allometry Regression
fit.size.null<-procD.lm(coords~logcsize,data=gdf.sup)
fit.size<-procD.lm(coords~logcsize*morph,groups=supCh$morph,data=gdf.sup)
A<-plotAllometry(fit.size,size=sup.land$Csize,logsz=TRUE, method="RegScore",pch=19,col=col.gp.sup)
B<-plotAllometry(fit.size,size=sup.land$Csize,logsz=TRUE, method="PredLine",pch=19,col=col.gp.sup)
C<-plotAllometry(fit.size,size=sup.land$Csize,logsz=TRUE, method="size.shape",pch=19,col=col.gp.sup)
D<-plotAllometry(fit.size,size=sup.land$Csize,logsz=TRUE, method="CAC",pch=19,col=col.gp.sup)


#####################################################
###### External shape data analysis ##########
#####################################################

###################################################################
################# External whole body ############################
###Reading data and general data preparation

dat<-readland.tps("Whole_body_LM.tps", negNA = TRUE, specID = c("imageID"))

#### Estimate missing landmarks
dat<-estimate.missing(dat,method = "TPS")


#### Flipping of individuals where it is needed.
x<-c(40,100)
Need_flipping_WB<-dat[,,x]
Flipped_WB<-rotate.coords(Need_flipping_WB, "flipX")
dat[,,x]<-Flipped_WB


#### Specifying which landmarks are sliding. LM =6,18,20,24,25,26,27
sliding<-matrix(c(5,6,7,17,18,19,18,19,20,26,24,25,24,25,26,25,26,27,26,27,24),ncol=3,byrow=TRUE)

####Performing Generlized Procrustes analysis
dat.lands<-gpagen(dat, curves = sliding)


### Adding character information
datCh<-read.csv("External_Ch.csv", stringsAsFactors = TRUE, header = TRUE, sep = ",")

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
palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
test<-plotAllometry(fit.dat2,gdf.dat$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.dat$morph)],xlab="Predictor")


spreds<-shape.predictor(gdf.dat$coords, test$PredLine, Intercept = TRUE, predmin = min(test$PredLine),predmax = max(test$PredLine))
M <- mshape(gdf.dat$coords)
plotRefToTarget(M, spreds$predmin, mag=1) #left
plotRefToTarget(M, spreds$predmax, mag=1) #right


## Plotting by morph (with size effects.)
palettec<-c("#74add1","#f46d43","#a50026","#313695")
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

hau<-readland.tps("Head_LM.tps", negNA = TRUE, specID = c("imageID"))

hau<-estimate.missing(hau,method = "TPS")


x<-c(40,100)
Need_flipping_head<-hau[,,x]
Flipped_head<-rotate.coords(Need_flipping_head, "flipX")
hau[,,x]<-Flipped_head

##### LM: 10, 11, 12, 13.
slidingH<-matrix(c(13,10,11,10,11,12,11,12,13,12,13,10),ncol=3,byrow=TRUE)

hau.lands<-gpagen(hau, curves = slidingH)


hauCh<-read.csv("External_Ch.csv", stringsAsFactors = TRUE, header = TRUE, sep = ",")

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

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
plotAllometry(fit.hau2,gdf.hau$logcsize,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.hau$morph)], xlab="Predictor") 

spreds<-shape.predictor(gdf.hau$coords, test$PredLine, Intercept = TRUE, predmin = min(test$PredLine),predmax = max(test$PredLine))
M <- mshape(gdf.hau$coords)
plotRefToTarget(M, spreds$predmin, mag=1) #left
plotRefToTarget(M, spreds$predmax, mag=1) #right


## Plotting by morph (with size effects.)
palettec<-c("#74add1","#f46d43","#a50026","#313695")
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


###################################################
################################################
###########################################
###################### Checking the relationships #######################

##### Whole Body ######

#Whole body
## Relationship between head centroid size and body centroid
#(only individuals with external photo)

## data Preparation
body<-readland.tps("Whole_Body_LM.tps", negNA = TRUE)
body1<-estimate.missing(body,method = "TPS")
sliding<-matrix(c(5,6,7,18,19,20,19,20,21,28,25,26,25,26,27,26,27,28,27,28,25),ncol=3,byrow=TRUE)
align<-gpagen(body1,curves=sliding)

pcran<-c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,35,36,37,38) #landmarks on postcranial shape
combinedtpshead<-body1[-pcran,,] #subset landmarks not in pcran
slidingh<-matrix(c(13,10,11,10,11,12,11,12,13,12,13,10),ncol=3,byrow=TRUE)
alignhead<-gpagen(combinedtpshead,curves=slidingh)

bodyxl<- read.csv("External_Ch.csv", stringsAsFactors = TRUE, header=TRUE)
bodyxl$loglength_cm<-log(bodyxl$length_cm)

Csizeb<-log(align$Csize)
Csizeh<-log(alignhead$Csize)
FL<-bodyxl$loglength_cm
age<-bodyxl$age

### Testing morph effects and allometry
fit.body0<-procD.lm(Csizeb ~ Csizeh,iter = 999, RRPP=TRUE)
fit.body1<-procD.lm(Csizeb ~ Csizeh + bodyxl$morph,iter = 999, RRPP=TRUE)
fit.body2<-procD.lm(Csizeb ~ Csizeh*bodyxl$morph,iter = 999, RRPP=TRUE)

anova(fit.body0,fit.body1)
anova(fit.body1,fit.body2)


palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
plotAllometry(fit.body2,Csizeb,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(bodyxl$morph)]) 


PW<-pairwise(fit.body2, fit.body1,groups = bodyxl$morph, covariate = Csizeh)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
PW<-pairwise(fit.body2, fit.body1,groups = bodyxl$morph, covariate = Csizeh)
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")


###Correlation between whole body and head
cor.test(Csizeb, Csizeh, method=c("pearson"))
plot(Csizeb, Csizeh, col=bodyxl$morph)

testdaat<-data.frame(Csizeb, Csizeh, FL, bodyxl$age, bodyxl$morph, bodyxl$sex)
m1<-ggplot(data = testdaat, aes(x=Csizeb, y=Csizeh, color=bodyxl.morph)) + geom_point() + labs(x="CSbody", y="CShead", color="Morph")
m1 + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))


### Correlation between head and FL
cor.test(FL, Csizeh, method=c("pearson"))
plot(FL, Csizeh)
summary(glm(Csizeb ~Csizeh * morph, data=bodyxl))
m2<-ggplot(data = testdaat, aes(x=FL, y=Csizeh, color=bodyxl.morph)) + geom_point() + labs(x="Fork length", y="CShead", color="Morph")
m2 + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))

#### Correlation between whole body and FL
cor.test(FL, Csizeb, method=c("pearson"))
plot(FL, Csizeb)
m3<-ggplot(data = testdaat, aes(x=FL, y=Csizeb, color=bodyxl.morph)) + geom_point() + labs(x="Fork length", y="CSbody", color="Morph")
m3 + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))


####### Effect of age
fit.body0<-procD.lm(FL ~ age,iter = 999, RRPP=TRUE)
fit.body1<-procD.lm(FL ~ age + bodyxl$morph,iter = 999, RRPP=TRUE)
fit.body2<-procD.lm(FL ~ age *bodyxl$morph,iter = 999, RRPP=TRUE)
anova(fit.body0,fit.body1)
anova(fit.body1, fit.body2)
summary(fit.body2)

newdata<-bodyxl[ !(bodyxl$num %in% c(12, 31, 36, 49, 126, 186)), ]
morph<-newdata$morph
age<-newdata$age

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB

plotAllometry(fit.body2,age,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(morph)])

plotAllometry(fit.body2,age,logsz = FALSE, method = c("RegScore"),pch = 19, col = palettething[as.numeric(morph)])

PW<-pairwise(fit.body2, fit.body1,groups = morph, covariate = age)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
PW<-pairwise(fit.body2, fit.body1,groups = morph, covariate = age)
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")


############################ Dentary #################################
####################################################################
## Relationship between Dentary centroid size and FL.

##Testing morph effects and allometry
fit.densize0<-procD.lm(logcsize ~ loglength, data = gdf.den, iter = 999, RRPP=TRUE)
fit.densize1<-procD.lm(logcsize ~ loglength + morph, data = gdf.den, iter = 999, RRPP=TRUE)
fit.densize2<-procD.lm(logcsize ~ loglength*morph, data = gdf.den, iter = 999, RRPP=TRUE)
fit.densize3<-procD.lm(logcsize ~ loglength*morph+age, data = gdf.den, iter = 999, RRPP=TRUE)
fit.densize4<-procD.lm(logcsize ~ loglength*morph*age, data = gdf.den, iter = 999, RRPP=TRUE)

anova(fit.densize0,fit.densize1)
anova(fit.densize1,fit.densize2)

gp.den<-interaction(gdf.den$morph)

PW2<-pairwise(fit.densize1, fit.densize0,groups = gp.den, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means

PW2<-pairwise(fit.densize1, fit.densize0,groups = gp.den, covariate = gdf.den$logcsize)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means


palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
plotAllometry(fit.densize2,gdf.den$loglength,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.den$morph)])

plotAllometry(fit.densize2,gdf.den$loglength,logsz = FALSE, method = c("RegScore"),pch = 19, col = palettething[as.numeric(gdf.den$morph)])

PW<-pairwise(fit.densize2, fit.densize1,groups = gdf.den$morph, covariate = gdf.den$loglength)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
PW<-pairwise(fit.densize2, fit.densize1,groups = gdf.den$morph, covariate = gdf.den$loglength)
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

#### Correlation between Dentary and FL
cor.test(gdf.den$logcsize, gdf.den$loglength, method=c("pearson"))
plot(gdf.den$logcsize, gdf.den$loglength)

datden<-data.frame(gdf.den$logcsize, gdf.den$loglength, gdf.den$morph, gdf.den$age)
names(datden)[1] <- "logcsize"
names(datden)[2] <- "loglength"
names(datden)[3] <- "morph"
names(datden)[4] <- "age" 

m4<-ggplot(data = datden, aes(x=loglength, y=logcsize, color=morph)) + geom_point() + labs(x="Log Fork length", y="Log CSden", color="Morph")
m4 + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))

lm1<-lm(logcsize~ loglength*morph, dat=datden)
anova(lm1)
summary(lm1)

m.lst <- lsmeans(lm1, "morph", var="loglength")
m.lst <- lstrends(lm1, "morph", var="loglength")
pairs(m.lst)

#### Effect of age
datdennew<-datden[!is.na(datden$age),]

fit.denage0<-procD.lm(logcsize ~ age,data = gdf.den,iter = 999, RRPP=TRUE)
fit.denage1<-procD.lm(logcsize ~ age + morph,data = gdf.den,iter = 999, RRPP=TRUE)
fit.denage2<-procD.lm(logcsize ~ age *morph,data = gdf.den,iter = 999, RRPP=TRUE)
anova(fit.denage0,fit.denage1)
anova(fit.denage1, fit.denage2)

summary(fit.denage1)

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB

plotAllometry(fit.denage2,datdennew$age,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(datdennew$morph)])


######################### Maxilla ###################################
###########################################################
## Relationship between Maxilla centroid size and FL.

#Testing morph effects and allometry
fit.maxsize0<-procD.lm(logcsize ~ loglength, data = gdf.max, iter = 999, RRPP=TRUE)
fit.maxsize1<-procD.lm(logcsize ~ loglength + morph, data = gdf.max, iter = 999, RRPP=TRUE)
fit.maxsize2<-procD.lm(logcsize ~ loglength*morph, data = gdf.max, iter = 999, RRPP=TRUE)
fit.maxsize3<-procD.lm(logcsize ~ loglength*morph+age, data = gdf.max, iter = 999, RRPP=TRUE)
fit.maxsize4<-procD.lm(logcsize ~ loglength*morph*age, data = gdf.max, iter = 999, RRPP=TRUE)

anova(fit.maxsize0,fit.maxsize1)
anova(fit.maxsize1,fit.maxsize2)


gp.max<-interaction(gdf.max$morph)
PW2<-pairwise(fit.maxsize1, fit.maxsize0,groups = gp.max, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means


palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
plotAllometry(fit.maxsize2,gdf.max$loglength,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.max$morph)])

plotAllometry(fit.maxsize2,gdf.max$loglength,logsz = FALSE, method = c("RegScore"),pch = 19, col = palettething[as.numeric(gdf.max$morph)])

PW<-pairwise(fit.maxsize2, fit.maxsize1,groups = gdf.max$morph, covariate = gdf.max$loglength)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

PW<-pairwise(fit.maxsize2, fit.maxsize1,groups = gdf.max$morph, covariate = gdf.max$loglength)
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")


##### Testing correlation
cor.test(gdf.max$logcsize, gdf.max$loglength, method=c("pearson"))
plot(gdf.max$logcsize, gdf.max$loglength)

datmax<-data.frame(gdf.max$logcsize, gdf.max$loglength, gdf.max$morph)
names(datmax)[1] <- "logcsize"
names(datmax)[2] <- "loglength"
names(datmax)[3] <- "morph"

m5<-ggplot(data = datmax, aes(x=loglength, y=logcsize, color=morph)) + geom_point() + labs(x="Log Fork length", y="Log CSmax", color="Morph")
m5 + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))

lm2<-lm(logcsize~loglength*morph, dat=datmax)
anova(lm2)
summary(lm2)

m.lst <- lstrends(lm2, "morph", var="loglength")
pairs(m.lst)

#### Effect of age 
fit.maxage0<-procD.lm(logcsize ~ age, data = gdf.max, iter = 999, RRPP=TRUE)
fit.maxage1<-procD.lm(logcsize ~ age + morph, data = gdf.max, iter = 999, RRPP=TRUE)
fit.maxage2<-procD.lm(logcsize ~ age*morph, data = gdf.max, iter = 999, RRPP=TRUE)

anova(fit.maxage0,fit.maxage1)
anova(fit.maxage1, fit.maxage2)
summary(fit.maxage2)


############################  Articular ##########################
###########################################################
## Relationship between Articular centroid size and FL.

#Testing morph effects and allometry
fit.artsize0<-procD.lm(logcsize ~ loglength, data = gdf.art, iter = 999, RRPP=TRUE)
fit.artsize1<-procD.lm(logcsize ~ loglength + morph, data = gdf.art, iter = 999, RRPP=TRUE)
fit.artsize2<-procD.lm(logcsize ~ loglength*morph, data = gdf.art, iter = 999, RRPP=TRUE)
fit.artsize3<-procD.lm(logcsize ~ loglength*morph+age, data = gdf.art, iter = 999, RRPP=TRUE)
fit.artsize4<-procD.lm(logcsize ~ loglength*morph*age, data = gdf.art, iter = 999, RRPP=TRUE)

anova(fit.artsize0,fit.artsize1)
anova(fit.artsize1,fit.artsize2)

gp.art<-interaction(gdf.art$morph)

PW2<-pairwise(fit.artsize1, fit.artsize0,groups = gp.art, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
plotAllometry(fit.artsize2,gdf.art$loglength,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.art$morph)])

plotAllometry(fit.artsize2,gdf.art$loglength,logsz = FALSE, method = c("RegScore"),pch = 19, col = palettething[as.numeric(gdf.art$morph)])


PW<-pairwise(fit.artsize2, fit.artsize1,groups = gdf.art$morph, covariate = gdf.art$loglength)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
PW<-pairwise(fit.artsize2, fit.artsize1,groups = gdf.art$morph, covariate = gdf.art$loglength)
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

### Testing correlation
cor.test(gdf.art$logcsize, gdf.art$loglength, method=c("pearson"))
plot(gdf.art$logcsize, gdf.art$loglength)

datart<-data.frame(gdf.art$logcsize, gdf.art$loglength, gdf.art$morph)
names(datart)[1] <- "logcsize"
names(datart)[2] <- "loglength"
names(datart)[3] <- "morph"

m6<-ggplot(data = datart, aes(x=loglength, y=logcsize, color=morph)) + geom_point() + labs(x="Log Fork length", y="Log CSart", color="Morph")
m6 + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))

lm3<-lm(logcsize~loglength*morph, dat=datart)
anova(lm3)
summary(lm3)

m.lst <- lstrends(lm3, "morph", var="loglength")
pairs(m.lst)

#### Effect of age
fit.artage0<-procD.lm(logcsize ~ age, data = gdf.art, iter = 999, RRPP=TRUE)
fit.artage1<-procD.lm(logcsize ~ age + morph, data = gdf.art, iter = 999, RRPP=TRUE)
fit.artage2<-procD.lm(logcsize ~ age*morph, data = gdf.art, iter = 999, RRPP=TRUE)

anova(fit.artage0,fit.artage1)
anova(fit.artage1, fit.artage2)

####################### Quadrate #####################################
####################################################################
## Relationship between Quadrate centroid size and FL.

#####Testing morph effects and allometry
fit.qusize0<-procD.lm(logcsize ~ loglength, data = gdf.qu, iter = 999, RRPP=TRUE)
fit.qusize1<-procD.lm(logcsize ~ loglength + morph, data = gdf.qu, iter = 999, RRPP=TRUE)
fit.qusize2<-procD.lm(logcsize ~ loglength*morph, data = gdf.qu, iter = 999, RRPP=TRUE)
fit.qusize3<-procD.lm(logcsize ~ loglength*morph+age, data = gdf.qu, iter = 999, RRPP=TRUE)
fit.qusize4<-procD.lm(logcsize ~ loglength*morph*age, data = gdf.qu, iter = 999, RRPP=TRUE)

anova(fit.qusize0,fit.qusize1)
anova(fit.qusize1,fit.qusize2)

gp.qu<-interaction(gdf.qu$morph)
PW2<-pairwise(fit.qusize1, fit.qusize0,groups = gp.qu, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
plotAllometry(fit.qusize2,gdf.qu$loglength,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.qu$morph)])

plotAllometry(fit.qusize2,gdf.qu$loglength,logsz = FALSE, method = c("RegScore"),pch = 19, col = palettething[as.numeric(gdf.qu$morph)])

PW<-pairwise(fit.qusize2, fit.qusize1,groups = gdf.qu$morph, covariate = gdf.qu$loglength)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
PW<-pairwise(fit.qusize2, fit.qusize1,groups = gdf.qu$morph, covariate = gdf.qu$loglength)
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

#### Testing correlation
cor.test(gdf.qu$logcsize, gdf.qu$loglength, method=c("pearson"))
plot(gdf.qu$logcsize, gdf.qu$loglength)

datqu<-data.frame(gdf.qu$logcsize, gdf.qu$loglength, gdf.qu$morph)
names(datqu)[1] <- "logcsize"
names(datqu)[2] <- "loglength"
names(datqu)[3] <- "morph"
m7<-ggplot(data = datqu, aes(x=loglength, y=logcsize, color=morph)) + geom_point() + labs(x="Log Fork length", y="Log CSqu", color="Morph")
m7 + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))

lm4<-lm(logcsize~loglength*morph, dat=datqu)
anova(lm4)
summary(lm4)

m.lst <- lstrends(lm4, "morph", var="loglength")
pairs(m.lst)

### Effect of age
fit.quage0<-procD.lm(logcsize ~ age, data = gdf.qu, iter = 999, RRPP=TRUE)
fit.quage1<-procD.lm(logcsize ~ age + morph, data = gdf.qu, iter = 999, RRPP=TRUE)
fit.quage2<-procD.lm(logcsize ~ age*morph, data = gdf.qu, iter = 999, RRPP=TRUE)

anova(fit.quage0,fit.quage1)
anova(fit.quage1, fit.quage2)

######################## Premaxilla ################################
#######################################################################
## Relationship between Premaxilla centroid size and FL.

##Testing morph effect and allometry.
fit.presize0<-procD.lm(logcsize ~ loglength, data = gdf.pre, iter = 999, RRPP=TRUE)
fit.presize1<-procD.lm(logcsize ~ loglength + morph, data = gdf.pre, iter = 999, RRPP=TRUE)
fit.presize2<-procD.lm(logcsize ~ loglength*morph, data = gdf.pre, iter = 999, RRPP=TRUE)
fit.presize3<-procD.lm(logcsize ~ loglength*morph+age, data = gdf.pre, iter = 999, RRPP=TRUE)
fit.presize4<-procD.lm(logcsize ~ loglength*morph*age, data = gdf.pre, iter = 999, RRPP=TRUE)

anova(fit.presize0,fit.presize1)
anova(fit.presize1,fit.presize2)

gp.pre<-interaction(gdf.pre$morph)
PW2<-pairwise(fit.presize1, fit.presize0,groups = gp.pre, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
plotAllometry(fit.presize2,gdf.pre$loglength,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.pre$morph)])

plotAllometry(fit.presize2,gdf.pre$loglength,logsz = FALSE, method = c("RegScore"),pch = 19, col = palettething[as.numeric(gdf.pre$morph)])

PW<-pairwise(fit.presize2, fit.presize1,groups = gdf.pre$morph, covariate = gdf.pre$loglength)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
PW<-pairwise(fit.presize2, fit.presize1,groups = gdf.pre$morph, covariate = gdf.pre$loglength)
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

### testing correlation
cor.test(gdf.pre$logcsize, gdf.pre$loglength, method=c("pearson"))
plot(gdf.pre$logcsize, gdf.pre$loglength)

datpre<-data.frame(gdf.pre$logcsize, gdf.pre$loglength, gdf.pre$morph)
names(datpre)[1] <- "logcsize"
names(datpre)[2] <- "loglength"
names(datpre)[3] <- "morph"
m8<-ggplot(data = datpre, aes(x=loglength, y=logcsize, color=morph)) + geom_point() + labs(x="Log Fork length", y="Log CSpre", color="Morph")
m8 + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))

lm5<-lm(logcsize~loglength*morph, dat=datpre)
anova(lm5)
summary(lm5)

m.lst <- lstrends(lm5, "morph", var="loglength")
pairs(m.lst)

### Effect of age
fit.preage0<-procD.lm(logcsize ~ age, data = gdf.pre, iter = 999, RRPP=TRUE)
fit.preage1<-procD.lm(logcsize ~ age + morph, data = gdf.pre, iter = 999, RRPP=TRUE)
fit.preage2<-procD.lm(logcsize ~ age*morph, data = gdf.pre, iter = 999, RRPP=TRUE)

anova(fit.preage0,fit.preage1)
anova(fit.preage1, fit.preage2)

########################## Supramaxilla #############################
####################################################################
## Relationship between Supramaxilla centroid size and FL.

#Testing morph effects and allometry
fit.supsize0<-procD.lm(logcsize ~ loglength, data = gdf.sup, iter = 999, RRPP=TRUE)
fit.supsize1<-procD.lm(logcsize ~ loglength + morph, data = gdf.sup, iter = 999, RRPP=TRUE)
fit.supsize2<-procD.lm(logcsize ~ loglength*morph, data = gdf.sup, iter = 999, RRPP=TRUE)
fit.supsize3<-procD.lm(logcsize ~ loglength*morph+age, data = gdf.sup, iter = 999, RRPP=TRUE)
fit.supsize4<-procD.lm(logcsize ~ loglength*morph*age, data = gdf.sup, iter = 999, RRPP=TRUE)

anova(fit.supsize0,fit.supsize1)
anova(fit.supsize1,fit.supsize2)

gp.sup<-interaction(gdf.sup$morph)
PW2<-pairwise(fit.supsize1, fit.supsize0,groups = gp.sup, covariate = NULL)
summary(PW2, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE) #difference between group means

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB
plotAllometry(fit.supsize2,gdf.sup$loglength,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(gdf.sup$morph)])

plotAllometry(fit.supsize2,gdf.sup$loglength,logsz = FALSE, method = c("RegScore"),pch = 19, col = palettething[as.numeric(gdf.sup$morph)])

PW<-pairwise(fit.supsize2, fit.supsize1,groups = gdf.sup$morph, covariate = gdf.sup$loglength)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
PW<-pairwise(fit.supsize2, fit.supsize1,groups = gdf.sup$morph, covariate = gdf.sup$loglength)
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")

#### Correlation test
cor.test(gdf.sup$logcsize, gdf.sup$loglength, method=c("pearson"))
plot(gdf.sup$logcsize, gdf.sup$loglength)

datsup<-data.frame(gdf.sup$logcsize, gdf.sup$loglength, gdf.sup$morph, gdf.sup$age)
names(datsup)[1] <- "logcsize"
names(datsup)[2] <- "loglength"
names(datsup)[3] <- "morph"
names(datsup)[4] <- "age"
m8<-ggplot(data = datsup, aes(x=loglength, y=logcsize, color=morph)) + geom_point() + labs(x="Log Fork length", y="Log CSsup", color="Morph")
m8 + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695"))

lm6<-lm(logcsize ~ loglength*morph, dat=datsup)
anova(lm6)
summary(lm6)

m.lst<-lstrends(lm6, "morph", var= "loglength")
pairs(m.lst)

datsupnew<-datsup[!is.na(datsup$age),]

### Effect of age
fit.supage0<-procD.lm(logcsize ~ age, data = gdf.sup, iter = 999, RRPP=TRUE)
fit.supage1<-procD.lm(logcsize ~ age + morph, data = gdf.sup, iter = 999, RRPP=TRUE)
fit.supage2<-procD.lm(logcsize ~ age*morph, data = gdf.sup, iter = 999, RRPP=TRUE)

anova(fit.supage0,fit.supage1)
anova(fit.supage1, fit.supage2)

palettething=c("#74add1", #LB
               "#f46d43", #PI
               "#a50026", #PL
               "#313695") #SB

plotAllometry(fit.supage2,datsupnew$age,logsz = FALSE, method = c("PredLine"),pch = 19, col = palettething[as.numeric(datsupnew$morph)])
plotAllometry(fit.supage2,datsupnew$age,logsz = FALSE, method = c("RegScore"),pch = 19, col = palettething[as.numeric(datsupnew$morph)])

PW<-pairwise(fit.supage1, fit.supage2,groups = datsupnew$morph, covariate = datsupnew$age)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE)
summary(PW, test.type = "VC", confidence = 0.95, stat.table = TRUE, show.vectors = FALSE, angle.type ="deg")








