
###### General examination of external data ##########
##Library packages 
library(readxl)
library(ggplot2)
library(Rmisc)
library(knitr)
library(patchwork)
library(dplyr)
library(hrbrthemes)


##Reading data
#If csv file
datallfish<- read.csv("Master_sheet.csv", stringsAsFactors = TRUE, header=TRUE)


# Making a subset for each morph
LB<-subset(datallfish, morph=="LB")
PI<-subset(datallfish, morph=="PI")
PL<-subset(datallfish, morph=="PL")
SB<-subset(datallfish, morph=="SB")



### To look for outliers
ggplot(datallfish, aes(x=log(length_cm), y=log(weight_g), color=morph)) +
  geom_point(size=3) +
  theme_ipsum() +
  geom_text(label=datallfish$id_me,size=2,hjust=0, vjust=0)

#### Data examination 
count(datallfish)

table(datallfish$morph)
table(datallfish$sex)
table(datallfish$morph, datallfish$sex)

####### Effect of fishing day #############
### Results from this part not shown in paper

table(datallfish$morph, datallfish$date_net_set)

## Testing difference in length between days.
testd1<-aov(length_cm ~ date_net_set*morph*sex, data = datallfish)
anova(testd1)
##### Difference was expect, need to test for each morph on its own.
testdLB<-aov(length_cm ~ date_net_set, data = LB)
anova(testdLB)
TukeyHSD(testdLB)

testdSB<-aov(length_cm ~ date_net_set, data = SB)
anova(testdSB)
TukeyHSD(testdSB)

testdPL<-aov(length_cm ~ date_net_set, data = PL)
anova(testdPL)
TukeyHSD(testdPL)

testdPI<-aov(length_cm ~ date_net_set, data = PI)
anova(testdPI)
TukeyHSD(testdPI)

## Testing difference in weight between days.

testdw<-aov(weight_g ~ date_net_set*sex*morph, data = datallfish)
anova(testdw)

testdwLB<-aov(weight_g ~ date_net_set, data = LB)
anova(testdwLB)
TukeyHSD(testdwLB)

testdwSB<-aov(weight_g ~ date_net_set, data = SB)
anova(testdwSB)
TukeyHSD(testdwSB)

testdwPL<-aov(weight_g ~ date_net_set, data = PL)
anova(testdwPL)
TukeyHSD(testdwPL)

testdwPI<-aov(weight_g ~ date_net_set, data = PI)
anova(testdwPI)
TukeyHSD(testdwPI)

## Testing difference in age between days.
testda<-aov(age ~ date_net_set*sex*morph, data = datallfish)
anova(testda)

testdaLB<-aov(age ~ date_net_set, data = LB)
anova(testdaLB)
TukeyHSD(testdaLB)

testdaSB<-aov(age ~ date_net_set, data = SB)
anova(testdaSB)
TukeyHSD(testdaSB)

testdaPL<-aov(age ~ date_net_set, data = PL)
anova(testdPL)
TukeyHSD(testdaPL)

testdaPI<-aov(age ~ date_net_set, data = PI)
anova(testdaPI)
TukeyHSD(testdaPI)

########################################################################
####### Length, weight and age differences in the population ############
######################################################################
###### All fishes together ######################

##Length 
mean(datallfish$length_cm)
sd(datallfish$length_cm)
median(datallfish$length_cm)

## Weight
mean(datallfish$weight_g)
sd(datallfish$weight_g)
median(datallfish$weight_g)

## Age 
mean(datallfish$age, na.rm = T)
sd(datallfish$age, na.rm = T)
median(datallfish$age, na.rm = T)

####### Morph differences #########

#Mean length 

datallfish %>% 
  dplyr::group_by(morph) %>% 
  dplyr::summarise(mean=mean(length_cm),sd=sd(length_cm), median = median(length_cm), count = n())

ggplot(data=datallfish, aes(x=morph, y=length_cm)) + geom_boxplot()

### Testing for differences between morphs
testl1<-aov(length_cm ~ morph, data = datallfish)
anova(testl1)
TukeyHSD(testl1)


#Mean Weight

datallfish %>% 
  dplyr::group_by(morph) %>% 
  dplyr::summarise(mean=mean(weight_g),sd=sd(weight_g), median = median(weight_g), count = n())


ggplot(data=datallfish, aes(x=morph, y=weight_g)) + geom_boxplot()

### Testing for differences between morphs
testw1<-aov(log(weight_g) ~ morph, data = datallfish)
anova(testw1)
TukeyHSD(testw1)


# Mean age

datallfish %>% 
  dplyr::group_by(morph) %>% 
  dplyr::summarise(mean=mean(age, na.rm = T),sd=sd(age, na.rm = T), median = median(age, na.rm = T), count = n())

ggplot(data=datallfish, aes(x=morph, y=age)) + geom_boxplot()

### Testing for differences between morphs
testa1<-aov(age ~ morph, data = datallfish)
anova(testa1)
TukeyHSD(testa1)


######## Sexual dimorphism ##########

## Length 

ggplot(data=datallfish, aes(x=sex, y=length_cm)) + geom_boxplot()

datallfish %>% 
  dplyr::group_by(sex) %>% 
  dplyr::summarise(mean=mean(length_cm),sd=sd(length_cm), median = median(length_cm), count = n())

### Testing for differences between the sexes (all morphs together)
testl2<-aov(length_cm ~ sex, data = datallfish)
anova(testl2)
TukeyHSD(testl2)


#Mean Weight
ggplot(data=datallfish, aes(x=sex, y=weight_g)) + geom_boxplot()

datallfish %>% 
  dplyr::group_by(sex) %>% 
  dplyr::summarise(mean=mean(weight_g),sd=sd(weight_g), median = median(weight_g), count = n())

### Testing for differences between the sexes (all morphs together)
testw2<-aov(log(weight_g) ~ sex, data = datallfish)
anova(testw2)
TukeyHSD(testw2)


# Mean age

datallfish %>% 
  dplyr::group_by(sex) %>% 
  dplyr::summarise(mean=mean(age, na.rm = T),sd=sd(age, na.rm = T), median = median(age, na.rm = T), count = n())

ggplot(data=datallfish, aes(x=sex, y=age)) + geom_boxplot()

### Testing for differences between morphs
testa2<-aov(age ~ sex, data = datallfish)
anova(testa2)
TukeyHSD(testa2)



######## Sexual dimorphism within morphs #######

# Length
datallfish %>% 
  dplyr::group_by(morph, sex) %>% 
  dplyr::summarise(mean=mean(length_cm),sd=sd(length_cm), median = median(length_cm), count = n())

ggplot(data=datallfish, aes(x=morph, y=length_cm, fill=sex)) + geom_boxplot()

### Testing for differences between the sexes
testl3<-aov(length_cm ~ morph+sex, data = datallfish)
testl4<-aov(length_cm ~ morph*sex, data = datallfish)
#anova(testl1, testl3) ### Not possible due to missing data
anova(testl2, testl3)
anova(testl3,testl4)
#Test 3 (morph + sex), best test
summary(testl3)
TukeyHSD(testl3)

summary(testl4)
TukeyHSD(testl4)



# Weight
datallfish %>% 
  dplyr::group_by(morph, sex) %>% 
  dplyr::summarise(mean=mean(weight_g),sd=sd(weight_g), median = median(weight_g), count = n())

ggplot(data=datallfish, aes(x=morph, y=weight_g, fill=sex)) + geom_boxplot()

### Testing for differences between the sexes
testw3<-aov(log(weight_g) ~ morph+sex, data = datallfish)
testw4<-aov(log(weight_g) ~ morph*sex, data = datallfish)
##anova(testw1, testw3) ## Not possible due to missing data
anova(testw2, testw3)
anova(testw3, testw4)
#Test 4 (morph * sex), best test

summary(testw3)
TukeyHSD(testw3)
summary(testw4)
TukeyHSD(testw4)

# Mean age

datallfish %>% 
  dplyr::group_by(morph, sex) %>% 
  dplyr::summarise(mean=mean(age, na.rm = T),sd=sd(age, na.rm = T), median = median(age, na.rm = T), count = n())

ggplot(data=datallfish, aes(x=morph, y=age, fill=sex)) + geom_boxplot()

### Testing for differences between the sexes
testa3<-aov(age ~ morph+sex, data = datallfish)
testa4<-aov(age ~ morph*sex, data = datallfish)
##anova(testa1,testa3) ###Not possible due to missing data
anova(testa2,testa3)
anova(testa3,testa4)
#Test 4 (morph*sex), best test

summary(testa3)
TukeyHSD(testa3)
summary(testa4)
TukeyHSD(testa4)

#########################################################
##############  Images preparation ######################

###By morph###
m3<-ggplot(datallfish, aes(x=age, fill=morph, col=morph)) +  geom_bar(alpha=0.5, position="identity")+ labs(x="Age", y="Count", fill="Morph")+  guides(fill = "none")+ facet_grid(morph ~ .) + guides(col = "none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()
ml<-ggplot(datallfish, aes(y= length_cm, x = morph, fill = morph)) + geom_boxplot() + labs(y="Length (cm)", x="Morph", fill="Morph") + guides(fill = "none")+ scale_fill_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + coord_flip() + geom_jitter(width=0.1,alpha=0.3)+ theme_bw()+ scale_x_discrete(limits=c("SB", "PL", "PI", "LB"))
mw<-ggplot(datallfish, aes(y= log(weight_g), x = morph, fill = morph)) + geom_boxplot() + labs(y="Log Weight (g)", x="Morph", fill="Morph")+ guides(fill = "none") + scale_fill_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + coord_flip() + geom_jitter(width=0.1,alpha=0.3)+ theme_bw()+ scale_x_discrete(limits=c("SB", "PL", "PI", "LB"))

png("6ages.png")
multiplot(m3,ml,mw, layout = matrix(c(1,2,3) ,nrow = 1,byrow = T))
dev.off()

### By morph by sex ###

mmsa<-ggplot(datallfish, aes(x=age, fill=morph, col=morph)) +  geom_bar(alpha=0.5, position="identity")+ labs(x="Age", y="Count", fill="Morph")+  guides(fill = "none")+ facet_grid(morph ~ sex) + guides(col = "none") + scale_color_manual(values=c("#74add1","#f46d43","#a50026","#313695")) + scale_fill_manual(values = c("#74add1","#f46d43","#a50026","#313695")) + theme_bw()
mmsl<-ggplot(datallfish, aes(y= length_cm, x = morph, fill = sex)) + geom_boxplot() + labs(y="Length (cm)", x="Morph", fill="Sex") + guides(fill = "none")+ scale_fill_manual(values=c("#af8dc3", "#7fbf7b")) + coord_flip() + geom_jitter(width=0.1,alpha=0.3)+ theme_bw()+ scale_x_discrete(limits=c("SB", "PL", "PI", "LB"))
mmsw<-ggplot(datallfish, aes(y= log(weight_g), x = morph, fill = sex)) + geom_boxplot() + labs(y="Log Weight (g)", x="Morph", fill="Sex")+ guides(fill = "none") + scale_fill_manual(values=c("#af8dc3", "#7fbf7b")) + coord_flip() + geom_jitter(width=0.1,alpha=0.3)+ theme_bw()+ scale_x_discrete(limits=c("SB", "PL", "PI", "LB"))

png("morphsex.png")
multiplot(mmsa,mmsl,mmsw, layout = matrix(c(1,2,3) ,nrow = 1,byrow = T))
dev.off()


pdf("morphsex.pdf")
multiplot(mmsa,mmsl,mmsw, layout = matrix(c(1,2,3) ,nrow = 1,byrow = T))
dev.off()
