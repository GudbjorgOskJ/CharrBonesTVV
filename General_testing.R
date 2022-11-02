
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
#gogn<- read.csv("Master_sheetGOJ.csv", stringsAsFactors = TRUE, header=TRUE)

##If excel
gogn<-read_excel("Master_sheetGOJ.xlsx", na="NA")
gogn$date_net_set<-as.factor(gogn$date_net_set)
gogn$date_net_up<-as.factor(gogn$date_net_up)
gogn$morph<-as.factor(gogn$morph)
gogn$sex<-as.factor(gogn$sex)
#gogn$age<-as.numeric(gogn$age)

# Making a subset for each morph
LB<-subset(gogn, morph=="LB")
PI<-subset(gogn, morph=="PI")
PL<-subset(gogn, morph=="PL")
SB<-subset(gogn, morph=="SB")

### To look for outliers
ggplot(gogn, aes(x=log(length_cm), y=log(weight_g), color=morph)) +
  geom_point(size=3) +
  theme_ipsum() +
  geom_text(label=gogn$`id-me`,size=2,hjust=0, vjust=0)

#### Data examination 
count(gogn)

table(gogn$morph)
table(gogn$sex)
table(gogn$morph, gogn$sex)

####### Effect of fishing day #############
### Results for this part not shown in paper.
### Should we just remove this part of the script for publishing

table(gogn$morph, gogn$date_net_set)

## Testing difference in length between days.
testd1<-aov(length_cm ~ date_net_set*morph*sex, data = gogn)
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

testdw<-aov(weight_g ~ date_net_set*sex*morph, data = gogn)
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
testda<-aov(age ~ date_net_set*sex*morph, data = gogn)
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

##########################################################################
####### Length, weight and age differences in the population ############
#########################################################################
###### All fishes together ######

##Length 
mean(gogn$length_cm)
sd(gogn$length_cm)
median(gogn$length_cm)

## Weight
mean(gogn$weight_g)
sd(gogn$weight_g)
median(gogn$weight_g)

## Age 
mean(gogn$age, na.rm = T)
sd(gogn$age, na.rm = T)
median(gogn$age, na.rm = T)

####### Morph differences #########

#Mean length 

gogn %>% 
  dplyr::group_by(morph) %>% 
  dplyr::summarise(mean=mean(length_cm),sd=sd(length_cm), median = median(length_cm), count = n())

ggplot(data=gogn, aes(x=morph, y=length_cm)) + geom_boxplot()

### Testing for differences between morphs
test1<-aov(length_cm ~ morph, data = gogn)
anova(test1)
TukeyHSD(test1)


#Mean Log Weight

gogn %>% 
  dplyr::group_by(morph) %>% 
  dplyr::summarise(mean=mean(weight_g),sd=sd(weight_g), median = median(weight_g), count = n())


ggplot(data=gogn, aes(x=morph, y=weight_g)) + geom_boxplot()

### Testing for differences between morphs
testw1<-aov(log(weight_g) ~ morph, data = gogn)
anova(testw1)
TukeyHSD(testw1)


# Mean age

gogn %>% 
  dplyr::group_by(morph) %>% 
  dplyr::summarise(mean=mean(age, na.rm = T),sd=sd(age, na.rm = T), median = median(age, na.rm = T), count = n())

ggplot(data=gogn, aes(x=morph, y=age)) + geom_boxplot()

### Testing for differences between morphs
testa1<-aov(age ~ morph, data = gogn)
anova(testa1)
TukeyHSD(testa1)


######## Sexual dimorphism ##########

## Length 

ggplot(data=gogn, aes(x=sex, y=length_cm)) + geom_boxplot()

gogn %>% 
  dplyr::group_by(sex) %>% 
  dplyr::summarise(mean=mean(length_cm),sd=sd(length_cm), median = median(length_cm), count = n())

### Testing for differences between the sexes (all morphs together)
test2<-aov(length_cm ~ sex, data = gogn)
anova(test2)
TukeyHSD(test2)


## Mean Log Weight
ggplot(data=gogn, aes(x=sex, y=weight_g)) + geom_boxplot()

gogn %>% 
  dplyr::group_by(sex) %>% 
  dplyr::summarise(mean=mean(weight_g),sd=sd(weight_g), median = median(weight_g), count = n())

### Testing for differences between the sexes (all morphs together)
testw2<-aov(log(weight_g) ~ sex, data = gogn)
anova(testw2)
TukeyHSD(testw2)


## Mean age

gogn %>% 
  dplyr::group_by(sex) %>% 
  dplyr::summarise(mean=mean(age, na.rm = T),sd=sd(age, na.rm = T), median = median(age, na.rm = T), count = n())

ggplot(data=gogn, aes(x=sex, y=age)) + geom_boxplot()

### Testing for differences between morphs
testa2<-aov(age ~ sex, data = gogn)
anova(testa2)
TukeyHSD(testa2)



######## Sexual dimorphism within morphs ########

# Length
gogn %>% 
  dplyr::group_by(morph, sex) %>% 
  dplyr::summarise(mean=mean(length_cm),sd=sd(length_cm), median = median(length_cm), count = n())

ggplot(data=gogn, aes(x=morph, y=length_cm, fill=sex)) + geom_boxplot()

### Testing for differences between the sexes
test3<-aov(length_cm ~ morph+sex, data = gogn)
test4<-aov(length_cm ~ morph*sex, data = gogn)
#anova(test1, test3) ### Not possible due to missing data
anova(test2, test3)
anova(test3,test4)
#Test 3 (morph + sex), best test
summary(test3)
TukeyHSD(test3)

summary(test4)
TukeyHSD(test4)



# Weight
gogn %>% 
  dplyr::group_by(morph, sex) %>% 
  dplyr::summarise(mean=mean(weight_g),sd=sd(weight_g), median = median(weight_g), count = n())

ggplot(data=gogn, aes(x=morph, y=weight_g, fill=sex)) + geom_boxplot()

### Testing for differences between the sexes
testw3<-aov(log(weight_g) ~ morph+sex, data = gogn)
testw4<-aov(log(weight_g) ~ morph*sex, data = gogn)
##anova(testw1, testw3) ## Not possible due to missing data
anova(testw2, testw3)
anova(testw3, testw4)
#Test 4 (morph * sex), best test

summary(testw3)
TukeyHSD(testw3)
summary(testw4)
TukeyHSD(testw4)

# Mean age

gogn %>% 
  dplyr::group_by(morph, sex) %>% 
  dplyr::summarise(mean=mean(age, na.rm = T),sd=sd(age, na.rm = T), median = median(age, na.rm = T), count = n())

ggplot(data=gogn, aes(x=morph, y=age, fill=sex)) + geom_boxplot()

### Testing for differences between the sexes
testa3<-aov(age ~ morph+sex, data = gogn)
testa4<-aov(age ~ morph*sex, data = gogn)
##anova(testa1,testa3) ###Not possible due to missing data
anova(testa2,testa3)
anova(testa3,testa4)
#Test 4 (morph*sex), best test

summary(testa3)
TukeyHSD(testa3)
summary(testa4)
TukeyHSD(testa4)


#########################################################################
################### Images preparation ##################################

###By morph###
m3<-ggplot(gogn, aes(x=age, fill=morph, col=morph)) +  geom_bar(alpha=0.5, position="identity")+ labs(x="Age", y="Count", fill="Morph")+  guides(fill = "none")+ facet_grid(morph ~ .) + guides(col = "none") + scale_color_manual(values=c("#00BA38","magenta4","#F8766D","#00BFC4")) + scale_fill_manual(values = c("#00BA38","magenta4","#F8766D","#00BFC4")) + theme_bw()
ml<-ggplot(gogn, aes(y= length_cm, x = morph, fill = morph)) + geom_boxplot() + labs(y="Length (cm)", x="Morph", fill="Morph") + guides(fill = "none")+ scale_fill_manual(values=c("#00BA38","magenta4","#F8766D","#00BFC4")) + coord_flip() + geom_jitter(width=0.1,alpha=0.3)+ theme_bw()+ scale_x_discrete(limits=c("SB", "PL", "PI", "LB"))
mw<-ggplot(gogn, aes(y= log(weight_g), x = morph, fill = morph)) + geom_boxplot() + labs(y="Log Weight (g)", x="Morph", fill="Morph")+ guides(fill = "none") + scale_fill_manual(values=c("#00BA38","magenta4","#F8766D","#00BFC4")) + coord_flip() + geom_jitter(width=0.1,alpha=0.3)+ theme_bw()+ scale_x_discrete(limits=c("SB", "PL", "PI", "LB"))

png("6ages.png")
multiplot(m3,ml,mw, layout = matrix(c(1,2,3) ,nrow = 1,byrow = T))
dev.off()

## By morph by sex ###
mmsa<-ggplot(gogn, aes(x=age, fill=morph, col=morph)) +  geom_bar(alpha=0.5, position="identity")+ labs(x="Age", y="Count", fill="Morph")+  guides(fill = "none")+ facet_grid(morph ~ sex) + guides(col = "none") + scale_color_manual(values=c("#00BA38","magenta4","#F8766D","#00BFC4")) + scale_fill_manual(values = c("#00BA38","magenta4","#F8766D","#00BFC4")) + theme_bw()
mmsl<-ggplot(gogn, aes(y= length_cm, x = morph, fill = sex)) + geom_boxplot() + labs(y="Length (cm)", x="Morph", fill="Sex") + guides(fill = "none")+ scale_fill_manual(values=c("#FF0066", "#3300FF")) + coord_flip() + geom_jitter(width=0.1,alpha=0.3)+ theme_bw()+ scale_x_discrete(limits=c("SB", "PL", "PI", "LB"))
mmsw<-ggplot(gogn, aes(y= log(weight_g), x = morph, fill = sex)) + geom_boxplot() + labs(y="Log Weight (g)", x="Morph", fill="Sex")+ guides(fill = "none") + scale_fill_manual(values=c("#FF0066", "#3300FF")) + coord_flip() + geom_jitter(width=0.1,alpha=0.3)+ theme_bw()+ scale_x_discrete(limits=c("SB", "PL", "PI", "LB"))

png("morphsex.png")
multiplot(mmsa,mmsl,mmsw, layout = matrix(c(1,2,3) ,nrow = 1,byrow = T))
dev.off()





