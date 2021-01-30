# MP Data

library (tidyverse)
library (ggfortify)
library(ggplot2)
library(emmeans)
library(stringr)
library(latex2exp)
library (lme4)
library(lmerTest)
library(car)
library(multcomp)
library(multcompView)
library(plyr)
library(dplyr)
library(lmodel2)
library(knitr)
library (corrplot)
library(ggpmisc)


setwd ("C:/Users/linds/Desktop/MPdata")
SH=read.csv("SHData.csv")
water = read.csv("WaterColumnAnova_MinusBlanks.csv") #This data has average #blank contamination removed

#ARCHIVED DATA - OLD VERSIONS
  #water = read.csv("WaterColumnAnova.csv")
  #MPperm3 = read.csv("NetAverageMatrix.csv")
  #OLD=read.csv("MPLabData.csv")
  #MP = read.csv("Truncated.csv")

head(water)
summary(water)
str(water)

###################################################
        #Standardizing Water Column Data#
                #Done in excel#

#Standardization process:
  #Flowmeter velocity (m/sec)
#Net opening area (width of net* water height in net) (m2)
        #Maximum net height = 0.305m,
        #Net width = 0.457m.
#Flow: m3/sec = net opening area * flowmeter velocity
#Flow: m3/min = net opening area * flowmeter velocity * 60
#MP/ min: total # MP / 30 min
#mp/m3: (MP/min)  / Flow


##############################################################
    #Removing blank contamination from Water Colmn Data#
              #Done in Excel#

#Average # MP for field blanks:
    #Dark Blue Fiber: 2
    #Black fiber: .66
    #Light blue fiber: .33
    #Gray fiber: .33
    #Red fiber: .33
#Average # MP for lab blanks:
    #Dark blue fiber: 3.33
    #Black fiber: 1
    #Clear fiber: .33

#For all water column samples, the following have been deleted from each sample:
    #Dark blue fiber: 5.33
    #Black Fiber: 1.66
    #Light blue fiber: .33
    #Gray fiber: .33
    #Red fiber: .33
    #Clear fiber: .33

#"WaterColumnAnova_MinusBlanks.csv" has minimum value = 0
#"WaterColumnAnova_MinusBlanks_Negatives.csv" has minimum values <0, including negative numbers


####################################################
            #NMDS PLOTS#
#Step 1: make two different dataframes 

wateravg = read.csv("NMDSData_MinusBlanks.csv")

data1 <- water[,7:35]
data1a <- water[,44]
cbind(data1, data1a)#just sample sites with mp data
data2 <- water[,3:6]# the rest of the data


data1 <- wateravg[,17:21]
data2 <- wateravg[,3:6]# the rest of the data

#Step 2: NMDS
library(vegan)

###PLOTS###
NMDS <- metaMDS(data1, distance = "bray", k=2)
#Data Visualization
co=c("red", "blue", "green", "orange", "black")  #Use different color for each site
shape = c(15, 18, 16, 17, 19) #Use different shape for each site
windows()
plot(NMDS$points, col=co[data2$Plot], pch=shape[data2$Plot],
     cex=1.2, main="MP Concentration Among Sites", xlab="axis 1", ylab="axis 2")
ordihull(NMDS,groups=data2$Plot,draw="polygon",col=c("red", "blue", "green","orange","gray90"),label=TRUE)
#ordispider(NMDS, groups = data2$Plot, label = TRUE) 
txt <- c('Site 1', 'Site 2', 'Site 3', 'Site 4', 'Site 5')#Add Legend
legend('bottomleft', txt, pch=c(15, 18, 16, 17, 19), col=c("red", "blue", "green", "orange", "black"), cex=1, bty="y")
#Bootstrapping
fit <- adonis(data1 ~ Plot, data=data2, permutations=999, method="bray")
fit
#Check homogeneity of multivariate dispersion
distance_data<- vegdist(data1)
anova(betadisper(distance_data, data2$Plot))



###MONTH###
NMDS <- metaMDS(data1, distance = "bray", k=2)
#Data Visualization
co=c("red", "blue", "green", "orange", "black")  #Use different color for each site
shape = c(15, 18, 16, 17, 19) #Use different shape for each site
windows()
plot(NMDS$points, col=co[data2$Time], pch=shape[data2$Time],
     cex=1.2, main="MP Concentration by Month", xlab="axis 1", ylab="axis 2")
ordihull(NMDS,groups=data2$Time,draw="polygon",col=c("red", "blue", "green","orange"),label=TRUE)
#ordispider(NMDS, groups = data2$Time, label = TRUE) 
txt <- c('April', 'May', 'June', 'July')#Add Legend
legend('bottomleft', txt, pch=c(15, 18, 16, 17, 19), col=c("red", "blue", "green", "orange", "black"), cex=1, bty="y")
#Bootstrapping
fit <- adonis(data1 ~ Time, data=data2, permutations=999, method="bray")
fit
#Check homogeneity of multivariate dispersion
distance_data<- vegdist(data1)
anova(betadisper(distance_data, data2$Time))



############################################################################################

                #Create Box and Whisker Plots#

# Calculate MP Concentration

#Find groups
water$Plot<-as.factor(water$Plot)
modelTime <-lm(MP_per_m3 ~ Time, data=water, family='binomial')
modelSite <- lm(MP_per_m3 ~ Plot, data=water, family='binomial')
Group1 <- emmeans (modelTime, ~Time, type='response')%>% cld(Letters = letters)
Group2 <- emmeans (modelSite, ~Plot, type='response')%>% cld(Letters = letters)
cld_dat = as.data.frame( cld(emmeans(modelTime, ~ Time),
                             Letters = letters ) ) 
cld_dat2 = as.data.frame( cld(emmeans(modelSite, ~ Plot),
                             Letters = letters ) )

#MP by Month
windows()
ggplot(water, aes(y=MP_per_m3, x=Time, color=Time)) + geom_boxplot() + ggtitle("Microplastic Concentration by Date") +
  geom_text(data = cld_dat, aes(y = 3.5, label = .group))


# MP by Plot
ggplot(water, aes(y=MP_per_m3, x=Plot, color=Plot)) + geom_boxplot() + ggtitle("Microplastic Concentration by Site") +
  geom_text(data = cld_dat2, aes(y = 3.5, label = .group))



#Flow Rate by Month
Flow = read.csv("Flow.csv")
head(Flow)
ggplot(Flow, aes(y=Flow_m3persec, x=Month, color=Month)) + geom_boxplot() + ggtitle("Average Water Volume per Month") +
  xlab("Month") + ylab("Average Flow (m3/sec)") 






####################################################################################################
                #SCATTERPLOTS#

#Scatterplot for MP Concentration per Water Volume
ggplot(water, aes(y=MP_per_min, x=flow_m3permin)) +geom_point(shape=1) +     geom_smooth(method=lm) +ggtitle("MP Concentration per Water Volume") +ylab("Total MP per minute") + xlab ("Water Volume (m3/minute)")
aggregate(MP_per_m3~Type, if.na=FALSE, data = water, FUN = mean)
aggregate(MP_per_m3~Type, data = water, FUN = sd)

#Scatterplot for MP Concentration per Fish Fork Length
fishmodel <-lm(Total~forklength, data=SH)
ggplot(SH, aes(y=Total, x=forklength)) +geom_point(shape=1) +     geom_smooth(method=lm) +ggtitle("MP Concentration by Fork Length") +ylab("Total MP") + xlab ("Fork Length (mm)") 
summary(fishmodel)$r.squared

#MP concentration over time compared across sites
ggplot(water, aes(y=MP_per_min, x=Time, color=Plot)) +geom_point(shape=1) +     geom_smooth(method=lm) +ggtitle("MP Concentration per Water Volume") +ylab("Total MP per minute") + xlab ("Water Volume (m3/minute)")






#############################################################
              #STATS#

#Site
# Compute the analysis of variance
site.aov <- aov(MP_per_m3 ~ Plot, data = water)
# Summary of the analysis
summary(site.aov)
TukeyHSD(site.aov, conf.level = 0.99)
plot(TukeyHSD(site.aov, conf.level = 0.99),las=1, col = "red")
#residuals
plot(site.aov)



#Month
# Compute the analysis of variance
time.aov <- aov(MP_per_m3 ~ Time, data = water)
# Summary of the analysis
summary(time.aov)
TukeyHSD(time.aov, conf.level = 0.99)
plot(TukeyHSD(time.aov, conf.level = 0.99),las=1, col = "red")
#residuals
plot(time.aov)






#Fish R Squared Trend
head(SH)
steelhead.lm = lm(Total ~ forklength, data=SH)
summary(steelhead.lm)$r.squared 
ggplot(SH, aes(x=forklength, y=Total)) + geom_point() +geom_smooth(method=lm)

steelheaddate.lm = lm(Total ~ Date, data=SH)
summary(steelheaddate.lm)$r.squared 
ggplot(SH, aes(x=Date, y=Total)) + geom_point() +geom_smooth(method=lm)



####################################################
        #SUBSETS OF DATA#
#To avoid spatial autocorrelation#

            #######################
          #BIG V MILL V SCOTT CREEK#
watercomparison = read.csv("WaterColumnAnova_MinusBlanks_RiverComparisons.csv") #This data has average #blank contamination removed
#Find groups
watercomparison$Plot<-as.factor(watercomparison$Plot)
modelTime <-lm(MP_per_m3 ~ Time, data=watercomparison, family='binomial')
modelSite <- lm(MP_per_m3 ~ Plot, data=watercomparison, family='binomial')
Group1 <- emmeans (modelTime, ~Time, type='response')%>% cld(Letters = letters)
Group2 <- emmeans (modelSite, ~Plot, type='response')%>% cld(Letters = letters)
cld_dat = as.data.frame( cld(emmeans(modelTime, ~ Time),
                             Letters = letters ) ) 
cld_dat2 = as.data.frame( cld(emmeans(modelSite, ~ Plot),
                              Letters = letters ) )
#MP by Month
windows()
ggplot(watercomparison, aes(y=MP_per_m3, x=Time, color=Time)) + geom_boxplot() + ggtitle("Microplastics collected by Date") +
  geom_text(data = cld_dat, aes(y = 3.5, label = .group))
# MP by Plot
ggplot(watercomparison, aes(y=MP_per_m3, x=Plot, color=Plot)) + geom_boxplot() + ggtitle("Microplastic community by Creek") +
  geom_text(data = cld_dat2, aes(y = 3.5, label = .group))

              #STATS#
        #Site#
# Compute the analysis of variance
site.aovcomp <- aov(MP_per_m3 ~ Plot, data = watercomparison)
# Summary of the analysis
summary(site.aovcomp)
TukeyHSD(site.aovcomp, conf.level = 0.99)
plot(TukeyHSD(site.aovcomp, conf.level = 0.99),las=1, col = "red")
#residuals
plot(site.aovcomp)

        #Month#
# Compute the analysis of variance
time.aovcomp <- aov(MP_per_m3 ~ Time, data = watercomparison)
# Summary of the analysis
summary(time.aovcomp)
TukeyHSD(time.aovcomp, conf.level = 0.99)
plot(TukeyHSD(time.aovcomp, conf.level = 0.99),las=1, col = "red")
#residuals
plot(time.aovcomp)


head(wateravg)
wateravg$Plot<-as.factor(wateravg$Plot)
ggplot(wateravg, aes(y=Fiber_Concentration, x=Plot, color=Time)) + geom_boxplot() + ggtitle("Fibers present in Samples") 









                  #############################
            #COMPARISON OF PLOTS ALONG THE MAINSTEM ONLY#
mainstem = read.csv("WaterColumnAnova_MinusBlanks_MainstemOnly.csv") #This data has average #blank contamination removed
#Find groups
mainstem$Plot<-as.factor(mainstem$Plot)
modelTime <-lm(MP_per_m3 ~ Time, data=mainstem, family='binomial')
modelSite <- lm(MP_per_m3 ~ Plot, data=mainstem, family='binomial')
Group1 <- emmeans (modelTime, ~Time, type='response')%>% cld(Letters = letters)
Group2 <- emmeans (modelSite, ~Plot, type='response')%>% cld(Letters = letters)
cld_dat = as.data.frame( cld(emmeans(modelTime, ~ Time),
                             Letters = letters ) ) 
cld_dat2 = as.data.frame( cld(emmeans(modelSite, ~ Plot),
                              Letters = letters ) )
#MP by Month
windows()
ggplot(mainstem, aes(y=MP_per_m3, x=Time, color=Time)) + geom_boxplot() + ggtitle("Microplastics collected by Date") +
  geom_text(data = cld_dat, aes(y = 3.5, label = .group))
# MP by Plot
ggplot(mainstem, aes(y=MP_per_m3, x=Plot, color=Plot)) + geom_boxplot() + ggtitle("Microplastic Community along Scott Creek") +
  geom_text(data = cld_dat2, aes(y = 3.5, label = .group))

                  #STATS#
        #Site#
# Compute the analysis of variance
site.aovscott <- aov(MP_per_m3 ~ Plot, data = mainstem)
# Summary of the analysis
summary(site.aovscott)
TukeyHSD(site.aovscott, conf.level = 0.99)
plot(TukeyHSD(site.aovscott, conf.level = 0.99),las=1, col = "red")
#residuals
plot(site.aovscott)

      #Month#
# Compute the analysis of variance
time.aovscott <- aov(MP_per_m3 ~ Time, data = mainstem)
# Summary of the analysis
summary(time.aovscott)
TukeyHSD(time.aovscott, conf.level = 0.99)
plot(TukeyHSD(time.aovscott, conf.level = 0.99),las=1, col = "red")
#residuals
plot(time.aovscott)




#########################################################
      #Plot MP Types over Sites and Months#
test=read.csv("test.csv")
require(tidyr)
head(test)
data_long <- gather(test, MP, Value, Fiber_Concentration:Styrofoam_Concentration, factor_key=TRUE)
head(data_long)

data_long$Plot<-as.factor(data_long$Plot)
ggplot(data_long, aes(Value, MP, fill=factor(Plot))) +
  geom_boxplot()
ggplot(data_long, aes(Value, MP, fill=factor(Time))) +
  geom_boxplot()




##########################################################################
