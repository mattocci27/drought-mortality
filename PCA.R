# Phisamai Maenpuen: pissamai@xtbg.ac.cn
# Kyle W. Tomlinson: kyle.tomlinson@xtbg.org.cn
# 2021-07-15
#  
# Code for PCA (Fig. 2)

rm(list=ls());ls()

library(FactoMineR)
library(missMDA)
library(factoextra)
library(tidyverse)
library(vegan)
library(dplyr)
library(ggpubr)

################################################################################################
################################# PCA group by Phenology #######################################
df1 <- read.csv("data/PCA_Log.csv", row.names = 1)
str(df1)
df1$SWCleaf <- as.numeric(df1$SWCleaf)

drops <- c("Abbr", "Lifeform", "Dieback", "Mortality")
df1 <- df1[ , !(names(df1) %in% drops)]
dim(df1)

names(df1)[1] <- expression("Phenology")
names(df1)[2] <- expression("MVL")
names(df1)[3] <- expression("Hv")
names(df1)[4] <- expression("K[s]")
names(df1)[5] <- expression("P[50]")
names(df1)[6] <- expression("P[88]")
names(df1)[7] <- expression(" ")
names(df1)[8] <- expression("rho[wood]")
names(df1)[9] <- expression("  ")
names(df1)[10] <- expression("   ")
names(df1)[11] <- expression("    ")
names(df1)[12] <- expression("     ")
names(df1)[13] <- expression("C[leaf]")
names(df1)[14] <- expression("epsilon")
names(df1)[15] <- expression("RWC[tlp]")
names(df1)[16] <- expression("LMA")
names(df1)[17] <- expression("rho[leaf]")
names(df1)[18] <- expression("HSM[50]")
names(df1)[19] <- expression("HSM[88]")
names(df1)[20] <- expression("HSM[tlp]")
names(df1)

nb_ALL <- estim_ncpPCA(df1[,2:20], ncp.max=5)
res_comp_Species <- imputePCA(df1[,2:20] , ncp=2, scale=T, method=c("Regularized"))
res_pca_Species <- PCA(res_comp_Species$completeObs)

write.csv(res_pca_Species$var$coord, "data/PCA_var.csv", row.names=T)
#write.csv(cbind(Traits[,1], res_pca_Species$ind$coord), "data/PCA_ind_Phenology.csv", row.names=T)

df <- facto_summarize(res_pca_Species, element = "var", 
                      result = c("coord", "contrib", "cos2"), axes = c(1,2))
colnames(df)[2:3] <-  c("x", "y")

# need this to calc scale
pca.ind <- get_pca_ind(res_pca_Species)
ind <- data.frame(pca.ind$coord[, c(1,2), drop=FALSE], stringsAsFactors = TRUE)
colnames(ind)<- c("x", "y")

# rescale variable coordinates
r <- min(
  (max(ind[,"x"])-min(ind[,"x"])/(max(df[,"x"])-min(df[,"x"]))),
  (max(ind[,"y"])-min(ind[,"y"])/(max(df[,"y"])-min(df[,"y"])))
)

# Greek letter with geom_text
fviz_pca_biplot(res_pca_Species, axes.linetype = "dotted", 
                col.var = "black", alpha.var = 0.5,
                col.ind = df1$Phenology, pointsize = 2, alpha.ind = 0.7, geom.ind = "point",
                palette = c("#FF8C00", "#008B00", "#1874CD"),
                #addEllipses = TRUE, ellipse.alpha = 0.05, ellipse.type = "convex",
                title="", arrowsize = 0.5, labelsize = 5, col.circle = "white",
                label = FALSE,
                mean.point = FALSE,
                repel = TRUE,   # Avoid text overlapping
                ggtheme=theme(legend.text = element_text(size = 10, face = "italic"),
                              legend.title = element_blank(),
                              strip.text = element_blank(),
                              axis.line=element_line(size=0.7),
                              axis.ticks.length = unit(0.15, "cm"),
                              axis.ticks=element_line(size=0.7),
                              axis.text=element_text(size = 13, face = "bold"),
                              axis.title=element_text(size=16, face="bold"),
                              strip.background = element_blank(),
                              panel.background = element_rect(fill = "white"),
                              panel.border = element_rect(fill = NA,
                                                          colour = "black",
                                                          size = 0.5),
                              panel.grid.major = element_line(colour = NA),
                              panel.grid.minor = element_line(colour = NA, size = 0.25)))+
                geom_text(data = df, 
                              aes(x = x * r * 0.85, y = y * r * 0.85, label = name), 
                              size = 5,
                              parse = TRUE, repel = TRUE) +
                labs(y= "Axis2 (17.39%)", x = "Axis1 (45.31%)") +
                coord_cartesian(xlim = c(-7.5, 7.5), ylim = c(-5.5, 5.5)) +
                scale_x_continuous(sec.axis = sec_axis(~ . * 1, labels = c("-1","-0.5","0","0.5","1"))) +
                scale_y_continuous(sec.axis = sec_axis(~ . * 1, labels = c("-1","-0.5","0","0.5","1"))) +
               theme(legend.position = c(0.15, 0.15)) 
  

ggsave("figs/PCA_Phenology3.tiff", height=5, width=6, dpi=600)
ggsave("figs/PCA_Phenology3.pdf", height=5, width=6)
dev.off()

##############################################################################################
################################ PCA group by Lifeform #######################################
df2 <- read.csv("data/PCA_Log.csv", row.names = 1)
str(df2)
df2$SWCleaf <- as.numeric(df2$SWCleaf)

drops <- c("Abbr", "Phenology", "Dieback", "Mortality")
df2 <- df2[ , !(names(df2) %in% drops)]
dim(df2)

names(df2)[1] <- expression("Lifeform")
names(df2)[2] <- expression("MVL")
names(df2)[3] <- expression("Hv")
names(df2)[4] <- expression("K[s]")
names(df2)[5] <- expression("P[50]")
names(df2)[6] <- expression("P[88]")
names(df2)[7] <- expression(" ")
names(df2)[8] <- expression("rho[wood]")
names(df2)[9] <- expression("  ")
names(df2)[10] <- expression("   ")
names(df2)[11] <- expression("    ")
names(df2)[12] <- expression("     ")
names(df2)[13] <- expression("C[leaf]")
names(df2)[14] <- expression("epsilon")
names(df2)[15] <- expression("RWC[tlp]")
names(df2)[16] <- expression("LMA")
names(df2)[17] <- expression("rho[leaf]")
names(df2)[18] <- expression("HSM[50]")
names(df2)[19] <- expression("HSM[88]")
names(df2)[20] <- expression("HSM[tlp]")
names(df2)

nb_ALL <- estim_ncpPCA(df2[,2:20], ncp.max=5)
res_comp_Species <- imputePCA(df2[,2:20] , ncp=2, scale=T, method=c("Regularized"))
res_pca_Species <- PCA(res_comp_Species$completeObs)

write.csv(cbind(df2[,1], res_pca_Species$ind$coord), "data/PCA_ind_Lifeform.csv", row.names=T)

df <- facto_summarize(res_pca_Species, element = "var", 
                      result = c("coord", "contrib", "cos2"), axes = c(1,2))
colnames(df)[2:3] <-  c("x", "y")

# need this to calc scale
pca.ind <- get_pca_ind(res_pca_Species)
ind <- data.frame(pca.ind$coord[, c(1,2), drop=FALSE], stringsAsFactors = TRUE)
colnames(ind)<- c("x", "y")

# rescale variable coordinates
r <- min(
  (max(ind[,"x"])-min(ind[,"x"])/(max(df[,"x"])-min(df[,"x"]))),
  (max(ind[,"y"])-min(ind[,"y"])/(max(df[,"y"])-min(df[,"y"])))
)

# Greek letter with geom_text
fviz_pca_biplot(res_pca_Species, axes.linetype = "dotted", 
                col.var = "black", alpha.var = 0.5,
                col.ind = df2$Lifeform, pointsize = 2, alpha.ind = 0.7, geom.ind = "point",
                palette = c("red", "dodgerblue4", "cyan"),
                #addEllipses = TRUE, ellipse.alpha = 0.05, ellipse.type = "convex",
                title="", arrowsize = 0.5, labelsize = 5, col.circle = "white",
                label = FALSE,
                mean.point = FALSE,
                repel = TRUE,   # Avoid text overlapping
                ggtheme=theme(legend.text = element_text(size = 10, face = "italic"),
                              legend.title = element_blank(),
                              strip.text = element_blank(),
                              axis.line=element_line(size=0.7),
                              axis.ticks.length = unit(0.15, "cm"),
                              axis.ticks=element_line(size=0.7),
                              axis.text=element_text(size = 13, face = "bold"),
                              axis.title=element_text(size=16, face="bold"),
                              strip.background = element_blank(),
                              panel.background = element_rect(fill = "white"),
                              panel.border = element_rect(fill = NA,
                                                          colour = "black",
                                                          size = 0.5),
                              panel.grid.major = element_line(colour = NA),
                              panel.grid.minor = element_line(colour = NA, size = 0.25)))+
  geom_text(data = df, 
            aes(x = x * r * 0.85, y = y * r * 0.85, label = name), 
            size = 5,
            parse = TRUE, repel = TRUE) +
  labs(y= "Axis2 (17.39%)", x = "Axis1 (45.31%)") +
  coord_cartesian(xlim = c(-7.5, 7.5), ylim = c(-5.5, 5.5)) +
  scale_x_continuous(sec.axis = sec_axis(~ . * 1, labels = c("-1","-0.5","0","0.5","1"))) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 1, labels = c("-1","-0.5","0","0.5","1"))) +
  theme(legend.position = c(0.1, 0.15))

ggsave("figs/PCA_Lifeform3.tiff", height=5, width=6, dpi=600)
ggsave("figs/PCA_Lifeform3.pdf", height=5, width=6)
dev.off()

#################################################################################################
###################################### LOGISTIC REGRESSION PLOTS COMBINED
df3 <- read.csv("data/PCA_Log.csv", row.names = 1)
str(df3)
df3$SWCleaf <- as.numeric(df3$SWCleaf)
df3$Lifeform <- as.factor(df3$Lifeform)
df3$Phenology <- as.factor(df3$Phenology)
summary(df3)
df3[,c(1,2)]

#############################################################
#### MAKE LIFEFORM PLOTS COMBINED
#fix colours
colours = c("Liana" = "#de2e19", "Shrub" = "#1f306f", "Tree" = "#61f0f8")
#fix symbols
shapes = c("Liana"=16, "Shrub"=17, "Tree"=15)

#subset to shrubdata only
groupdata <- df3
dim(groupdata)

#subset trait data
traits <- groupdata[,c(4:22)]  # the data for traits
summary(traits)
dim(traits)

#impute missing values
nb_ALL <- estim_ncpPCA(traits, ncp.max=5)
grouptraits <- imputePCA(traits , ncp=2, scale=T, method=c("Regularized"))

## run PCA 
res_pca_Species <- PCA(grouptraits$completeObs)

fviz_pca_biplot(res_pca_Species, repel = TRUE,col.var = "black", col.ind = "black")

#########################################
#LOGISTIC REGRESSION

#store the old data
olddata <- data.frame(groupdata$Mortality,groupdata$Dieback,groupdata$Phenology,groupdata$Lifeform,res_pca_Species$ind$coord[,c(1,2)])
names(olddata) <- c("Mortality","Dieback","Phenology","Lifeform","PCA1","PCA2")
summary(olddata)

# A. mortality
#PCA 1
glm1q <- glm(cbind(Mortality,100-Mortality)~PCA1*Lifeform,olddata,family="quasibinomial")
anova(glm1q,test='Chisq')
glm2q <- glm(cbind(Mortality,100-Mortality)~PCA1+Lifeform,olddata,family="quasibinomial")
anova(glm2q,test='Chisq')
glm3q <- glm(cbind(Mortality,100-Mortality)~PCA1,olddata,family="quasibinomial")
anova(glm3q,test='Chisq')
glm4q <- glm(cbind(Mortality,100-Mortality)~Lifeform,olddata,family="quasibinomial")
anova(glm4q,test='Chisq')

#PCA2
#second, test mortality against PCA axis2
glm1q2 <- glm(cbind(Mortality,100-Mortality)~PCA2*Lifeform,olddata,family="quasibinomial")
anova(glm1q2,test='Chisq')
glm2q2 <- glm(cbind(Mortality,100-Mortality)~PCA2+Lifeform,olddata,family="quasibinomial")
anova(glm2q2,test='Chisq')
glm3q2 <- glm(cbind(Mortality,100-Mortality)~PCA2,olddata,family="quasibinomial")
anova(glm3q2,test='Chisq')
glm4q2 <- glm(cbind(Mortality,100-Mortality)~Lifeform,olddata,family="quasibinomial")
anova(glm4q2,test='Chisq')

#make mortality plot

## build new data frame to house predictions
##it should cover the same range as the values you have recorded in oberevations
newdata <- data.frame(olddata$PCA1,olddata$Lifeform); names(newdata) <- c("PCA1","Lifeform")
newdata <- arrange(newdata,PCA1)
newdata
summary(newdata)

## make the predictions - notice the type='link' argument
preds <- predict(glm1q, newdata=newdata, se.fit=T)
#preds <- predict(glm1q, newdata=newdata,type= 'link', se.fit=T)

## add CI predictions to new data frame in the transformed range 
# as n = 50 for original data, using z=1.96 to get the 95% CIs is reasonable
newdata$fit_T <- preds$fit
newdata$upr_T <- preds$fit+preds$se.fit*1.96
newdata$lwr_T <- preds$fit-preds$se.fit*1.96

#back-transform them with plogis to get response range
newdata$fit <- plogis(newdata$fit_T)
newdata$upr <- plogis(newdata$upr_T) 
newdata$lwr <- plogis(newdata$lwr_T)

head(newdata) #note the extra columns in the dataframe

ggplot(data=olddata, aes(x=PCA1,col=Lifeform,shape=Lifeform)) + 
  geom_point(aes(y = Mortality/100), size=4) + 
  geom_smooth(data=newdata,aes(x=PCA1, y=fit, ymin=lwr, ymax=upr),stat='identity')+ 
  theme_classic(base_size=24) + 
  ylab("Top-kill ratio") +
  geom_ribbon(data=newdata,aes(ymin = lwr, ymax = upr, fill =Lifeform), 
              alpha = 0.5, col = "transparent", show.legend = FALSE)+
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+ 
  scale_shape_manual(values = shapes)

#only shrubs are significant
newdata <- subset(newdata,Lifeform=="Shrub")

## plot raw data and model prediction in the response range
#for this we will use ggplot
g1.1 <- ggplot(data=olddata, aes(x=PCA1,col=Lifeform,shape=Lifeform)) + 
  geom_point(aes(y = Mortality/100), size=3, show.legend = FALSE) + 
  geom_smooth(data=newdata,aes(x=PCA1, y=fit, ymin=lwr, ymax=upr),stat='identity', show.legend = FALSE, lwd = 1.2)+  
  theme_classic(base_size=15) + 
  ylab("Top-kill ratio") +   
  geom_ribbon(data=newdata,aes(ymin = lwr, ymax = upr, fill =Lifeform), 
              alpha = 0.3, col = "transparent", show.legend = FALSE)+ 
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+ 
  scale_shape_manual(values = shapes)
g1.1

#################################
# B. Dieback

#PCA 1

glm1q <- glm(cbind(Dieback,100-Dieback)~PCA1*Lifeform,olddata,family="quasibinomial")
anova(glm1q,test='Chisq')
glm2q <- glm(cbind(Dieback,100-Dieback)~PCA1+Lifeform,olddata,family="quasibinomial")
anova(glm2q,test='Chisq')
glm3q <- glm(cbind(Dieback,100-Dieback)~PCA1,olddata,family="quasibinomial")
anova(glm3q,test='Chisq')
glm4q <- glm(cbind(Dieback,100-Dieback)~Lifeform,olddata,family="quasibinomial")
anova(glm4q,test='Chisq')

#PCA2

#second, test mortality against PCA axis2
glm1q2 <- glm(cbind(Dieback,100-Dieback)~PCA2*Lifeform,olddata,family="quasibinomial")
anova(glm1q2,test='Chisq')
glm2q2 <- glm(cbind(Dieback,100-Dieback)~PCA2+Lifeform,olddata,family="quasibinomial")
anova(glm2q2,test='Chisq')
glm3q2 <- glm(cbind(Dieback,100-Dieback)~PCA2,olddata,family="quasibinomial")
anova(glm3q2,test='Chisq')
glm4q2 <- glm(cbind(Dieback,100-Dieback)~Lifeform,olddata,family="quasibinomial")
anova(glm4q2,test='Chisq')

#################################

#make dieback plot

## build new data frame to house predictions
##it should cover the same range as the values you have recorded in oberevations
newdata <- data.frame(olddata$PCA1,olddata$Lifeform); names(newdata) <- c("PCA1","Lifeform")
newdata <- arrange(newdata,PCA1)
newdata
summary(newdata)

## make the predictions - notice the type='link' argument
preds <- predict(glm1q, newdata=newdata, se.fit=T)
#preds <- predict(glm1q, newdata=newdata,type= 'link', se.fit=T)

## add CI predictions to new data frame in the transformed range 
# as n = 50 for original data, using z=1.96 to get the 95% CIs is reasonable
newdata$fit_T <- preds$fit
newdata$upr_T <- preds$fit+preds$se.fit*1.96
newdata$lwr_T <- preds$fit-preds$se.fit*1.96

#back-transform them with plogis to get response range
newdata$fit <- plogis(newdata$fit_T)
newdata$upr <- plogis(newdata$upr_T) 
newdata$lwr <- plogis(newdata$lwr_T)

head(newdata) #note the extra columns in the dataframe
str(newdata)

ggplot(data=olddata, aes(x=PCA1,col=Lifeform,shape=Lifeform)) + 
  geom_point(aes(y = Dieback/100), size=4, show.legend = FALSE) + 
  geom_smooth(data=newdata,aes(x=PCA1, y=fit, ymin=lwr, ymax=upr),stat='identity', show.legend = FALSE)+  
  theme_classic(base_size=24) + 
  ylab("Branch dieback ratio") +   
  geom_ribbon(data=newdata,aes(ymin = lwr, ymax = upr, fill =Lifeform), 
              alpha = 0.5, col = "transparent", show.legend = FALSE)+ 
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+ 
  scale_shape_manual(values = shapes)

## plot raw data and model prediction in the response range
#for this we will use ggplot
g1.2 <- ggplot(data=olddata, aes(x=PCA1,col=Lifeform,shape=Lifeform)) + 
  geom_point(aes(y = Dieback/100), size=3, show.legend = FALSE) + 
  geom_smooth(data=newdata,aes(x=PCA1, y=fit, ymin=lwr, ymax=upr),stat='identity', show.legend = FALSE, lwd = 1.2)+  
  theme_classic(base_size=15) + 
  ylab("Branch dieback ratio") +   
  geom_ribbon(data=newdata,aes(ymin = lwr, ymax = upr, fill =Lifeform), 
              alpha = 0.3, col = "transparent", show.legend = FALSE)+ 
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+ 
  scale_shape_manual(values = shapes)
g1.2

#############################################################
#### MAKE LEAF HABIT PLOTS COMBINED
#fix colours
colours = c("Deciduous" = "#FF8C00", "Evergreen" = "#008B00", "Semi-deciduous" = "#1874CD")
#fix symbols
shapes = c("Deciduous"=16, "Evergreen"=17, "Semi-deciduous"=15)

#subset to shrubdata only
groupdata <- df3
dim(groupdata)

#subset trait data
traits <- groupdata[,c(4:22)]  # the data for traits
summary(traits)
dim(traits)

#impute missing values
nb_ALL <- estim_ncpPCA(traits, ncp.max=5)
grouptraits <- imputePCA(traits , ncp=2, scale=T, method=c("Regularized"))

## run PCA 
res_pca_Species <- PCA(grouptraits$completeObs)

fviz_pca_biplot(res_pca_Species, repel = TRUE,col.var = "black", col.ind = "black")

#########################################
#LOGISTIC REGRESSION

#store the old data
olddata <- data.frame(groupdata$Mortality,groupdata$Dieback,groupdata$Phenology,groupdata$Lifeform,res_pca_Species$ind$coord[,c(1,2)])
names(olddata) <- c("Mortality","Dieback","Phenology","Lifeform","PCA1","PCA2")
summary(olddata)

# A. mortality
#PCA 1
glm1q <- glm(cbind(Mortality,100-Mortality)~PCA1*Phenology,olddata,family="quasibinomial")
anova(glm1q,test='Chisq')
glm2q <- glm(cbind(Mortality,100-Mortality)~PCA1+Phenology,olddata,family="quasibinomial")
anova(glm2q,test='Chisq')
glm3q <- glm(cbind(Mortality,100-Mortality)~PCA1,olddata,family="quasibinomial")
anova(glm3q,test='Chisq')
glm4q <- glm(cbind(Mortality,100-Mortality)~Phenology,olddata,family="quasibinomial")
anova(glm4q,test='Chisq')

#PCA2
#second, test mortality against PCA axis2
glm1q2 <- glm(cbind(Mortality,100-Mortality)~PCA2*Phenology,olddata,family="quasibinomial")
anova(glm1q2,test='Chisq')
glm2q2 <- glm(cbind(Mortality,100-Mortality)~PCA2+Phenology,olddata,family="quasibinomial")
anova(glm2q2,test='Chisq')
glm3q2 <- glm(cbind(Mortality,100-Mortality)~PCA2,olddata,family="quasibinomial")
anova(glm3q2,test='Chisq')
glm4q2 <- glm(cbind(Mortality,100-Mortality)~Phenology,olddata,family="quasibinomial")
anova(glm4q2,test='Chisq')

#make mortality plot
## build new data frame to house predictions
##it should cover the same range as the values you have recorded in oberevations
newdata <- data.frame(olddata$PCA1,olddata$Phenology); names(newdata) <- c("PCA1","Phenology")
newdata <- arrange(newdata,PCA1)
newdata
summary(newdata)

## make the predictions - notice the type='link' argument
preds <- predict(glm1q, newdata=newdata, se.fit=T)
#preds <- predict(glm1q, newdata=newdata,type= 'link', se.fit=T)

## add CI predictions to new data frame in the transformed range 
# as n = 50 for original data, using z=1.96 to get the 95% CIs is reasonable
newdata$fit_T <- preds$fit
newdata$upr_T <- preds$fit+preds$se.fit*1.96
newdata$lwr_T <- preds$fit-preds$se.fit*1.96

#back-transform them with plogis to get response range
newdata$fit <- plogis(newdata$fit_T)
newdata$upr <- plogis(newdata$upr_T) 
newdata$lwr <- plogis(newdata$lwr_T)

head(newdata) #note the extra columns in the dataframe

ggplot(data=olddata, aes(x=PCA1,col=Phenology,shape=Phenology)) + 
  geom_point(aes(y = Mortality/100), size=4, show.legend = FALSE) + 
  geom_smooth(data=newdata,aes(x=PCA1, y=fit, ymin=lwr, ymax=upr),stat='identity', show.legend = FALSE)+  
  theme_classic(base_size=24) + 
  ylab("Top-kill ratio") +   
  geom_ribbon(data=newdata,aes(ymin = lwr, ymax = upr, fill =Phenology), 
              alpha = 0.5, col = "transparent", show.legend = FALSE)+ 
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+ 
  scale_shape_manual(values = shapes)

#subset
newdata <- subset(newdata,Phenology=="Semi-deciduous")

## plot raw data and model prediction in the response range
#for this we will use ggplot
g2.1 <- ggplot(data=olddata, aes(x=PCA1,col=Phenology,shape=Phenology)) + 
  geom_point(aes(y = Mortality/100), size=3, show.legend = FALSE) + 
  geom_smooth(data=newdata,aes(x=PCA1, y=fit, ymin=lwr, ymax=upr),stat='identity', show.legend = FALSE, lwd = 1.2)+  
  theme_classic(base_size=15) + 
  ylab("Top-kill ratio") +   
  geom_ribbon(data=newdata,aes(ymin = lwr, ymax = upr, fill =Phenology), 
              alpha = 0.3, col = "transparent", show.legend = FALSE)+ 
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+ 
  scale_shape_manual(values = shapes)
g2.1

#################################
# B. Dieback
#PCA 1
glm1q <- glm(cbind(Dieback,100-Dieback)~PCA1*Phenology,olddata,family="quasibinomial")
anova(glm1q,test='Chisq')
glm2q <- glm(cbind(Dieback,100-Dieback)~PCA1+Phenology,olddata,family="quasibinomial")
anova(glm2q,test='Chisq')
glm3q <- glm(cbind(Dieback,100-Dieback)~PCA1,olddata,family="quasibinomial")
anova(glm3q,test='Chisq')
glm4q <- glm(cbind(Dieback,100-Dieback)~Phenology,olddata,family="quasibinomial")
anova(glm4q,test='Chisq')

#PCA2
#second, test mortality against PCA axis2
glm1q2 <- glm(cbind(Dieback,100-Dieback)~PCA2*Phenology,olddata,family="quasibinomial")
anova(glm1q2,test='Chisq')
glm2q2 <- glm(cbind(Dieback,100-Dieback)~PCA2+Phenology,olddata,family="quasibinomial")
anova(glm2q2,test='Chisq')
glm3q2 <- glm(cbind(Dieback,100-Dieback)~PCA2,olddata,family="quasibinomial")
anova(glm3q2,test='Chisq')
glm4q2 <- glm(cbind(Dieback,100-Dieback)~Phenology,olddata,family="quasibinomial")
anova(glm4q2,test='Chisq')

#################################
#make dieback plot
## build new data frame to house predictions
##it should cover the same range as the values you have recorded in oberevations
newdata <- data.frame(olddata$PCA1,olddata$Phenology); names(newdata) <- c("PCA1","Phenology")
newdata <- arrange(newdata,PCA1)
newdata
summary(newdata)

## make the predictions - notice the type='link' argument
preds <- predict(glm1q, newdata=newdata, se.fit=T)
#preds <- predict(glm1q, newdata=newdata,type= 'link', se.fit=T)

## add CI predictions to new data frame in the transformed range 
# as n = 50 for original data, using z=1.96 to get the 95% CIs is reasonable
newdata$fit_T <- preds$fit
newdata$upr_T <- preds$fit+preds$se.fit*1.96
newdata$lwr_T <- preds$fit-preds$se.fit*1.96

#back-transform them with plogis to get response range
newdata$fit <- plogis(newdata$fit_T)
newdata$upr <- plogis(newdata$upr_T) 
newdata$lwr <- plogis(newdata$lwr_T)

head(newdata) #note the extra columns in the dataframe
str(newdata)

ggplot(data=olddata, aes(x=PCA1,col=Phenology,shape=Phenology)) + 
  geom_point(aes(y = Dieback/100), size=4, show.legend = FALSE) + 
  geom_smooth(data=newdata,aes(x=PCA1, y=fit, ymin=lwr, ymax=upr),stat='identity', show.legend = FALSE)+  
  theme_classic(base_size=24) + 
  ylab("Branch dieback ratio") +   
  geom_ribbon(data=newdata,aes(ymin = lwr, ymax = upr, fill =Phenology), 
              alpha = 0.5, col = "transparent", show.legend = FALSE)+ 
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+ 
  scale_shape_manual(values = shapes)


## plot raw data and model prediction in the response range
#for this we will use ggplot
g2.2 <- ggplot(data=olddata, aes(x=PCA1,col=Phenology,shape=Phenology)) + 
  geom_point(aes(y = Dieback/100), size=3, show.legend = FALSE) + 
  geom_smooth(data=newdata,aes(x=PCA1, y=fit, ymin=lwr, ymax=upr),stat='identity', show.legend = FALSE, lwd = 1.2)+  
  theme_classic(base_size=15) + 
  ylab("Branch dieback ratio") +   
  geom_ribbon(data=newdata,aes(ymin = lwr, ymax = upr, fill =Phenology), 
              alpha = 0.3, col = "transparent", show.legend = FALSE)+ 
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+ 
  scale_shape_manual(values = shapes)
g2.2

####### combine the plots  ########

PCA_relationship <- ggarrange(g2.1,g2.2,g1.1,g1.2, ncol = 2, nrow = 2)
PCA_relationship

ggsave("figs/PCA_relationship3.tiff", height=10, width=10, dpi=600)
ggsave("figs/PCA_relationship3.pdf", height=10, width=10)
dev.off()
