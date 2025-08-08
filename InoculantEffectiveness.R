
### LOAD PACKAGES ###
library(multcomp)
library(tidyverse)
library(codyn)
library(vegan)
library(abdiv)
library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)
library(olsrr)
library(car)
library(MuMIn)


# Read in Nodule number and weight data #
Nodules <- read.csv("GreenhouseNoduleData_Su2024.csv") %>%
  mutate(NoduleNumber=ifelse(NoduleNumber=="DEAD", NA, NoduleNumber),
         TotalNoduleWeight=ifelse(TotalNoduleWeight=="DEAD"|TotalNoduleWeight=="NotDetected", NA, TotalNoduleWeight),
         # removing certain individuals according to the notes column #
         NoduleNumber = ifelse(Individual=="312"|Individual=="540"|Individual=="131"|Individual=="191"|
                                 Individual=="276"|Individual=="537"|Individual=="552"|Individual=="911"|
                                 Individual=="912"|Individual=="921"|Individual=="922"|
                                 Individual=="924"|Individual=="960"|Individual=="1000"|
                                 Individual=="1014"|Individual=="1041"|Individual=="1053"|
                                 Individual=="1126"|Individual=="1163"|Individual=="1164"|
                                 Individual=="1174"|Individual=="1006"|Individual=="1073"|
                                 Individual=="1074"|Individual=="316"|Individual=="457"|
                                 Individual=="458"|Individual=="515"|Individual=="544"|
                                 Individual=="902"|Individual=="903"|Individual=="904"|
                                 Individual=="925"|Individual=="926"|Individual=="928"|
                                 Individual=="976"|Individual=="995"|Individual=="996"|
                                 Individual=="1033"|Individual=="1034"|Individual=="1035"|
                                 Individual=="1045"|Individual=="1057"|Individual=="1058"|
                                 Individual=="1059"|Individual=="1060"|Individual=="1069"|
                                 Individual=="1070"|Individual=="1072"|Individual=="1081"|
                                 Individual=="1117"|Individual=="1118"|Individual=="1119"|
                                 Individual=="1120"|Individual=="1132", NA, NoduleNumber),
         TotalNoduleWeight = ifelse(Individual=="312"|Individual=="540"|Individual=="131"|Individual=="191"|
                                      Individual=="276"|Individual=="537"|Individual=="552"|Individual=="911"|
                                      Individual=="912"|Individual=="921"|Individual=="922"|
                                      Individual=="924"|Individual=="960"|Individual=="1000"|
                                      Individual=="1014"|Individual=="1041"|Individual=="1053"|
                                      Individual=="1126"|Individual=="1163"|Individual=="1164"|
                                      Individual=="1174"|Individual=="1006"|Individual=="1073"|
                                      Individual=="1074"|Individual=="316"|Individual=="457"|
                                      Individual=="458"|Individual=="515"|Individual=="544"|
                                      Individual=="1006"|Individual=="902"|Individual=="903"|
                                      Individual=="904"|Individual=="925"|Individual=="926"|
                                      Individual=="928"|Individual=="976"|Individual=="995"|
                                      Individual=="996"|Individual=="1033"|Individual=="1034"|
                                      Individual=="1035"|Individual=="1045"|Individual=="1057"|
                                      Individual=="1058"|Individual=="1059"|Individual=="1060"|
                                      Individual=="1069"|Individual=="1070"|Individual=="1072"|
                                      Individual=="1081"|Individual=="1117"|Individual=="1118"|
                                      Individual=="1119"|Individual=="1120"|Individual=="1132", NA, TotalNoduleWeight)) %>%
  filter(Treatment!="S2") %>%
  mutate(NoduleNumber = as.numeric(NoduleNumber)) %>%
  mutate(TotalNoduleWeight = as.numeric(TotalNoduleWeight))


# Read in ANPP data #
Biomass <- read.csv("Greenhouse_LegumeBiomass_Fall2024.csv")

data<- merge(Nodules, Biomass, by=c("Species", "Treatment", "Individual"), all=T) %>%
  mutate(TotalNoduleWeight = ifelse(NoduleNumber==0, NA, TotalNoduleWeight))


# visualizing outliers #
ggplot(data=data2, aes(x=Treatment, y=ANPP))+
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 3) +
  facet_wrap(~Species)
ggplot(data=data2, aes(x=Treatment, y=BNPP))+
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 3) +
  facet_wrap(~Species)
ggplot(data=data, aes(x=Treatment, y=NoduleNumber))+
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 3)
ggplot(data=data, aes(x=Treatment, y=TotalNoduleWeight))+
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 3)
ggplot(data=data, aes(x=Treatment, y=Soil_NO3))+
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 3) +
  facet_wrap(~Species)
ggplot(data=data, aes(x=Treatment, y=Soil_NH4))+
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 3) +
  facet_wrap(~Species)
ggplot(data=data, aes(x=Treatment, y=PercN))+
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 3) +
  facet_wrap(~Species)


# Read in soil N #
SoilN <- read.csv("Soil_N_InoculantEffectiveness.csv")


# filtering out outliers for each response variable within each treatment group #
data1 <- merge(data, SoilN, by=c("Individual"), all=T) %>%
  mutate(Soil_NO3 = ifelse(Soil_NO3=="<0.1", 0.01, Soil_NO3)) 

# Read in Leaf N #
LeafN <- read.csv("Leaf_N_InoculantEffectiveness.csv")

data2 <- merge (data1, LeafN, by=c("Individual"), all=T)%>%
  filter(Treatment!="NA") %>%
  group_by(Treatment, Species, Individual) %>%
  filter(Individual!=25 & Individual!=333 & Individual!=345 & Individual!=947 & Individual!=474 & Individual!=559 & 
           Individual!=124 & Individual!=256 & Individual!=483 & Individual!=1083 & Individual!=1084 & 
           Individual!=247 & Individual!=296 & Individual!=558 & Individual!=991 & Individual!=217 & 
           Individual!=261 & Individual!=593 & Individual!=350 & Individual!=481 & Individual!=106 &
           Individual!=164 & Individual!=181 & Individual!=156 & Individual!=125 & Individual!=1167)

data2$Soil_NO3 <- as.numeric(data2$Soil_NO3)


# Read in GCLog #
EthAreaHeight <- read.csv("Greenhouse_GClog.csv")

# Read in trt data #
Metadata <- read.csv("Greenhouse Master_forR.csv")

#merge data #
data2.5 <- merge(EthAreaHeight, Metadata, by=c("Individual"))
data3 <- merge(data2.5, data2, by=c("Individual", "Species", "Treatment"))

# get average Ethylene area per individual #
# standard 'curve' for 100ppm ethylene is "y = 41.278x - 1.7315", so if 5ml standard (100ppm) injected, curve area is 204.6585 #
AvgIndArea <- data3 %>%
  group_by(Species, Treatment, Individual) %>%
  summarise(AvgEthArea = mean(Ethylene_Area, na.rm=T), AvgEthPPM = (100*AvgEthArea)/204.6585) %>%
  filter(Treatment!="S2")

######################################################################################################
# Summary Statistics #
summary_table <- data2 %>%
  group_by(Treatment, Species) %>%
  summarise(mean_ANPP = mean(ANPP, na.rm = TRUE),mean_BNPP = mean(BNPP, na.rm = TRUE),
            
    mean_ANPP = mean(ANPP, na.rm = TRUE),
    median_ANPP = median(ANPP, na.rm = TRUE),
    n_ANPP = sum(!is.na(ANPP)),
    se_ANPP = sd(ANPP, na.rm = TRUE) / sqrt(n_ANPP),
    
    mean_BNPP = mean(BNPP, na.rm = TRUE),
    median_BNPP = median(BNPP, na.rm = TRUE),
    n_BNPP = sum(!is.na(BNPP)),
    se_BNPP = sd(BNPP, na.rm = TRUE) / sqrt(n_BNPP),
    
    mean_NoduleNumber = mean(NoduleNumber, na.rm = TRUE),
    median_NoduleNumber = median(NoduleNumber, na.rm = TRUE),
    n_NoduleNumber = sum(!is.na(NoduleNumber)),
    se_NoduleNumber = sd(NoduleNumber, na.rm = TRUE) / sqrt(n_NoduleNumber),
    
    mean_TotalNoduleWeight = mean(TotalNoduleWeight, na.rm = TRUE),
    median_TotalNoduleWeight = median(TotalNoduleWeight, na.rm = TRUE),
    n_TotalNoduleWeight = sum(!is.na(TotalNoduleWeight)),
    se_TotalNoduleWeight = sd(TotalNoduleWeight, na.rm = TRUE) / sqrt(n_TotalNoduleWeight),
    
    mean_Soil_NO3 = mean(Soil_NO3, na.rm = TRUE),
    median_Soil_NO3 = median(Soil_NO3, na.rm = TRUE),
    n_Soil_NO3 = sum(!is.na(Soil_NO3)),
    se_Soil_NO3 = sd(Soil_NO3, na.rm = TRUE) / sqrt(n_Soil_NO3),
    
    mean_Soil_NH4 = mean(Soil_NH4, na.rm = TRUE),
    median_Soil_NH4 = median(Soil_NH4, na.rm = TRUE),
    n_Soil_NH4 = sum(!is.na(Soil_NH4)),
    se_Soil_NH4 = sd(Soil_NH4, na.rm = TRUE) / sqrt(n_Soil_NH4),
    
    mean_PercN = mean(PercN, na.rm = TRUE),
    median_PercN = median(PercN, na.rm = TRUE),
    n_PercN = sum(!is.na(PercN)),
    se_PercN = sd(PercN, na.rm = TRUE) / sqrt(n_PercN)
  )

summary_table2 <- AvgIndArea %>%
  group_by(Treatment) %>%
  summarise(
    mean_PPM = mean(AvgEthPPM, na.rm = TRUE),
    median_PPM = median(AvgEthPPM, na.rm = TRUE),
    n_PPM = sum(!is.na(AvgEthPPM)),
    se_PPM = sd(AvgEthPPM, na.rm = TRUE) / sqrt(n_PPM))
##############################################################################################
# Define custom colors
my_colors <- c("Control" = "#704020", "S1" = "#8B8C64", "S3" = "#d17200")


# Analyses #

hist(data2$ANPP)
res_anpp <- aov(log(ANPP) ~ Treatment*Species, data = data2)
resanpp <- residuals(res_anpp, type="pearson")
plot(resanpp)
shapiro.test(residuals(res_anpp))
leveneTest(log(ANPP) ~ Treatment*Species, data = data2)

summary(res_anpp)
anpp_emm <- emmeans(res_anpp, ~ Treatment|Species, adjust="BH") 
pairs(anpp_emm)
###########################################################################


hist(data2$BNPP)
res_bnpp <- aov(log(BNPP) ~ Treatment*Species, data = data2)
resbnpp <- residuals(res_bnpp, type="pearson")
plot(rsbnpp)
shapiro.test(residuals(res_bnpp))
leveneTest(log(BNPP) ~ Treatment*Species, data = data2)

summary(res_bnpp)
bnpp_emm <- emmeans(res_bnpp, ~ Treatment|Species, adjust="BH") 
pairs(bnpp_emm)





########################################################################################################################
hist(data2$NoduleNumber)
res_nod <- aov(log1p(NoduleNumber) ~ Treatment*Species, data=data2)

resnod <- residuals(res_nod, type="pearson")
plot(resnod)
shapiro.test(residuals(res_nod))
leveneTest(log1p(NoduleNumber) ~ Treatment*Species, data = data2)

summary(res_nod)
nod_emm <- emmeans(res_nod, ~ Treatment|Species, adjust="BH") 
pairs(nod_emm)

#####################################################################################################################

hist(data2$TotalNoduleWeight)
res_weight <- aov(log(TotalNoduleWeight) ~ Treatment*Species, data=data2)

resweight <- residuals(res_weight, type="pearson")
plot(resweight)
shapiro.test(residuals(res_weight))
leveneTest(log(TotalNoduleWeight) ~ Treatment*Species, data = data2)

summary(res_weight)
weight_emm <- emmeans(res_weight, ~ Treatment|Species, adjust="BH") 
pairs(weight_emm)

#################################################################################################################
hist(data2$Soil_NO3)
res_NO3 <- aov(Soil_NO3 ~ Treatment*Species, data=data2)

resNO3 <- residuals(res_NO3, type="pearson")
plot(resNO3)
shapiro.test(residuals(res_NO3))
leveneTest(Soil_NO3 ~ Treatment*Species, data = data2)

summary(res_NO3)
NO3_emm <- emmeans(res_NO3, ~ Treatment|Species, adjust="BH") 
pairs(NO3_emm)

############################################################################################################

hist(data2$Soil_NH4)
res_NH4 <- aov(log(Soil_NH4) ~ Treatment*Species, data=data2)

resNH4 <- residuals(res_NH4, type="pearson")
plot(resNH4)
shapiro.test(residuals(res_NH4))
leveneTest(Soil_NH4 ~ Treatment*Species, data = data2)

summary(res_NH4)
NH4_emm <- emmeans(res_NH4, ~ Treatment|Species, adjust="BH")
pairs(NH4_emm)

#####################################################

hist(data2$PercN)
res_PercN <- aov(PercN ~ Treatment*Species, data=data2)

resPercN <- residuals(res_PercN, type="pearson")
plot(resPercN)
shapiro.test(residuals(res_PercN))
leveneTest(PercN ~ Treatment*Species, data = data2)


summary(res_PercN)
PercN_emm <- emmeans(res_PercN, ~ Treatment|Species, adjust="BH") 
pairs(PercN_emm)






############### BOXPLOTS ########################################################################

#save all as 1600x1600


ggplot(data = subset(data2, Species == "BA"), 
       aes(x = Treatment, y = ANPP, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1.0) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Shoot Biomass (g)") +
  xlab("Treatment Strain") +
  annotate("text", x = 1, y = 0.145, label = "a", size = 30) +
  annotate("text", x = 2, y = 0.330, label = "c", size = 30) +
  annotate("text", x = 3, y = 0.230, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "CN "), 
       aes(x = Treatment, y = ANPP, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Shoot Biomass (g)") +
  xlab("Treatment Strain") + ylim(0,0.3) +
  annotate("text", x = 1, y = 0.090, label = "a", size = 30) +
  annotate("text", x = 2, y = 0.265, label = "b", size = 30) +
  annotate("text", x = 3, y = 0.280, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "LH"), 
       aes(x = Treatment, y = ANPP, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Shoot Biomass (g)") +
  xlab("Treatment Strain") + 
  annotate("text", x = 1, y = 0.160, label = "a", size = 30) +
  annotate("text", x = 2, y = 0.345, label = "b", size = 30) +
  annotate("text", x = 3, y = 0.375, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))




ggplot(data = subset(data2, Species == "BA"), 
       aes(x = Treatment, y = BNPP, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1.0) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Root Biomass (g)") +
  xlab("Treatment Strain") + 
  annotate("text", x = 1, y = 1.210, label = "c", size = 30) +
  annotate("text", x = 2, y = 1.250, label = "b", size = 30) +
  annotate("text", x = 3, y = 0.550, label = "a", size = 30) +
  ylim(0,1.25) + 
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "CN "), 
       aes(x = Treatment, y = BNPP, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1.0) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Root Biomass (g)") +
  xlab("Treatment Strain") + 
  annotate("text", x = 1, y = 0.25, label = "a", size = 30) +
  annotate("text", x = 2, y = 0.35, label = "b", size = 30) +
  annotate("text", x = 3, y = 0.60, label = "c", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "LH"), 
       aes(x = Treatment, y = BNPP, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1.0) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Root Biomass (g)") +
  xlab("Treatment Strain") + 
  annotate("text", x = 1, y = 0.15, label = "a", size = 30) +
  annotate("text", x = 2, y = 0.25, label = "c", size = 30) +
  annotate("text", x = 3, y = 0.30, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))





ggplot(data = subset(data2, Species == "BA"), 
       aes(x = Treatment, y = NoduleNumber, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Nodule Number") +
  xlab("Treatment Strain") + 
  ylim(0,50) + 
  annotate("text", x = 1, y = 20, label = "a", size = 30) +
  annotate("text", x = 2, y = 48, label = "b", size = 30) +
  annotate("text", x = 3, y = 43, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "CN "), 
       aes(x = Treatment, y = NoduleNumber, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Nodule Number") +
  xlab("Treatment Strain") + ylim(0,100) +
  annotate("text", x = 1, y = 30, label = "a", size = 30) +
  annotate("text", x = 2, y = 100, label = "c", size = 30) +
  annotate("text", x = 3, y = 80, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "LH"), 
       aes(x = Treatment, y = NoduleNumber, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Nodule Number") +
  xlab("Treatment Strain") + 
  annotate("text", x = 1, y = 25, label = "a", size = 30) +
  annotate("text", x = 2, y = 48, label = "c", size = 30) +
  annotate("text", x = 3, y = 50, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))




ggplot(data = subset(data2, Species == "BA"), 
       aes(x = Treatment, y = TotalNoduleWeight, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Total Nodule Weight (g)") +
  xlab("Treatment Strain") +
  ylim(0, 0.08) + 
  annotate("text", x = 1, y = 0.045, label = "b", size = 30) +
  annotate("text", x = 2, y = 0.075, label = "b", size = 30) +
  annotate("text", x = 3, y = 0.035, label = "a", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "CN "), 
       aes(x = Treatment, y = TotalNoduleWeight, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Total Nodule Weight (g)") +
  xlab("Treatment Strain") +
  ylim(0,0.085) +
  annotate("text", x = 1, y = 0.02, label = "a", size = 30) +
  annotate("text", x = 2, y = 0.085, label = "b", size = 30) +
  annotate("text", x = 3, y = 0.085, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))
                            
ggplot(data = subset(data2, Species == "LH"), 
       aes(x = Treatment, y = TotalNoduleWeight, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Total Nodule Weight (g)") +
  xlab("Treatment Strain") +
  annotate("text", x = 1, y = 0.030, label = "a", size = 30) +
  annotate("text", x = 2, y = 0.060, label = "c", size = 30) +
  annotate("text", x = 3, y = 0.065, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))



ggplot(data = subset(data2, Species == "BA" & !is.na(Soil_NO3)), 
       aes(x = Treatment, y = Soil_NO3, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Soil Nitrate (ppm)") +
  xlab("Treatment Strain") +
  annotate("text", x = 1, y = 0.90, label = "a", size = 30) +
  annotate("text", x = 2, y = 1.00, label = "b", size = 30) +
  annotate("text", x = 3, y = 1.05, label = "b", size = 30) +
  ylim(0,1.05) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "CN " & !is.na(Soil_NO3)), 
       aes(x = Treatment, y = Soil_NO3, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Soil Nitrate (ppm)") +
  xlab("Treatment Strain") +
  annotate("text", x = 1, y = 0.55, label = "a", size = 30) +
  annotate("text", x = 2, y = 1.05, label = "b", size = 30) +
  annotate("text", x = 3, y = 0.95, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "LH" & !is.na(Soil_NO3)), 
       aes(x = Treatment, y = Soil_NO3, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Soil Nitrate (ppm)") +
  xlab("Treatment Strain") + 
  scale_y_continuous(breaks = c(0.0, 0.25, 0.50, 0.75, 1.00, 1.20)) + 
  annotate("text", x = 1, y = 0.95, label = "a", size = 30) +
  annotate("text", x = 2, y = 1.10, label = "b", size = 30) +
  annotate("text", x = 3, y = 1.25, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))



ggplot(data = subset(data2, Species == "BA" & !is.na(Soil_NH4)), 
       aes(x = Treatment, y = Soil_NH4, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Soil Ammonium (ppm)") +
  xlab("Treatment Strain") + ylim(0,5)+
  annotate("text", x = 1, y = 3.50, label = "a", size = 30) +
  annotate("text", x = 2, y = 4.40, label = "b", size = 30) +
  annotate("text", x = 3, y = 5.00, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "CN " & !is.na(Soil_NH4)), 
       aes(x = Treatment, y = Soil_NH4, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Soil Ammonium (ppm)") +
  xlab("Treatment Strain") +
  annotate("text", x = 1, y = 4.50, label = "a", size = 30) +
  annotate("text", x = 2, y = 4.00, label = "a", size = 30) +
  annotate("text", x = 3, y = 4.05, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "LH" & !is.na(Soil_NH4)), 
       aes(x = Treatment, y = Soil_NH4, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Soil Ammonium (ppm)") +
  xlab("Treatment Strain") +
  annotate("text", x = 1, y = 4.00, label = "a", size = 30) +
  annotate("text", x = 2, y = 3.60, label = "a", size = 30) +
  annotate("text", x = 3, y = 3.60, label = "a", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))



ggplot(data = subset(data2, Species == "BA"), 
       aes(x = Treatment, y = PercN, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Leaf Tissue Nitrogen (%)") +
  xlab("Treatment Strain") +
  annotate("text", x = 1, y = 3.15, label = "a", size = 30) +
  annotate("text", x = 2, y = 4.05, label = "b", size = 30) +
  annotate("text", x = 3, y = 3.95, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "CN "), 
       aes(x = Treatment, y = PercN, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Leaf Tissue Nitrogen (%)") +
  xlab("Treatment Strain") +
  annotate("text", x = 1, y = 2.90, label = "a", size = 30) +
  annotate("text", x = 2, y = 4.75, label = "c", size = 30) +
  annotate("text", x = 3, y = 4.15, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

ggplot(data = subset(data2, Species == "LH"), 
       aes(x = Treatment, y = PercN, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Leaf Tissue Nitrogen (%)") +
  xlab("Treatment Strain") + ylim(0,3) +
  annotate("text", x = 1, y = 2.75, label = "a", size = 30) +
  annotate("text", x = 2, y = 2.75, label = "b", size = 30) +
  annotate("text", x = 3, y = 2.80, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Commercial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))

#######################################################

Summary_Table2 <- AvgIndArea %>%
  group_by(Treatment, Species) %>%
  summarise(mean_Eth = mean(AvgEthArea, na.rm = TRUE),
          median_Eth = median(AvgEthArea, na.rm = TRUE),
          n_Eth = sum(!is.na(AvgEthArea)),
          se_Eth = sd(AvgEthArea, na.rm = TRUE) / sqrt(n_Eth))



hist(AvgIndArea$AvgEthPPM)
res_EthPPM <- aov(log1p(AvgEthPPM) ~ Treatment*Species, data=AvgIndArea)

resEthPPM <- residuals(res_EthPPM, type="pearson")
plot(resEthPPM)
shapiro.test(residuals(res_EthPPM))
leveneTest(log1p(AvgEthPPM) ~ Treatment*Species, data = AvgIndArea)


summary(res_EthPPM)
EthPPM_emm <- emmeans(res_EthPPM, ~ Treatment, adjust="BH") 
pairs(EthPPM_emm)


ggplot(data = AvgIndArea, 
       aes(x = Treatment, y = AvgEthPPM, color = Treatment)) +
  geom_boxplot(aes(group = Treatment), fill = NA, outlier.shape = NA, size = 6) +  # Boxplot outline only
  stat_summary(fun = mean, aes(group = Treatment), geom = "crossbar", width = 0.75, # Match the boxplot width
               color = "black", size = 1) +
  geom_jitter(width = 0.2, size = 10, alpha = 0.7) + # Raw points
  ylab("Ethylene Produced (ppm)") +
  xlab("Treatment Strain") +
  annotate("text", x = 1, y = 15.00, label = "a", size = 30) +
  annotate("text", x = 2, y = 160.00, label = "c", size = 30) +
  annotate("text", x = 3, y = 95.00, label = "b", size = 30) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(labels = c("Control", "Native", "Comemrcial")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 65),
        axis.text.x = element_text(size = 70),
        axis.text.y = element_text(size = 70),
        legend.position = "none",
        axis.ticks.length = unit(0.1, "inch"))



#########################################################################
# OLD GRAPHS #

ggplot(data=barGraphStats(data=subset(data2, Species=="BA"),variable="ANPP",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Shoot Biomass (g)") +
  xlab("Treatment") +
  annotate("text", x = 1, y = 0.075, label = "a", size=25) +
  annotate("text", x = 2, y = 0.125, label = "b", size=25) +
  annotate("text", x = 3, y = 0.095, label = "c", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(data2, Species=="CN "),variable="ANPP",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Shoot Biomass (g)") +
  xlab("Treatment") +
  scale_y_continuous(breaks = c(0.00, 0.03, 0.06, 0.09)) + 
  annotate("text", x = 1, y = 0.04, label = "a", size=25) +
  annotate("text", x = 2, y = 0.10, label = "b", size=25) +
  annotate("text", x = 3, y = 0.10, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(data2, Species=="LH"),variable="ANPP",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Shoot Biomass (g)") +
  xlab("Treatment") +
  scale_y_continuous(breaks = c(0.00, 0.03, 0.06, 0.09)) + 
  annotate("text", x = 1, y = 0.042, label = "a", size=25) +
  annotate("text", x = 2, y = 0.105, label = "b", size=25) +
  annotate("text", x = 3, y = 0.09, label = "c", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(data2, Species=="BA"),variable="NoduleNumber",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Nodule Number") +
  xlab("Treatment") +
  scale_y_continuous(breaks = c(0, 4, 8, 12)) + 
  annotate("text", x = 1, y = 2, label = "a", size=25) +
  annotate("text", x = 2, y = 13, label = "b", size=25) +
  annotate("text", x = 3, y = 11, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(data2, Species=="CN "),variable="NoduleNumber",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Nodule Number") +
  xlab("Treatment") +
  ylim(0, 51) +
  annotate("text", x = 1, y = 6, label = "a", size=25) +
  annotate("text", x = 2, y = 48, label = "b", size=25) +
  annotate("text", x = 3, y = 38, label = "c", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(data2, Species=="LH"),variable="NoduleNumber",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Nodule Number") +
  xlab("Treatment") +
  ylim(0,15) +
  annotate("text", x = 1, y = 3, label = "a", size=25) +
  annotate("text", x = 2, y = 14, label = "b", size=25) +
  annotate("text", x = 3, y = 12, label = "c", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(data2, Species=="BA"),variable="TotalNoduleWeight",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Total Nodule Weight (g)") +
  xlab("Treatment") +
  ylim(0.00,0.03) + 
  annotate("text", x = 1, y = 0.026, label = "a", size=25) +
  annotate("text", x = 2, y = 0.028, label = "a", size=25) +
  annotate("text", x = 3, y = 0.015, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(data2, Species=="CN "),variable="TotalNoduleWeight",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Total Nodule Weight (g)") +
  xlab("Treatment") +
  ylim(0.00,0.04) + 
  annotate("text", x = 1, y = 0.013, label = "a", size=25) +
  annotate("text", x = 2, y = 0.032, label = "b", size=25) +
  annotate("text", x = 3, y = 0.033, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(data2, Species=="LH"),variable="TotalNoduleWeight",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Total Nodule Weight (g)") +
  xlab("Treatment") +
  ylim(0.00,0.03) + 
  annotate("text", x = 1, y = 0.013, label = "a", size=25) +
  annotate("text", x = 2, y = 0.023, label = "b", size=25) +
  annotate("text", x = 3, y = 0.018, label = "a", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(BiomassSoilN, Species=="BA"),variable="Soil_NO3",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Soil Nitrate (ppm)") +
  xlab("Treatment") +
  ylim(0.0, 0.65) + 
  annotate("text", x = 1, y = 0.25, label = "a", size=25) +
  annotate("text", x = 2, y = 0.58, label = "b", size=25) +
  annotate("text", x = 3, y = 0.62, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(BiomassSoilN, Species=="CN "),variable="Soil_NO3",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Soil Nitrate (ppm)") +
  xlab("Treatment") +
  ylim(0.0, 0.65) + 
  annotate("text", x = 1, y = 0.21, label = "a", size=25) +
  annotate("text", x = 2, y = 0.60, label = "b", size=25) +
  annotate("text", x = 3, y = 0.605, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(BiomassSoilN, Species=="LH"),variable="Soil_NO3",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Soil Nitrate (ppm)") +
  xlab("Treatment") +
  ylim(0.0, 0.65) + 
  annotate("text", x = 1, y = 0.33, label = "a", size=25) +
  annotate("text", x = 2, y = 0.59, label = "b", size=25) +
  annotate("text", x = 3, y = 0.64, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=BiomassSoilN, variable="Soil_NH4",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Soil Ammonium (ppm)") +
  xlab("Treatment") +
  annotate("text", x = 1, y = 0.24, label = "a", size=25) +
  annotate("text", x = 2, y = 0.59, label = "a", size=25) +
  annotate("text", x = 3, y = 0.62, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(BiomassSoilLeafN, Species=="BA"),variable="PercN",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Leaf Nitrogen (%)") +
  xlab("Treatment") +
  ylim(0,3) +
  annotate("text", x = 1, y = 1.6, label = "a", size=25) +
  annotate("text", x = 2, y = 2.9, label = "b", size=25) +
  annotate("text", x = 3, y = 3.0, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(BiomassSoilLeafN, Species=="CN "),variable="PercN",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Leaf Nitrogen (%)") +
  xlab("Treatment") +
  annotate("text", x = 1, y = 2.3, label = "a", size=25) +
  annotate("text", x = 2, y = 3.8, label = "b", size=25) +
  annotate("text", x = 3, y = 3.5, label = "c", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=subset(BiomassSoilLeafN, Species=="LH"),variable="PercN",byFactorNames=c("Treatment")), aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Leaf Nitrogen (%)") +
  xlab("Treatment") +
  ylim(0,3) +
  annotate("text", x = 1, y = 1.6, label = "a", size=25) +
  annotate("text", x = 2, y = 2.2, label = "b", size=25) +
  annotate("text", x = 3, y = 2.2, label = "b", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 50),axis.text.x=element_text(size = 40), 
        legend.position="none",axis.text.y=element_text(size = 50),axis.ticks.length=unit(0.1,"inch"))
# 1600 x 1600 #

ggplot(data=barGraphStats(data=AvgIndArea, variable="AvgEthPPM",byFactorNames=c("Treatment")),
       aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_bar(stat='identity', position="dodge", width=0.75) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  ylab("Ethylene Produced (ppm)") +
  xlab("Treatment") +
  ylim(0,15) +
  annotate("text", x = 1, y = 2, label = "a", size=25) +
  annotate("text", x = 2, y = 15, label = "b", size=25) +
  annotate("text", x = 3, y = 9, label = "c", size=25) +
  scale_fill_manual(values = c("#704020", "#8B8C64", "#d17200")) +
  scale_x_discrete(labels=c("Control", "Typical Native Strain", "Restoration Strain")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 30), legend.position="none",
        axis.text.y=element_text(size = 40))
# 1400 x 1400 #