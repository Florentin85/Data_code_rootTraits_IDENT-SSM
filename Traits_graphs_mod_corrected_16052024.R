
#################################################################################
#### Title: 2nd chapter graphics and statistical models after laptop was ########
#### stolen and I could not re-create my results because of removed outliers ####
#################################################################################

########################################################
############# Date: 01.07.2023 #########################
############# Author: Florentin C. Jaeger ##############
########################################################

# Wipe the environment that I originally worked with in "Pilot_CV_jul13_2023.R"
# and continue with only the below workspace objects:
rm(list = setdiff(ls(), c("absorptive", "pilot", "prop_root_frac",
                          "collinearity_matrix", "new_names_traits",
                          "collinearity_matrix_z", "coll.mat.z.pca",
                          "pca_trait_var", "cor_matrix_traits",
                          "absorptive_cluster", "lm_probing_inter_rld",
                          "cat_plot_lmTraits")))

#### packages ####

install.packages("ggplot2")
install.packages("ggfortify")
install.packages("plyr") # Always load "plyr" before "dplyr"!!!
install.packages("tidyverse")
install.packages("viridis") # Color palette for handicapped people
install.packages("readxl")
install.packages("psych")
install.packages("vegan")
install.packages("factoextra")
install.packages("lme4")
install.packages("pbkrtest")
install.packages("lmerTest")
install.packages("effects")
install.packages("LMERConvenienceFunctions")
install.packages("sjPlot")
install.packages("emmeans")
install.packages("MuMIn")
# install.packages("pca3d")
# install.packages("magick")
install.packages("AICcmodavg")
install.packages("rJava") # "rJava" is needed for "xlsx"! sudo apt-get install default-jre / sudo apt-get install default-jdk / sudo R CMD javareconf
install.packages("xlsx")
install.packages("plotrix")
install.packages("car")
install.packages("cowplot")
install.packages("bestNormalize")
install.packages("ggcorrplot")
install.packages("Hmisc")
install.packages("interactions")
install.packages("visreg")
install.packages("car")
install.packages("bestNormalize")
install.packages("cowplot")

# library(raster) # coefficient of variation. Important to load before the tydiverse packages --> errors!!!
library(plyr) # Always load plyr before dplyr!
# detach("package:plyr", unload = TRUE)
library(ggplot2)
# library(ggfortify)
# library(viridis)
library(tidyverse)
library(magrittr) # To have pipe operator!
library(lme4)
library(pbkrtest)
library(lmerTest)
# library(effects)
# library(LMERConvenienceFunctions)
# library(sjPlot)
library(emmeans)
library(MuMIn)
# library(pca3d)
# library(magick)
# library(AICcmodavg)
library(car)
library(MASS)
library(bestNormalize)
# library(rJava)
# library(xlsx)
# library(plotrix)
library(cowplot)
# library(sjmisc)
# library(sjlabelled)
library(readxl)
library(psych) # supports ordination, factor analysis and structural equation modeling (suggests lavaan)!
library(vegan)
library(factoextra)
library(ggcorrplot)
library(Hmisc)
library(interactions)
library(visreg)
library(sjPlot)
library(effects)
# library(data.table)
# detach("package:data.table", unload = TRUE)

options(tibble.print_max = Inf) # set tibble options to print all rows of df

#### -------------------------------- remove RTD outlier from absorptive ----------------------------- ####

# Create copy to not meddle with the original (outlier included):
absorptive_rtd_outRem <- absorptive
# row = 151, col = 18

# Replace the rtd outlier value with NA
absorptive_rtd_outRem$rtd_g_cm3[151] <- NA

#### ------------------------------ Select variables for PCA ------------------------------- ####

collinearity_matrix <- subset(absorptive, select = c(frb_mg_cm3,
                                                     srl_m_g,
                                                     RLD_cm_cm3,
                                                     avgdiam_mm,
                                                     rootvol_aff_class_cm3,
                                                     rtd_g_cm3,
                                                     ntips_length_n_cm,
                                                     nforks_length_n_cm,
                                                     rdmc_dry_g_fresh_g,
                                                     sra_cm2_g))
# renaming needs as.data.frame:
collinearity_matrix <- as.data.frame(collinearity_matrix)

# rename variables for PCA:
# (It is possible, but very tedious to re-name variables within the PCA call!)
new_names_traits <- c("FRB", "SRL", "RLD", "Diam", "Vol",
                      "RTD", "RBD", "NF", "RDMC", "SRA")

# Directly assigning new names to colnames()
colnames(collinearity_matrix) <- new_names_traits

#### ----------------------------- Principal Component Analysis ----------------------------- ####

# different units need to be standardized for PCA (vegan package)
# Better to use method standardize because I want to give equal importance to all "species" and
# remove any scale-related bias in the analysis.
# hellinger on the other hand reduces the impact of rare species and emphasizes the importance
# of dominant species. It works well for data with many zeros (sparse data) and is often used
# in ecological distance metrics. --> 3rd (fungi) chapter (ASV tables)!
collinearity_matrix_z <- decostand(collinearity_matrix, method = "standardize")
apply(collinearity_matrix_z, 2, mean) # the data is now centered (means ~ 0)
apply(collinearity_matrix_z, 2, sd) # the data is now scaled (standard deviations = 1)

coll.mat.z.pca <- prcomp(collinearity_matrix_z) # Principal Components Analysis (PCA)
summary(coll.mat.z.pca)

### Simple PCA exploration plot with data points and traits:
biplot(coll.mat.z.pca)
# plot of eigenvalues ordered from largest to the smallest.
# The number of component is determined at the point, beyond which the remaining
# eigenvalues are all relatively small and of comparable size (Jollife 2002,
# Peres-Neto, Jackson, and Somers (2005)).
### Scree plot: ###
fviz_eig(coll.mat.z.pca, addlabels = TRUE, ylim = c(0, 50))
# The first axis do not decline very rapidly!

# Extract the results, for variables, from a PCA:
pca_trait_var <- get_pca_var(coll.mat.z.pca)
pca_trait_var$coord # coordinates of variables to create a scatter plot
pca_trait_var$cos2 # represents the quality of representation for variables on
# the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
pca_trait_var$contrib # contains the contributions (in percentage) of the variables
# to the principal components. The contribution of a variable (var) to a given principal
# component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).

### Total cos2 of variables on Dim.1 and Dim.2 (plot):
fviz_cos2(coll.mat.z.pca, choice = "var", axes = 1:2)
# FRB, RLD, NF, RBD, and Diam are most important.

### PCA graph of variables ###
# Color by cos2 values: quality on the factor map
# A high cos2 indicates a good representation of the variable on the principal component
# (variable positioned further outward).
# A low cos2 indicates that the variable is not perfectly represented by the PCs.
# In this case the variable is close to the center.
# For a given variable, the sum of the cos2 on all the principal components is equal to one
# If a variable is perfectly represented by only two principal components (Dim.1 & Dim.2),
# the sum of the cos2 on these two PCs is equal to one.
fviz_pca_var(coll.mat.z.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, title = "Principal Component Analysis - Root traits")


#### -------------------------------- pearson correlations -------------------------------- ####

# Calculate Pearson correlation matrix
cor_matrix_traits <- cor(collinearity_matrix, method = "pearson")
# ggplot2 needs as.data.frame:
cor_matrix_traits <- as.data.frame(cor_matrix_traits)

# Pearson's correlation coefficients in a heatmap:
ggcorrplot(cor_matrix_traits, hc.order = TRUE, type = "lower",
           lab = TRUE)

#### ------------------------------ Hirarchical variable clustering ------------------------------- ####

absorptive_cluster <- as.data.frame(absorptive)
names(absorptive_cluster)[50]<-"FRB"
names(absorptive_cluster)[13]<-"SRL"
names(absorptive_cluster)[16]<-"Diam"
names(absorptive_cluster)[17]<-"Vol"
names(absorptive_cluster)[18]<-"RTD"
names(absorptive_cluster)[22]<-"RBD"
names(absorptive_cluster)[24]<-"RDMC"
names(absorptive_cluster)[66]<-"SRA"
names(absorptive_cluster)[67]<-"RLD"

plot(varclus(~., data = absorptive_cluster[, -c(1,2,3,4,5,6,7,8,9,10,11,12,14,15,
                                                19,20,21,23,25,26,27,28,29,30,31,
                                                32,33,34,35,36,37,38,39,40,41,42,
                                                43,44,45,46,47,48,49,51,52,53,54,
                                                55,56,57,58,59,60,61,62,63,64,65)]),
     las = 1, cex.lab = 1.5)

# Reduced set of traits (variables) for 2nd paper manuscript:
plot(varclus(~., data = absorptive_cluster[, c(22, 66, 13, 16, 24, 18, 17, 50, 67)]),
     las = 1, cex.lab = 1.5)


#### ------------------------ 3-way-interaction-visualization --------------------------- ####

# https://interactions.jacob-long.com/reference/cat_plot.html#ref-usage --> from R package
# "interactions" from Jacob A. Long:
lm_probing_inter_rld <- lm(RLD_cm_cm3 ~ species * Soil_depth * h2o, data = absorptive)

cat_plot_lmTraits <- cat_plot(lm_probing_inter_rld, pred = species,
                              modx = Soil_depth, mod2 = h2o,
                              interval = TRUE, plot.points = TRUE,
                              legend.main = "Soil depth",
                              mod2.labels = c("H2O = High", "H2O = Low"),
                              main.title = "Three-way interactions based on linear regression")

# Customize the axis and legend labels using labs() from ggplot2
# label_value() should only displays the value of a factor. I don't know why it does not work!?
# I tried changing the facet_grid labels with expression, but that did also not work easily!
# Only do it if the reviewer asks for it!
cat_plot_lmTraits + labs(x = expression(),
                         y = expression(Root ~ length ~ density ~ "("*cm ~ cm^"-3"*")")) +
  theme(axis.text.x = element_text(face = "italic", angle = 25, hjust = 1))

#### ---------------------------------- Graphs "outlier" detection ---------------------------------- ####

### RBD ###
ggplot(data = absorptive, aes(y = ntips_length_n_cm, x = species, fill = h2o)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7),
               outlier.color = "white", notch = TRUE) + # Outliers from boxplot fnc. and geom_jitter()
  # show as double values if outlier color is not specified as white against the white background!
  geom_jitter() +
  labs(x = "tree species", y = "ntips_length_n_cm") +
  scale_fill_viridis_d(alpha = 1, begin = 0.4, end = 0.9, direction = 1, option = "B") + theme_light() +
  theme(aspect.ratio = 0.75, text = element_text(size = 16), axis.text = element_text(size = 14),
        axis.text.x = element_text(face = "italic", angle = 15, hjust = 1))

absorptive %>% group_by(ntips_length_n_cm) %>% summarise(n = n()) %>% filter(n>1)
# I have got quite a few non-unique values in my "ntips_length_n_cm" (later found out this was
# due to double plotting of boxplot outliers and geom_jitter() points)!
View(absorptive) # df not ordered according to h2o and species.
subset(absorptive$ntips_length_n_cm, absorptive$species == "Larix laricina")
# There is only one very high values = 7.99 --> so for some reason geom_jitter()
# displays this wrongly (Yes, double display as mentioned above)!

### RTD ###
ggplot(data = absorptive, aes(y = rtd_g_cm3, x = species, fill = h2o)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7),
               outlier.color = "white", notch = TRUE) +
  geom_jitter() +
  labs(x = "tree species", y = "rtd_g_cm3") +
  scale_fill_viridis_d(alpha = 1, begin = 0.4, end = 0.9, direction = 1, option = "B") + theme_light() +
  theme(aspect.ratio = 0.75, text = element_text(size = 16), axis.text = element_text(size = 14),
        axis.text.x = element_text(face = "italic", angle = 15, hjust = 1))

absorptive$rtd_g_cm3[151] # 18.617 g_cm-3
absorptive[151,] # 8_2_PIGL_0_5_High ==> This sample shows a very small volume and
# a very high weight (from 100 % fragments added) in the "RootsMonoSSM2018_forR_f100prop.xlsx"
# file. This sample also has mycelium present, and this could have increased the weight
# to volume ratio. The rootvol_aff_class_cm3 = 0.057 --> very low compared to the other
# samples! This might cause the density to be very high.
# --> Delete this outlier because it could have been an error coming from mycelium present in the sample!

### RTD after outlier (18.617) removed (NA) ###
ggplot(data = absorptive_rtd_outRem, aes(y = rtd_g_cm3, x = species, fill = h2o)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7),
               outlier.color = "white", notch = TRUE, na.rm = TRUE) +
  geom_jitter() + # "na.rm = TRUE" screws up jitter, better leave it out!
  labs(x = "tree species", y = "rtd_g_cm3") +
  scale_fill_viridis_d(alpha = 1, begin = 0.4, end = 0.9, direction = 1, option = "B") + theme_light() +
  theme(aspect.ratio = 0.75, text = element_text(size = 16), axis.text = element_text(size = 14),
        axis.text.x = element_text(face = "italic", angle = 15, hjust = 1))

absorptive_rtd_outRem$rtd_g_cm3>=2.5
absorptive_rtd_outRem$rtd_g_cm3[86] # 2.806 --> Check this value in excel table from the lab
# and check if mycelium was present. If not, I don't have further justification for removal!
absorptive_rtd_outRem[86,] # 6_13_PIST_05_Low ==> No mycelium present in excel table
# ("RootsMonoSSM2018_forR_f100prop.xlsx") --> No justification to remove!

### SRA ###
ggplot(data = absorptive, aes(y = sra_cm2_g, x = species, fill = h2o)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7),
               outlier.color = "white", notch = TRUE) +
  geom_jitter() +
  labs(x = "tree species", y = "sra_cm2_g") +
  scale_fill_viridis_d(alpha = 1, begin = 0.4, end = 0.9, direction = 1, option = "B") + theme_light() +
  theme(aspect.ratio = 0.75, text = element_text(size = 16), axis.text = element_text(size = 14),
        axis.text.x = element_text(face = "italic", angle = 15, hjust = 1))
# There's one extreme outlier around 2500 cm-2 g! --> Check it out!
# It's most likely Acer saccharum Low?
absorptive$sra_cm2_g >= 2400
absorptive$sra_cm2_g[146] # 2502.233 cm-2 g
absorptive[146,] # 3, 21, Acer, Low, 0-5
# I checked in the excel table "RootsMonoSSM2018_forR_f100prop.xlsx" and there was no mycelium present.
# Furthermore, the other variables seem normal, so there's no reason to remove this value.
# ==> No true outlier!

#### --------------------------------------- Summarize -------------------------------------------- ####

sum_absor <- absorptive %>%
  dplyr::select(species, h2o, Soil_depth, srl_m_g, length_aff_cm, avgdiam_mm,
                rootvol_aff_class_cm3, rtd_g_cm3, ntips_length_n_cm, rdmc_dry_g_fresh_g,
                frb_mg_cm3, FRB_g_m2, coarse_fine_ratio_gg, surfarea_cm2,
                sra_cm2_g, RLD_cm_cm3) %>%
  dplyr::group_by(species, h2o, Soil_depth) %>%
  na.omit() %>%
  dplyr::summarise(SRL_m_g = mean (srl_m_g),
                   SRL_m_g_se = (sd(srl_m_g) / sqrt(n())),
                   LENGTH_cm = mean (length_aff_cm),
                   LENGTH_cm_se = (sd(length_aff_cm) / sqrt(n())),
                   Diam_mm = mean (avgdiam_mm),
                   Diam_mm_se = (sd(avgdiam_mm) / sqrt(n())),
                   Vol_cm3 = mean (rootvol_aff_class_cm3),
                   Vol_cm3_se = (sd(rootvol_aff_class_cm3) / sqrt(n())),
                   RTD_g_cm3 = mean (rtd_g_cm3),
                   RTD_g_cm3_se = (sd(rtd_g_cm3) / sqrt(n())),
                   RBD_n_cm = mean (ntips_length_n_cm),
                   RBD_n_cm_se = (sd(ntips_length_n_cm) / sqrt(n())),
                   RDMC_g_g = mean (rdmc_dry_g_fresh_g),
                   RDMC_g_g_se = (sd(rdmc_dry_g_fresh_g) / sqrt(n())),
                   FRB_mg_m3 = mean (frb_mg_cm3),
                   FRB_mg_m3_se = (sd(frb_mg_cm3) / sqrt(n())),
                   FRB_G_M2 = mean (FRB_g_m2),
                   FRB_G_M2_se = (sd(FRB_g_m2) / sqrt(n())),
                   COARSE_FINE = mean (coarse_fine_ratio_gg),
                   COARSE_FINE_se = (sd(coarse_fine_ratio_gg) / sqrt(n())),
                   SURFAREA_cm2 = mean (surfarea_cm2),
                   SURFAREA_cm2_se = (sd(surfarea_cm2) / sqrt(n())),
                   SRA_cm2_g = mean (sra_cm2_g),
                   SRA_cm2_se = (sd(sra_cm2_g) / sqrt(n())),
                   RLD_CM_CM3 = mean (RLD_cm_cm3),
                   RLD_CM_CM3_se = (sd(RLD_cm_cm3) / sqrt(n())),) %>%
  dplyr::mutate(group = interaction(species, h2o, species = "_"))

#### ------------------------------- Summarize RTD with one outlier removed ------------------------------- ####

sum_absor_rtd_outRem <- absorptive_rtd_outRem %>%
  dplyr::select(species, h2o, Soil_depth, rtd_g_cm3) %>%
  dplyr::group_by(species, h2o, Soil_depth) %>%
  na.omit() %>%
  dplyr::summarise(RTD_g_cm3 = mean (rtd_g_cm3),
                   RTD_g_cm3_se = (sd(rtd_g_cm3) / sqrt(n())),) %>%
  dplyr::mutate(group = interaction(species, h2o, species = "_"))


### Check the results for standard error of the mean in dplyr code above:
SEM_min_NAs <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

SEM_min_NAs(absorptive$srl_m_g[which(absorptive$species == "Betula papyrifera" &
                                       absorptive$h2o == "High" &
                                       absorptive$Soil_depth == "0-5")], na.rm = TRUE)
# 6.892085 --> correct!
SEM_min_NAs(absorptive$srl_m_g[which(absorptive$species == "Betula papyrifera" &
                                       absorptive$h2o == "High" &
                                       absorptive$Soil_depth == "5-10")], na.rm = TRUE)
# 20.47002 --> correct!
# Check the result of the above SEM function:
sd(absorptive$srl_m_g[which(absorptive$species == "Betula papyrifera" &
                              absorptive$h2o == "High" &
                              absorptive$Soil_depth == "0-5")], na.rm = TRUE)
# 13.78417
(13.78417/sqrt(4)) # 6.892085 --> correct! 4 = number of degrees of freedom per
# species, soil depth, and water treatment grouping.

#### ---------------- summarize to species and H2O treatment level (12 means & se's) ------------------ ####

# (summarize RTD separately because of 1 outlier removed!)
sum_absor_spH2o <- absorptive %>%
  dplyr::select(species, h2o, ntips_length_n_cm, RLD_cm_cm3, srl_m_g,
                sra_cm2_g, rdmc_dry_g_fresh_g, avgdiam_mm, length_aff_cm,
                rootvol_aff_class_cm3) %>%
  dplyr::group_by(species, h2o) %>%
  na.omit() %>%
  dplyr::summarise(RBD_n_cm = mean (ntips_length_n_cm),
                   RBD_n_cm_se = (sd(ntips_length_n_cm) / sqrt(n())),
                   RLD_CM_CM3 = mean (RLD_cm_cm3),
                   RLD_CM_CM3_se = (sd(RLD_cm_cm3) / sqrt(n())),
                   SRL_m_g = mean (srl_m_g),
                   SRL_m_g_se = (sd(srl_m_g) / sqrt(n())),
                   SRA_cm2_g = mean (sra_cm2_g),
                   SRA_cm2_se = (sd(sra_cm2_g) / sqrt(n())),
                   RDMC_g_g = mean (rdmc_dry_g_fresh_g),
                   RDMC_g_g_se = (sd(rdmc_dry_g_fresh_g) / sqrt(n())),
                   Diam_mm = mean (avgdiam_mm),
                   Diam_mm_se = (sd(avgdiam_mm) / sqrt(n())),
                   LENGTH_cm = mean (length_aff_cm),
                   LENGTH_cm_se = (sd(length_aff_cm) / sqrt(n())),
                   Vol_cm3 = mean (rootvol_aff_class_cm3),
                   Vol_cm3_se = (sd(rootvol_aff_class_cm3) / sqrt(n())),) %>%
  dplyr::mutate(group = interaction(species, h2o, species = "_"))
# SRA is a lot lower than before! How can this be? I did not find an error in the code. Maybe, I took
# the mean of the mean before, or did not call summarise from dplyr?!

### summarize RTD ###
sum_absor_spH2o_rtd_outRem <- absorptive_rtd_outRem %>%
  dplyr::select(species, h2o, rtd_g_cm3) %>%
  dplyr::group_by(species, h2o) %>%
  na.omit() %>%
  dplyr::summarise(RTD_g_cm3 = mean (rtd_g_cm3),
                   RTD_g_cm3_se = (sd(rtd_g_cm3) / sqrt(n())),) %>%
  dplyr::mutate(group = interaction(species, h2o, species = "_"))


#### --------------------------------------- vertical barplots ---------------------------------------- ####

# (Source: ChatGPT and Google) If you have a continuous, numerical variable, it is typically more
# appropriate to use geom_col() rather than geom_bar(). Here's why: geom_col(): This geometry is specifically
# designed for creating bar plots with a continuous, numerical variable. The height of the bars directly
# represents the values of the variable.
# By default, geom_bar() calculates the count or frequency of observations within each category and
# represents them as bar heights. If your data is already summarised or includes values for
# y (height of the bars), use geom_col(). If, however, you want ggplot() to count up the number of
# rows in your dataset, use geom_bar(). geom_bar() basically runs the count() function and plots it.

# It's possible to add the data points to the barplot and eliminate the weakness
# of the error bars and bars to hide the individual data values (geom_jitter{ggplot2}).
# However, I would need to be on the "absorptive" data.frame level, and not the "sum_abs"
# data.frame level (geom_point only shows one point for the means)!
# geom_point(aes(y = RTD_g_cm3), color = "red", size = 3)
# In general it is better to use boxplots with jitter() and notch = TRUE

# Root branching density (n cm-1):
ggplot(data = sum_absor, aes(x = Soil_depth, y = RBD_n_cm, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = RBD_n_cm - RBD_n_cm_se,
                     ymax = RBD_n_cm + RBD_n_cm_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Root ~ branching ~ density ~ "("*n ~ cm^"-1"*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.03, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Root length density (cm cm-3):
# This variable stayed the same after laptop was stolen (old data, outliers included)
ggplot(data = sum_absor, aes(x = Soil_depth, y = RLD_CM_CM3, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = RLD_CM_CM3 - RLD_CM_CM3_se,
                     ymax = RLD_CM_CM3 + RLD_CM_CM3_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Root ~ length ~ density ~ "("*cm ~ cm^"-3"*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.03, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Specific root length (m g-1):
ggplot(data = sum_absor, aes(x = Soil_depth, y = SRL_m_g, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = SRL_m_g - SRL_m_g_se,
                     ymax = SRL_m_g + SRL_m_g_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Specific ~ root ~ length ~ "("*m ~ g^"-1"*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.03, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Root tissue density (g cm-3):
ggplot(data = sum_absor, aes(x = Soil_depth, y = RTD_g_cm3, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = RTD_g_cm3 - RTD_g_cm3_se,
                     ymax = RTD_g_cm3 + RTD_g_cm3_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Root ~ tissue ~ density ~ "("*g ~ cm^"-3"*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.03, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

### RTD outlier removed (18.617 outlier) ###
ggplot(data = sum_absor_rtd_outRem, aes(x = Soil_depth, y = RTD_g_cm3, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = RTD_g_cm3 - RTD_g_cm3_se,
                     ymax = RTD_g_cm3 + RTD_g_cm3_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Root ~ tissue ~ density ~ "("*g ~ cm^"-3"*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.53, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

### supplementary figures ### :
# Specific root area (cm-2 g):
ggplot(data = sum_absor, aes(x = Soil_depth, y = SRA_cm2_g, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = SRA_cm2_g - SRA_cm2_se,
                     ymax = SRA_cm2_g + SRA_cm2_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Specific ~ root ~ area ~ "("*cm^"-2" ~ g*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.03, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Root diameter (mm-1):
ggplot(data = sum_absor, aes(x = Soil_depth, y = Diam_mm, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = Diam_mm - Diam_mm_se,
                     ymax = Diam_mm + Diam_mm_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Root ~ diameter ~ "("*mm^"-1"*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.03, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))
# Root diameter stayed the same as before laptop was stolen!

# Root length (cm-1):
ggplot(data = sum_absor, aes(x = Soil_depth, y = LENGTH_cm, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = LENGTH_cm - LENGTH_cm_se,
                     ymax = LENGTH_cm + LENGTH_cm_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Root ~ length ~ "("*cm^"-1"*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.03, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))
# Root length stayed the same as from before laptop was stolen!

# Root dry matter content (g g-1):
ggplot(data = sum_absor, aes(x = Soil_depth, y = RDMC_g_g, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = RDMC_g_g - RDMC_g_g_se,
                     ymax = RDMC_g_g + RDMC_g_g_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Root ~ dry ~ matter ~ content ~ "("*g ~ g^"-1"*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.36, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))
# Stayed the same as before laptop was stolen!

# Root volume (cm-3):
ggplot(data = sum_absor, aes(x = Soil_depth, y = Vol_cm3, fill = h2o)) +
  geom_col(width = 0.7, position = position_dodge(0.7)) +
  facet_grid(.~ species) +
  geom_linerange(aes(ymin = Vol_cm3 - Vol_cm3_se,
                     ymax = Vol_cm3 + Vol_cm3_se),
                 position = position_dodge(width = 0.7)) +
  labs(x = expression(Soil ~ depth ~ "("*cm*")"),
       y = expression(Root ~ volume ~ "("*cm^"-3"*")")) +
  labs(fill = expression(H[2]*O)) +
  scale_fill_manual(values = c("#33CCFF", "#ff9900")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        legend.justification = c(0, 1), legend.position = c(0.03, 0.95),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))
# Stayed the same as before laptop was stolen!


#### ------------------------------ transformations ------------------------------ ####

### log transformation log(RBD_n_cm1 + 1) ###
hist(absorptive$ntips_length_n_cm, breaks = 20) # Positive skew!
# log(x+1) for positively skewed data.
ntips_length_n_cm_log1 <- log(absorptive$ntips_length_n_cm + 1)
hist(absorptive$ntips_length_n_cm_log1, breaks = 20) # Much better!
absorptive <- cbind(absorptive, ntips_length_n_cm_log1)
names(absorptive)[68] <- "ntips_length_n_cm_log1"

### log(x + 1) on RLD_cm_cm3 ###
hist(absorptive$RLD_cm_cm3, breaks = 20) # Positive skew!
# log(x+1) for positively skewed data.
RLD_cm_cm3_log1 <- log(absorptive$RLD_cm_cm3 + 1)
hist(RLD_cm_cm3_log1, breaks = 20) # Much better!
absorptive <- cbind(absorptive, RLD_cm_cm3_log1)
names(absorptive)[69] <- "RLD_cm_cm3_log1"

### log(x + 1) on SRL_m_g ###
hist(absorptive$srl_m_g, breaks = 20) # Positive skew, but not very strongly!
# Two extreme high outliers.
SRL_m_g_log1 <- log(absorptive$srl_m_g + 1)
hist(SRL_m_g_log1, breaks = 20) # Better but not very good!
absorptive <- cbind(absorptive, SRL_m_g_log1)
names(absorptive)[70] <- "SRL_m_g_log1"

### transformation log(RTD + 1) ###
hist(absorptive_rtd_outRem$rtd_g_cm3, breaks = 20) # Strong positive skew, and
# 4 (untrue) high outliers!
RTD_g_cm3_log1 <- log(absorptive_rtd_outRem$rtd_g_cm3 + 1)
hist(RTD_g_cm3_log1, breaks = 20) # It did not really resolve the right skew!
absorptive_rtd_outRem <- cbind(absorptive_rtd_outRem, RTD_g_cm3_log1)
names(absorptive_rtd_outRem)[71] <- "RTD_g_cm3_log1"

### RTD log10(x+1) ###
RTD_g_cm3_log10 <- log10(absorptive_rtd_outRem$rtd_g_cm3 + 1)
hist(RTD_g_cm3_log10, breaks = 20) # It did not really resolve the right skew!
absorptive_rtd_outRem <- cbind(absorptive_rtd_outRem, RTD_g_cm3_log10)
names(absorptive_rtd_outRem)[72] <- "RTD_g_cm3_log10"

### RTD sqrt(x) ###
RTD_g_cm3_sqrt <- sqrt(absorptive_rtd_outRem$rtd_g_cm3)
hist(RTD_g_cm3_sqrt, breaks = 20) # Yes, it resolved the right skewness and
# made rtd a little more normal!
absorptive_rtd_outRem <- cbind(absorptive_rtd_outRem, RTD_g_cm3_sqrt)
names(absorptive_rtd_outRem)[73] <- "RTD_g_cm3_sqrt"

### RTD boxcox transformation: ###
# boxcox can't handle zero, negative, or NA values, add constant:
rtd_g_cm3_bc001 <- (absorptive_rtd_outRem$rtd_g_cm3 + 0.001)
absorptive_rtd_outRem <- cbind(absorptive_rtd_outRem, rtd_g_cm3_bc001)
names(absorptive_rtd_outRem)[74] <- "rtd_g_cm3_bc001"
shapiro.test(absorptive_rtd_outRem$rtd_g_cm3_bc001) # p-value = < 2.2e-16
# It's better to use visuals rather than shapiro.test!

# Boxcox is an intense transformation. Only use it when log, sqrt or Poisson
# distribution do not fix the problem!
bestNormalize::boxcox(absorptive_rtd_outRem$rtd_g_cm3_bc001)
# Estimated statistics:
# lambda = 0.04698855
# mean (before standardization) = -1.993871
# sd (before standardization) = 0.8896522

### Supplementary ###
### log(x + 1) on sra_cm2_g ###
hist(absorptive$sra_cm2_g, breaks = 20) # One extreme positive outlier!
# log(x+1) for positively skewed data.
SRA_cm2_g_log1 <- log(absorptive$sra_cm2_g + 1)
hist(SRA_cm2_g_log1, breaks = 20) # Much better!
absorptive <- cbind(absorptive, SRA_cm2_g_log1)
names(absorptive)[71] <- "SRA_cm2_g_log1"

### sqrt(x + 1) on sra_cm2_g ###
hist(absorptive$sra_cm2_g, breaks = 20) # One extreme positive outlier!
# sqrt(x+1) for positively skewed data.
SRA_cm2_g_sqrt1 <- sqrt(absorptive$sra_cm2_g + 1)
hist(SRA_cm2_g_sqrt1, breaks = 20) # Does not deal as well with the positive outlier as log()!
absorptive <- cbind(absorptive, SRA_cm2_g_sqrt1)
names(absorptive)[72] <- "SRA_cm2_g_sqrt1"

### log(x + 1) on avgdiam_mm ###
hist(absorptive$avgdiam_mm, breaks = 20) # Looks normal!
avgdiam_mm_log1 <- log(absorptive$avgdiam_mm + 1)
hist(avgdiam_mm_log1, breaks = 20) # Looks almost the same as the original variable!
absorptive <- cbind(absorptive, avgdiam_mm_log1)
names(absorptive)[74] <- "avgdiam_mm_log1"

### log(x+1) on length_aff_cm ###
hist(absorptive$length_aff_cm, breaks = 20) # Positive skew!
length_aff_cm_log1 <- log(absorptive$length_aff_cm + 1)
hist(length_aff_cm_log1, breaks = 20) # sqrt transformation seems to be better in this case!
absorptive <- cbind(absorptive, length_aff_cm_log1)
names(absorptive)[74] <- "length_aff_cm_log1"

### sqrt(x+1) on length_aff_cm ###
hist(absorptive$length_aff_cm, breaks = 20) # Positive skew!
length_aff_cm_sqrt1 <- sqrt(absorptive$length_aff_cm + 1)
hist(length_aff_cm_sqrt1, breaks = 20) # sqrt transformation seems to be better in this case!
absorptive <- cbind(absorptive, length_aff_cm_sqrt1)
names(absorptive)[76] <- "length_aff_cm_sqrt1"

### log(x+1) on rdmc_dry_g_fresh_g ###
hist(absorptive$rdmc_dry_g_fresh_g, breaks = 20) # Looks actually quite normal, slight positive skew!
rdmc_dry_g_fresh_g_log1 <- log(absorptive$rdmc_dry_g_fresh_g + 1)
hist(rdmc_dry_g_fresh_g_log1, breaks = 20) # Logarithmic transformation looks good!
absorptive <- cbind(absorptive, rdmc_dry_g_fresh_g_log1)
names(absorptive)[75] <- "rdmc_dry_g_fresh_g_log1"

### sqrt(x+1) on rdmc_dry_g_fresh_g ###
hist(absorptive$rdmc_dry_g_fresh_g, breaks = 20) # Looks actually quite normal, slight positive skew!
rdmc_dry_g_fresh_g_sqrt1 <- sqrt(absorptive$rdmc_dry_g_fresh_g + 1)
hist(rdmc_dry_g_fresh_g_sqrt1, breaks = 20) # Both log and sqrt look good!
absorptive <- cbind(absorptive, rdmc_dry_g_fresh_g_sqrt1)
names(absorptive)[77] <- "rdmc_dry_g_fresh_g_sqrt1"

### log(x+1) on rootvol_aff_class_cm3 ###
hist(absorptive$rootvol_aff_class_cm3, breaks = 20) # Positive skew!
rootvol_aff_class_cm3_log1 <- log(absorptive$rootvol_aff_class_cm3 + 1)
hist(rootvol_aff_class_cm3_log1, breaks = 20) # Looks like normalization was successful!
absorptive <- cbind(absorptive, rootvol_aff_class_cm3_log1)
names(absorptive)[78] <- "rootvol_aff_class_cm3_log1"

### sqrt(x+1) on rootvol_aff_class_cm3 ###
hist(absorptive$rootvol_aff_class_cm3, breaks = 20) # Positive skew!
rootvol_aff_class_cm3_sqrt1 <- sqrt(absorptive$rootvol_aff_class_cm3 + 1)
hist(rootvol_aff_class_cm3_sqrt1, breaks = 20) # Looks like normalization was successful!
absorptive <- cbind(absorptive, rootvol_aff_class_cm3_sqrt1)
names(absorptive)[79] <- "rootvol_aff_class_cm3_sqrt1"


#### --------------- Restricted maximum likelihood-based mixed models (MMs) ------------------ ####

### RBD_n_cm1 ###
lmer_rbd_spDh2o <- lme4::lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                              + (1|bloc_pairs)
                              + (1|bloc_pairs : species)
                              + (1|bloc_pairs : Soil_depth)
                              + (1|bloc_pairs : h2o)
                              + (1|block_plot),
                              data = absorptive, REML = TRUE)
# Use REML = TRUE if I have random effects present in my data!

summary(lmer_rbd_spDh2o)
car::Anova(lmer_rbd_spDh2o, type = 3)
library(visreg)
# Warning messages:
# 1: In printHypothesis(L, rhs, names(b)) :
#  one or more coefficients in the hypothesis include
# arithmetic operators in their names;
# the printed representation of the hypothesis will be omitted
# [02/08/2023] ChatGPT: In R, when you specify a hypothesis using the car::Anova() function,
# it checks the coefficients' names in your model to define the hypothesis. If any coefficient
# names include arithmetic operators (e.g., *, :, +, -, etc.), it can lead to issues with the printed
# representation of the hypothesis, and thus, the hypothesis will be omitted from the output.
# To resolve this warning message, you need to ensure that there are no arithmetic operators in
# the names of the coefficients in your hypothesis. One common reason for having arithmetic operators
# in the names is when you create interactions between variables in your model formula,
# as you have done in your lmer() function.
# I think I can safely ignore this warning message!!!
anova(lmer_rbd_spDh2o)
plot(lmer_rbd_spDh2o) # Plot looks funnel shaped. Variance increases with the mean.
# Try a log transformation!

### log(x + 1) transformed RBD_n_cm1 ###
lmer_rbd_spDh2o_log1 <- lme4::lmer(ntips_length_n_cm_log1 ~ species * Soil_depth * h2o
                                   + (1|bloc_pairs)
                                   + (1|bloc_pairs : species)
                                   + (1|bloc_pairs : Soil_depth)
                                   + (1|bloc_pairs : h2o)
                                   + (1|block_plot),
                                   data = absorptive, REML = TRUE)

plot(lmer_rbd_spDh2o_log1) # Looks much better!
summary(as_lmerModLmerTest(lmer_rbd_spDh2o_log1)) # Yes, produces results!
anova(as_lmerModLmerTest(lmer_rbd_spDh2o_log1), type = 3)
car::Anova(lmer_rbd_spDh2o_log1, type = 3)
qqnorm(resid(lmer_rbd_spDh2o_log1))
qqline(resid(lmer_rbd_spDh2o_log1)) # That residuals leave fitted line at both ends is normal!
hist(resid(lmer_rbd_spDh2o_log1), breaks = 100)
car::vif(lmer_rbd_spDh2o_log1)
# As a rule of thumb, a vif score over 5 is a problem.
# A score over 10 should be remedied (and you should consider dropping the problematic
# variable from the regression model or creating an index of all the closely related variables).
# Okay, the variance inflation is huge here (in the 1000s)! I am not sure if I should believe this?
# According to Stephane the high VIF is not a problem, because I work with exploratory models and
# not predictive models! I do not want to predict values of a different sample population
# with my estimates. Instead, I just want to know if my predictors have a significant effect
# on my response.
r.squaredGLMM(lmer_rbd_spDh2o_log1) # R2m = 0.4791647, R2c = 0.4847937 (Good fit!)

fixef(lmer_rbd_spDh2o_log1)
ranef(lmer_rbd_spDh2o_log1)
# [24.07.2023] Stephane Daigle said again that the zeros in my bloc : species interactions
# are not a problem, and that I should keep working with the full model!
coef(lmer_rbd_spDh2o_log1)
VarCorr(lmer_rbd_spDh2o_log1)
plot(lmer_rbd_spDh2o_log1, form = resid(.) ~ fitted(.) | species)
plot(lmer_rbd_spDh2o_log1, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_rbd_spDh2o_log1, form = resid(.) ~ fitted(.) | h2o)
plot(lmer_rbd_spDh2o_log1, form = ntips_length_n_cm_log1 ~ resid(.))
plot(ranef(lmer_rbd_spDh2o_log1))
lattice::qqmath(ranef(lmer_rbd_spDh2o_log1))

### Model visualization ###
# Two-way interaction (species * h2o):
plot(allEffects(lmer_rbd_spDh2o_log1), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_rbd_spDh2o_log1))
# Summary statistics:
plot_model(lmer_rbd_spDh2o_log1, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_rbd_spDh2o_log1, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Betula, Quercus, Pinus, and (Picea) show the strongest response!
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_rbd_spDh2o_log1, type = "eff", terms = c("species"),
           show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_rbd_spDh2o_log1, type = "eff", terms = c("h2o"),
           show.data = T) + theme_bw()
# Overall effect of water treatment per species:
emmip(lmer_rbd_spDh2o_log1, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)
# Betula and Quercus show strongest response!
# Plot by h2o over species (Not too useful!):
visreg(lmer_rbd_spDh2o_log1, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Both h2o = high and low follow the same pattern!
# Plot by species over h2o (Better visualization!):
visreg(lmer_rbd_spDh2o_log1, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)

### emmeans (Tukey’s Honest Significant Difference Test) ###

# Carsten F. Dormann: Parametrische Statistik: Verteilungen, maximum likelihood
# und GLM in R
# S. 201: Wenn wir aber partout eine post-hoc-Test durchführen wollen, dann müssen wir
# zumindest für die Anzahl Vergleiche korrigieren, um so tatsächlich bei 5 % Irrtums-
# wahrscheinlichkeit zu bleiben (Day & Quinn 1989, stellen ein hervorragende Über-
# sicht zur Verfügung). Dazu sind zwei Ansätze verbreitet: die Bonferroni-Korrektur
# und Tukey’s honest significant difference Test.
# Tukey’s Honest Significant Difference Test (S. 202):
# Tukey’s HSD ist noch strikter in der Korrektur der paarweisen Vergleiche. Er ist
# sozusagen der akzeptiere Maximalstandard. Wenn post-hoc-Vergleiche im Tukey’s
# HSD immer noch signifikant sind, dann gibt es daran nichts zu rütteln.
# Technisch gesehen ist Tukey’s HSD eine Variante des paarweisen t-Tests, nur
# wird statt der t-Verteilung eine studentised range distribution benutzt, weshalb der
# Test auch als Tukey’s range test bezeichnet wird.
# https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html
# pairwise comparisons default to adjust = "tukey", i.e., the Tukey HSD method.
# An adjustment method that is usually appropriate is Bonferroni; however, it can be quite conservative.
# --> see Nakagawa (2004): A farewell to Bonferroni: the problems of low statistical power and publication bias
# https://stackoverflow.com/questions/65002727/is-using-adjust-tukey-in-emmeans-equivalent-to-a-tukey-hsd-test
# Tukey adjustment: That adjustment is only appropriate for a single set of pairwise comparisons.
# If you specify adjust = "tukey" for non-pairwise comparisons or arbitrary contrasts,
# it will overrule you and use the "sidak" adjustment instead.
# I do pairwise comparisons on h2o High vs. Low, so emmeans should use the Tukey HSD!
# [28.07.2023] ChatGPT: Keep in mind that the decision to conduct post-hoc tests
# should be guided by your research questions and hypotheses. If you have specific
# research questions about the comparisons between certain groups, you may perform
# appropriate post-hoc tests even if the three-way interaction is not significant.
# However, it is essential to be cautious about multiple testing and consider
# controlling for Type I error rates (e.g., using Bonferroni correction) to avoid
# false-positive results.
# The above sources and reasoning provides the argument to use in my papers for
# conducting post-hoc tests in the absence of significant main effects (type 3 ANOVA)!

emm_lmer_rbd_spDh2o_log1 <- emmeans::emmeans(lmer_rbd_spDh2o_log1,
                                             specs = list(pairwise ~ h2o | species | Soil_depth),
                                             type = "response", adjust = "tukey")
# Results make sense given the figure.
# Below plot is overplotted!
plot(emm_lmer_rbd_spDh2o_log1, comparisons = TRUE)

emmip(lmer_rbd_spDh2o_log1, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = ntips_length_n_cm_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ branching ~ density)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))

### RLD_cm_cm3 ###
lmer_rld_spDh2o <- lme4::lmer(RLD_cm_cm3 ~ species * Soil_depth * h2o
                              + (1|bloc_pairs)
                              + (1|bloc_pairs : species)
                              + (1|bloc_pairs : Soil_depth)
                              + (1|bloc_pairs : h2o)
                              + (1|block_plot),
                              data = absorptive, REML = TRUE)

summary(lmer_rld_spDh2o)
car::Anova(lmer_rld_spDh2o, type = 3)
plot(lmer_rld_spDh2o) # Funnel shaped! This comes from the positively skewed (right skewed)
# response variable RLD_cm_cm3.

### log(x + 1) transformed RLD_cm_cm3 ###
lmer_rld_spDh2o_log1 <- lme4::lmer(RLD_cm_cm3_log1 ~ species * Soil_depth * h2o
                                   + (1|bloc_pairs)
                                   + (1|bloc_pairs : species)
                                   + (1|bloc_pairs : Soil_depth)
                                   + (1|bloc_pairs : h2o)
                                   + (1|block_plot),
                                   data = absorptive, REML = TRUE)

plot(lmer_rld_spDh2o_log1) # Much better than without log(x+1) transformation!
summary(as_lmerModLmerTest(lmer_rld_spDh2o_log1)) # Yes, produces results!
anova(as_lmerModLmerTest(lmer_rld_spDh2o_log1), type = 3)
car::Anova(lmer_rld_spDh2o_log1, type = 3)
qqnorm(resid(lmer_rld_spDh2o_log1))
qqline(resid(lmer_rld_spDh2o_log1))
hist(resid(lmer_rld_spDh2o_log1), breaks = 100)
car::vif(lmer_rld_spDh2o_log1)
r.squaredGLMM(lmer_rld_spDh2o_log1) # R2m = 0.691563, R2c = 0.7563457

fixef(lmer_rld_spDh2o_log1)
ranef(lmer_rld_spDh2o_log1)
coef(lmer_rld_spDh2o_log1)
VarCorr(lmer_rld_spDh2o_log1)
plot(lmer_rld_spDh2o_log1, form = resid(.) ~ fitted(.) | species)
plot(lmer_rld_spDh2o_log1, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_rld_spDh2o_log1, form = resid(.) ~ fitted(.) | h2o)
plot(lmer_rld_spDh2o_log1, form = RLD_cm_cm3_log1 ~ resid(.))
plot(ranef(lmer_rld_spDh2o_log1))
lattice::qqmath(ranef(lmer_rld_spDh2o_log1))

### Model visualization ###
# Two-way interaction (species * h2o):
plot(allEffects(lmer_rld_spDh2o_log1), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Only Larix looks significant in 0-5 cm!
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_rld_spDh2o_log1))
# Summary statistics:
plot_model(lmer_rld_spDh2o_log1, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_rld_spDh2o_log1, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# In the overall plots, Pinus looks actually more significant!
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_rld_spDh2o_log1, type = "eff", terms = c("species"),
           show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_rld_spDh2o_log1, type = "eff", terms = c("h2o"),
           show.data = T) + theme_bw()
# Overall effect of water treatment per species:
emmip(lmer_rld_spDh2o_log1, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)

# Plot by h2o over species (Not too useful!):
visreg(lmer_rld_spDh2o_log1, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Both h2o = high and low follow the same pattern!
# Plot by species over h2o (Better visualization!):
visreg(lmer_rld_spDh2o_log1, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_rld_spDh2o_log1 <- emmeans::emmeans(lmer_rld_spDh2o_log1,
                                             specs = list(pairwise ~ h2o | species | Soil_depth),
                                             type = "response", adjust = "tukey")
# Results make sense given the figure, and stayed the same from before laptop was stolen!

emmip(lmer_rld_spDh2o_log1, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = RLD_cm_cm3_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ length ~ density)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))

### SRL_m_g1 ###
lmer_srl_spDh2o <- lme4::lmer(srl_m_g ~ species * Soil_depth * h2o
                              + (1|bloc_pairs)
                              + (1|bloc_pairs : species)
                              + (1|bloc_pairs : Soil_depth)
                              + (1|bloc_pairs : h2o)
                              + (1|block_plot),
                              data = absorptive, REML = TRUE)

summary(lmer_srl_spDh2o)
car::Anova(lmer_srl_spDh2o, type = 3)
plot(lmer_srl_spDh2o) # Yikes! Narrow variance band, and then increases
# exponentially at high mean values.

### log(x + 1) transformed SRL_m_g1 ###
lmer_srl_spDh2o_log1 <- lme4::lmer(SRL_m_g_log1 ~ species * Soil_depth * h2o
                                   + (1|bloc_pairs)
                                   + (1|bloc_pairs : species)
                                   + (1|bloc_pairs : Soil_depth)
                                   + (1|bloc_pairs : h2o)
                                   + (1|block_plot),
                                   data = absorptive, REML = TRUE)

plot(lmer_srl_spDh2o_log1)
summary(as_lmerModLmerTest(lmer_srl_spDh2o_log1)) # Yes, produces results!
anova(as_lmerModLmerTest(lmer_srl_spDh2o_log1), type = 3)
car::Anova(lmer_srl_spDh2o_log1, type = 3)
# Stephane Daigle [03.08.2023]: ... just a few outliers is not a problem if the rest of
# the distribution is homogeneous, which is the case for this model after transformation
# (lmer_srl_spDh2o_log1). So yes, you can proceed with the transformed response variable.
qqnorm(resid(lmer_srl_spDh2o_log1))
qqline(resid(lmer_srl_spDh2o_log1))
hist(resid(lmer_srl_spDh2o_log1), breaks = 100)
car::vif(lmer_srl_spDh2o_log1)
r.squaredGLMM(lmer_srl_spDh2o_log1) # R2m = 0.5543036, R2c = 0.5798121

fixef(lmer_srl_spDh2o_log1)
ranef(lmer_srl_spDh2o_log1)
coef(lmer_srl_spDh2o_log1)
VarCorr(lmer_srl_spDh2o_log1)
plot(lmer_srl_spDh2o_log1, form = resid(.) ~ fitted(.) | species)
plot(lmer_srl_spDh2o_log1, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_srl_spDh2o_log1, form = resid(.) ~ fitted(.) | h2o)
plot(lmer_srl_spDh2o_log1, form = SRL_m_g_log1 ~ resid(.)) # Two outliers!
plot(ranef(lmer_srl_spDh2o_log1))
lattice::qqmath(ranef(lmer_srl_spDh2o_log1))

### Model visualization ###
# Two-way interaction (species * h2o):
plot(allEffects(lmer_srl_spDh2o_log1), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Larix looks significant in 0-5 cm and 20-30 cm depth, but shows opposite pattern!
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_srl_spDh2o_log1))
# Summary statistics:
plot_model(lmer_srl_spDh2o_log1, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_srl_spDh2o_log1, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_srl_spDh2o_log1, type = "eff", terms = c("species"),
           show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_srl_spDh2o_log1, type = "eff", terms = c("h2o"),
           show.data = T) + theme_bw()
# Overall effect of water treatment per species:
emmip(lmer_srl_spDh2o_log1, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)

# Plot by h2o over species (Not too useful!):
visreg(lmer_srl_spDh2o_log1, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Both h2o = high and low follow the same pattern!
# Plot by species over h2o (Better visualization!):
visreg(lmer_srl_spDh2o_log1, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_srl_spDh2o_log1 <- emmeans::emmeans(lmer_srl_spDh2o_log1,
                                             specs = list(pairwise ~ h2o | species | Soil_depth),
                                             type = "response", adjust = "tukey")

emmip(lmer_srl_spDh2o_log1, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = SRL_m_g_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Specific ~ root ~ length)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))
# Indeed, Larix in 0-5 and 20-30 shows the strongest response!

### RTD_g_cm3 outlier (18.617) removed ###
lmer_rtd_spDh2o <- lme4::lmer(rtd_g_cm3 ~ species * Soil_depth * h2o
                              + (1|bloc_pairs)
                              + (1|bloc_pairs : species)
                              + (1|bloc_pairs : Soil_depth)
                              + (1|bloc_pairs : h2o)
                              + (1|block_plot),
                              data = absorptive_rtd_outRem, REML = TRUE)
# I introduced one NA for the extreme outlier value. na.action --> The default action
# (na.omit, inherited from the 'factory fresh' value of getOption("na.action"))
# strips any observations with any missing values in any variables.
summary(lmer_rtd_spDh2o)
car::Anova(lmer_rtd_spDh2o, type = 3)
plot(lmer_rtd_spDh2o) # No, this is not okay!

### RTD log(x+1) transformed
lmer_rtd_spDh2o_log1 <- lme4::lmer(RTD_g_cm3_log1 ~ species * Soil_depth * h2o
                                   + (1|bloc_pairs)
                                   + (1|bloc_pairs : species)
                                   + (1|bloc_pairs : Soil_depth)
                                   + (1|bloc_pairs : h2o)
                                   + (1|block_plot),
                                   data = absorptive_rtd_outRem, REML = TRUE)

summary(lmer_rtd_spDh2o_log1)
car::Anova(lmer_rtd_spDh2o_log1, type = 3)
plot(lmer_rtd_spDh2o_log1) # No, problem of increasing variance with the mean persists!

### RTD log10(x+1) transformation:
lmer_rtd_spDh2o_log10 <- lme4::lmer(RTD_g_cm3_log10 ~ species * Soil_depth * h2o
                                    + (1|bloc_pairs)
                                    + (1|bloc_pairs : species)
                                    + (1|bloc_pairs : Soil_depth)
                                    + (1|bloc_pairs : h2o)
                                    + (1|block_plot),
                                    data = absorptive_rtd_outRem, REML = TRUE)
# Oh wow, why does R not give a model singularity warning now!?
summary(lmer_rtd_spDh2o_log10)
car::Anova(lmer_rtd_spDh2o_log10, type = 3)
plot(lmer_rtd_spDh2o_log10) # Nope, the variance still increases with the mean!

### RTD sqrt(x) transformation (email 25.09.2023: Stephane Daigle wrote that the SQRT transformation was the best,
# so stick with this transformation! If there are no modes in my residuals, the sqrt transformation is okay,
# even if the variance increases slightly with the mean.):
# Use model "lmer_rtd_spDh2o_sqrt" for the 2nd paper:
lmer_rtd_spDh2o_sqrt <- lme4::lmer(RTD_g_cm3_sqrt ~ species * Soil_depth * h2o
                                   + (1|bloc_pairs)
                                   + (1|bloc_pairs : species)
                                   + (1|bloc_pairs : Soil_depth)
                                   + (1|bloc_pairs : h2o)
                                   + (1|block_plot),
                                   data = absorptive_rtd_outRem, REML = TRUE)

summary(lmer_rtd_spDh2o_sqrt)
plot(lmer_rtd_spDh2o_sqrt) # Stephane said this is okay!
anova(as_lmerModLmerTest(lmer_rtd_spDh2o_sqrt), type = 3)
summary(lmer_rtd_spDh2o_sqrt)
qqnorm(resid(lmer_rtd_spDh2o_sqrt))
qqline(resid(lmer_rtd_spDh2o_sqrt))
hist(resid(lmer_rtd_spDh2o_sqrt), breaks = 100)
car::vif(lmer_rtd_spDh2o_sqrt)
fixef(lmer_rtd_spDh2o_sqrt)
ranef(lmer_rtd_spDh2o_sqrt)
coef(lmer_rtd_spDh2o_sqrt)
VarCorr(lmer_rtd_spDh2o_sqrt)

r.squaredGLMM(lmer_rtd_spDh2o_sqrt) # R2m = 0.1971338, R2c = 0.4825986
plot(lmer_rtd_spDh2o_sqrt, form = resid(.) ~ fitted(.) | species, na.omit = TRUE) # Does not work with NA!
plot(lmer_rtd_spDh2o_sqrt, form = resid(.) ~ fitted(.) | Soil_depth) # Does not work with NA!
plot(lmer_rtd_spDh2o_sqrt, form = resid(.) ~ fitted(.) | h2o) # Does not work with NA!
plot(lmer_rtd_spDh2o_sqrt, form = RTD_g_cm3_sqrt ~ resid(.)) # Does not work with NA!

plot(ranef(lmer_rtd_spDh2o_sqrt))
lattice::qqmath(ranef(lmer_rtd_spDh2o_sqrt))

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_rtd_spDh2o_sqrt <- emmeans::emmeans(lmer_rtd_spDh2o_sqrt,
                                             specs = list(pairwise ~ h2o | species | Soil_depth),
                                             type = "response", adjust = "tukey")
# Results make sense by looking at the graph!

emmip(lmer_rtd_spDh2o_sqrt, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = RTD_g_cm3_sqrt), data = absorptive_rtd_outRem,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ tissue ~ density)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))


### boxcox transformed "rtd_g_cm3" (constant 0.001 added before to "rtd_g_cm3_bc001"
# to find the optimal lambda value):
# I accidentally used the sqrt transformed rtd for the boxcox model, which was not right.
# The sqrt transformed variable performed equally well, without the mistake of sqrt+boxcox!
lmer_rtd_spDh2o_bc <- lme4::lmer(((rtd_g_cm3 ^ 0.04698855 - 1) / 0.04698855) ~
                                   species * Soil_depth * h2o
                                 + (1|bloc_pairs)
                                 + (1|bloc_pairs : species)
                                 + (1|bloc_pairs : Soil_depth)
                                 + (1|bloc_pairs : h2o)
                                 + (1|block_plot),
                                 data = absorptive_rtd_outRem, REML = TRUE)

plot(lmer_rtd_spDh2o_bc)
summary(as_lmerModLmerTest(lmer_rtd_spDh2o_bc)) # Yes, produces results!
anova(as_lmerModLmerTest(lmer_rtd_spDh2o_bc), type = 3)
summary(lmer_rtd_spDh2o_bc)
car::Anova(lmer_rtd_spDh2o_bc, type = 3) # Impressive how robust the lmer are against violations
# from normality! All anova calls produced very similar results, with p-values only differing by
# a few 0.1 units. However, I guess Post-hoc analysis becomes more problematic!?
qqnorm(resid(lmer_rtd_spDh2o_bc))
qqline(resid(lmer_rtd_spDh2o_bc))
hist(resid(lmer_rtd_spDh2o_bc), breaks = 100)
car::vif(lmer_rtd_spDh2o_bc)
fixef(lmer_rtd_spDh2o_bc)
ranef(lmer_rtd_spDh2o_bc)
coef(lmer_rtd_spDh2o_bc)
VarCorr(lmer_rtd_spDh2o_bc)

r.squaredGLMM(lmer_rtd_spDh2o_bc) # R2m = 0.1969692, R2c = 0.5362702
plot(lmer_rtd_spDh2o_bc, form = resid(.) ~ fitted(.) | species, na.omit = TRUE) # Does not work with NA!
plot(lmer_rtd_spDh2o_bc, form = resid(.) ~ fitted(.) | Soil_depth) # Does not work with NA!
plot(lmer_rtd_spDh2o_bc, form = resid(.) ~ fitted(.) | h2o) # Does not work with NA!
plot(lmer_rtd_spDh2o_bc, form = rtd_g_cm3 ~ resid(.)) # Does not work with NA!

plot(ranef(lmer_rtd_spDh2o_bc))
lattice::qqmath(ranef(lmer_rtd_spDh2o_bc))

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_rtd_spDh2o_bc <- emmeans::emmeans(lmer_rtd_spDh2o_bc,
                                           specs = list(pairwise ~ h2o | species | Soil_depth),
                                           type = "response", adjust = "tukey")

emmip(lmer_rtd_spDh2o_bc, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = rtd_g_cm3), data = absorptive_rtd_outRem,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ tissue ~ density)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))
# Triangles (data points) all don't fit because of the boxcox transformation!

#### ------------ Restricted maximum likelihood-based mixed models (MMs) for supplementary ---------- ####

## Variables in supplementary include: 1.) SRA, 2.) Diam, 3.) Root length, 4.) RDMC, 5.) Root Volume;

### sra_cm2_g ###
lmer_sra_spDh2o <- lme4::lmer(sra_cm2_g ~ species * Soil_depth * h2o
                              + (1|bloc_pairs)
                              + (1|bloc_pairs : species)
                              + (1|bloc_pairs : Soil_depth)
                              + (1|bloc_pairs : h2o)
                              + (1|block_plot),
                              data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_sra_spDh2o) # Not good --> ChatGPT recommended a logarithmic, square root, or inverse transformation!

### sqrt(sra_cm2_g + 1) ###
lmer_sra_spDh2o_sqrt1 <- lme4::lmer(SRA_cm2_g_sqrt1 ~ species * Soil_depth * h2o
                                    + (1|bloc_pairs)
                                    + (1|bloc_pairs : species)
                                    + (1|bloc_pairs : Soil_depth)
                                    + (1|bloc_pairs : h2o)
                                    + (1|block_plot),
                                    data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_sra_spDh2o_sqrt1) # Logarithmic transformation looks slightly better!
# log(x+1) produces a more random cloud and keeps values closer together.

### log(sra_cm2_g + 1) ###
# Use "lmer_sra_spDh2o_log1" for the second paper!
lmer_sra_spDh2o_log1 <- lme4::lmer(SRA_cm2_g_log1 ~ species * Soil_depth * h2o
                                   + (1|bloc_pairs)
                                   + (1|bloc_pairs : species)
                                   + (1|bloc_pairs : Soil_depth)
                                   + (1|bloc_pairs : h2o)
                                   + (1|block_plot),
                                   data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_sra_spDh2o_log1) # Much better!

summary(as_lmerModLmerTest(lmer_sra_spDh2o_log1)) # Yes, produces results!
anova(as_lmerModLmerTest(lmer_sra_spDh2o_log1), type = 3)
qqnorm(resid(lmer_sra_spDh2o_log1))
qqline(resid(lmer_sra_spDh2o_log1)) # Looks very good!
hist(resid(lmer_sra_spDh2o_log1), breaks = 100)
car::vif(lmer_sra_spDh2o_log1) # Actualluy a good Vif this time!
r.squaredGLMM(lmer_sra_spDh2o_log1) # R2m = 0.520537, R2c = 0.5497926

fixef(lmer_sra_spDh2o_log1)
ranef(lmer_sra_spDh2o_log1)
coef(lmer_sra_spDh2o_log1)
VarCorr(lmer_sra_spDh2o_log1)
plot(lmer_sra_spDh2o_log1, form = resid(.) ~ fitted(.) | species)
plot(lmer_sra_spDh2o_log1, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_sra_spDh2o_log1, form = resid(.) ~ fitted(.) | h2o)
plot(lmer_sra_spDh2o_log1, form = SRA_cm2_g_log1 ~ resid(.))
plot(ranef(lmer_sra_spDh2o_log1))
lattice::qqmath(ranef(lmer_sra_spDh2o_log1))

### Model visualization ###
# Two-way interaction (species * h2o):
plot(allEffects(lmer_sra_spDh2o_log1), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_sra_spDh2o_log1))
# Summary statistics:
plot_model(lmer_sra_spDh2o_log1, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_sra_spDh2o_log1, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_sra_spDh2o_log1, type = "eff", terms = c("species"),
           show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_sra_spDh2o_log1, type = "eff", terms = c("h2o"),
           show.data = T) + theme_bw()
# Overall effect of water treatment per species:
emmip(lmer_sra_spDh2o_log1, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)

# Plot by h2o over species (Not too useful!):
visreg(lmer_sra_spDh2o_log1, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# A few (untrue) outliers.
# Plot by species over h2o (Better visualization!):
visreg(lmer_sra_spDh2o_log1, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_sra_spDh2o_log1 <- emmeans::emmeans(lmer_sra_spDh2o_log1,
                                             specs = list(pairwise ~ h2o | species | Soil_depth),
                                             type = "response", adjust = "tukey")

emmip(lmer_sra_spDh2o_log1, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = SRA_cm2_g_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Specific ~ root ~ area)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))
# Response directions are kind of inconclusive!

### ------------------------------------------- avgdiam_mm ---------------------------------------------- ###

### log(avgdiam_mm + 1) ###
lmer_avgdiam_spDh2o_log1 <- lme4::lmer(avgdiam_mm_log1 ~ species * Soil_depth * h2o
                                       + (1|bloc_pairs)
                                       + (1|bloc_pairs : species)
                                       + (1|bloc_pairs : Soil_depth)
                                       + (1|bloc_pairs : h2o)
                                       + (1|block_plot),
                                       data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_avgdiam_spDh2o_log1) # No change, so keep the original avgdiam_mm!

# Use "lmer_avgdiam_spDh2o" for the second paper!
lmer_avgdiam_spDh2o <- lme4::lmer(avgdiam_mm ~ species * Soil_depth * h2o
                                  + (1|bloc_pairs)
                                  + (1|bloc_pairs : species)
                                  + (1|bloc_pairs : Soil_depth)
                                  + (1|bloc_pairs : h2o)
                                  + (1|block_plot),
                                  data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_avgdiam_spDh2o) # Looks almost okay I would say, try log(x+1) nonetheless (tried above)!

summary(as_lmerModLmerTest(lmer_avgdiam_spDh2o)) # Yes, produces results!
anova(as_lmerModLmerTest(lmer_avgdiam_spDh2o), type = 3)
qqnorm(resid(lmer_avgdiam_spDh2o))
qqline(resid(lmer_avgdiam_spDh2o)) # So so!
hist(resid(lmer_avgdiam_spDh2o), breaks = 100)
car::vif(lmer_avgdiam_spDh2o) # Actualluy a good Vif this time!
r.squaredGLMM(lmer_avgdiam_spDh2o) # R2m = 0.7922844, R2c = 0.8099255

fixef(lmer_avgdiam_spDh2o)
ranef(lmer_avgdiam_spDh2o)
coef(lmer_avgdiam_spDh2o)
VarCorr(lmer_avgdiam_spDh2o)
plot(lmer_avgdiam_spDh2o, form = resid(.) ~ fitted(.) | species)
plot(lmer_avgdiam_spDh2o, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_avgdiam_spDh2o, form = resid(.) ~ fitted(.) | h2o)
plot(lmer_avgdiam_spDh2o, form = avgdiam_mm ~ resid(.))
plot(ranef(lmer_avgdiam_spDh2o))
lattice::qqmath(ranef(lmer_avgdiam_spDh2o)) # Why don't they have CI bands?
# Perhaps because I used the original avgdiam_mm?

### Model visualization ###
# Two-way interaction (species * h2o):
plot(allEffects(lmer_avgdiam_spDh2o), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_avgdiam_spDh2o))
# Summary statistics:
plot_model(lmer_avgdiam_spDh2o, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_avgdiam_spDh2o, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_avgdiam_spDh2o, type = "eff", terms = c("species"),
           show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_avgdiam_spDh2o, type = "eff", terms = c("h2o"),
           show.data = T) + theme_bw()
# Overall effect of water treatment per species:
emmip(lmer_avgdiam_spDh2o, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)

# Plot by h2o over species (Not too useful!):
visreg(lmer_avgdiam_spDh2o, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# A few (untrue) outliers.
# Plot by species over h2o (Better visualization!):
visreg(lmer_avgdiam_spDh2o, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Larix shows the greatest difference between water treatments!

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_avgdiam_spDh2o <- emmeans::emmeans(lmer_avgdiam_spDh2o,
                                            specs = list(pairwise ~ h2o | species | Soil_depth),
                                            type = "response", adjust = "tukey")

emmip(lmer_avgdiam_spDh2o, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = avgdiam_mm), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ diameter)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))

### ----------------------------------------------- Root length ----------------------------------------------- ###

### length_aff_cm ###
lmer_length_spDh2o <- lme4::lmer(length_aff_cm ~ species * Soil_depth * h2o
                                 + (1|bloc_pairs)
                                 + (1|bloc_pairs : species)
                                 + (1|bloc_pairs : Soil_depth)
                                 + (1|bloc_pairs : h2o)
                                 + (1|block_plot),
                                 data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_length_spDh2o) # Funnel shape!

### length_aff_cm_log1 ###
lmer_length_spDh2o_log1 <- lme4::lmer(length_aff_cm_log1 ~ species * Soil_depth * h2o
                                      + (1|bloc_pairs)
                                      + (1|bloc_pairs : species)
                                      + (1|bloc_pairs : Soil_depth)
                                      + (1|bloc_pairs : h2o)
                                      + (1|block_plot),
                                      data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_length_spDh2o_log1) # Funnel shape, variance decreases with increasing mean.
# Try sqrt transformation instead!

### length_aff_cm_sqrt1 ###
# Use "lmer_length_spDh2o_sqrt1" for second paper!
lmer_length_spDh2o_sqrt1 <- lme4::lmer(length_aff_cm_sqrt1 ~ species * Soil_depth * h2o
                                       + (1|bloc_pairs)
                                       + (1|bloc_pairs : species)
                                       + (1|bloc_pairs : Soil_depth)
                                       + (1|bloc_pairs : h2o)
                                       + (1|block_plot),
                                       data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_length_spDh2o_sqrt1) # Much better than logarithmic transformation!

summary(as_lmerModLmerTest(lmer_length_spDh2o_sqrt1)) # Yes, produces results!
anova(as_lmerModLmerTest(lmer_length_spDh2o_sqrt1), type = 3)
qqnorm(resid(lmer_length_spDh2o_sqrt1))
qqline(resid(lmer_length_spDh2o_sqrt1)) # Good!
hist(resid(lmer_length_spDh2o_sqrt1), breaks = 100)
car::vif(lmer_length_spDh2o_sqrt1) # Actualluy a good Vif!
r.squaredGLMM(lmer_length_spDh2o_sqrt1) # R2m = 0.6787534, R2c = 0.7538268

fixef(lmer_length_spDh2o_sqrt1)
ranef(lmer_length_spDh2o_sqrt1)
coef(lmer_length_spDh2o_sqrt1)
VarCorr(lmer_length_spDh2o_sqrt1)
plot(lmer_length_spDh2o_sqrt1, form = resid(.) ~ fitted(.) | species)
plot(lmer_length_spDh2o_sqrt1, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_length_spDh2o_sqrt1, form = resid(.) ~ fitted(.) | h2o)
plot(lmer_length_spDh2o_sqrt1, form = length_aff_cm_sqrt1 ~ resid(.))
plot(ranef(lmer_length_spDh2o_sqrt1))
lattice::qqmath(ranef(lmer_length_spDh2o_sqrt1)) # First CI band leaves zero!

### Model visualization ###
# Two-way interaction (species * h2o):
plot(allEffects(lmer_length_spDh2o_sqrt1), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Larix and Picea show significantly higher length under h2o high!
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_length_spDh2o_sqrt1))
# Summary statistics:
plot_model(lmer_length_spDh2o_sqrt1, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_length_spDh2o_sqrt1, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_length_spDh2o_sqrt1, type = "eff", terms = c("species"),
           show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_length_spDh2o_sqrt1, type = "eff", terms = c("h2o"),
           show.data = T) + theme_bw()
# Overall effect of water treatment per species:
emmip(lmer_length_spDh2o_sqrt1, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)

# Plot by h2o over species (Not too useful!):
visreg(lmer_length_spDh2o_sqrt1, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# A few (untrue) outliers.
# Plot by species over h2o (Better visualization!):
visreg(lmer_length_spDh2o_sqrt1, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Larix and Picea show the greatest difference between water treatments!

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_length_spDh2o_sqrt1 <- emmeans::emmeans(lmer_length_spDh2o_sqrt1,
                                                 specs = list(pairwise ~ h2o | species | Soil_depth),
                                                 type = "response", adjust = "tukey")

emmip(lmer_length_spDh2o_sqrt1, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = length_aff_cm_sqrt1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ length)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))

### ----------------------------------------------- RDMC ----------------------------------------------- ###

### rdmc_dry_g_fresh_g ###
lmer_rdmc_spDh2o <- lme4::lmer(rdmc_dry_g_fresh_g ~ species * Soil_depth * h2o
                               + (1|bloc_pairs)
                               + (1|bloc_pairs : species)
                               + (1|bloc_pairs : Soil_depth)
                               + (1|bloc_pairs : h2o)
                               + (1|block_plot),
                               data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_rdmc_spDh2o) # Funnel shape again, try logarithmic and sqrt transformation again!

### rdmc_dry_g_fresh_g_sqrt1 ###
lmer_rdmc_spDh2o_sqrt1 <- lme4::lmer(rdmc_dry_g_fresh_g_sqrt1 ~ species * Soil_depth * h2o
                                     + (1|bloc_pairs)
                                     + (1|bloc_pairs : species)
                                     + (1|bloc_pairs : Soil_depth)
                                     + (1|bloc_pairs : h2o)
                                     + (1|block_plot),
                                     data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_rdmc_spDh2o_sqrt1)
# Both sqrt and log transformations produced very similar results, so stick with logarithmic transformation!

### rdmc_dry_g_fresh_g_log1 ###
# Use "lmer_rdmc_spDh2o_log1" for the second paper!
lmer_rdmc_spDh2o_log1 <- lme4::lmer(rdmc_dry_g_fresh_g_log1 ~ species * Soil_depth * h2o
                                    + (1|bloc_pairs)
                                    + (1|bloc_pairs : species)
                                    + (1|bloc_pairs : Soil_depth)
                                    + (1|bloc_pairs : h2o)
                                    + (1|block_plot),
                                    data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_rdmc_spDh2o_log1) # Model residuals still have a funnel shape.
# However, the posthoc test results make sense, so keep this model and
# don't waste time with the boxcox transformation.

summary(as_lmerModLmerTest(lmer_rdmc_spDh2o_log1)) # Yes, produces results!
anova(as_lmerModLmerTest(lmer_rdmc_spDh2o_log1), type = 3)
qqnorm(resid(lmer_rdmc_spDh2o_log1))
qqline(resid(lmer_rdmc_spDh2o_log1)) # Bad, it matches the funnel shape that is still visible!
# However, I do not want to invest the time to apply a boxcox transformation to this variable at this point!
# (only supplementary)
hist(resid(lmer_rdmc_spDh2o_log1), breaks = 100)
car::vif(lmer_rdmc_spDh2o_log1) # Actualluy a good Vif!
r.squaredGLMM(lmer_rdmc_spDh2o_log1) # R2m = 0.4562386, R2c = 0.5398303 --> Still good!

fixef(lmer_rdmc_spDh2o_log1)
ranef(lmer_rdmc_spDh2o_log1)
coef(lmer_rdmc_spDh2o_log1)
VarCorr(lmer_rdmc_spDh2o_log1)
plot(lmer_rdmc_spDh2o_log1, form = resid(.) ~ fitted(.) | species)
plot(lmer_rdmc_spDh2o_log1, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_rdmc_spDh2o_log1, form = resid(.) ~ fitted(.) | h2o)
plot(lmer_rdmc_spDh2o_log1, form = rdmc_dry_g_fresh_g_log1 ~ resid(.))
plot(ranef(lmer_rdmc_spDh2o_log1))
lattice::qqmath(ranef(lmer_rdmc_spDh2o_log1)) # First CI band leaves zero!

### Model visualization ###
# Two-way interaction (species * h2o):
plot(allEffects(lmer_rdmc_spDh2o_log1), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Quercus shows consistently higher rdmc under low water!
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_rdmc_spDh2o_log1))
# Summary statistics:
plot_model(lmer_rdmc_spDh2o_log1, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_rdmc_spDh2o_log1, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Quercus shows much higher rdmc under low water!
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_rdmc_spDh2o_log1, type = "eff", terms = c("species"),
           show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_rdmc_spDh2o_log1, type = "eff", terms = c("h2o"),
           show.data = T) + theme_bw()
# Overall effect of water treatment per species:
emmip(lmer_rdmc_spDh2o_log1, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)

# Plot by h2o over species (Not too useful!):
visreg(lmer_rdmc_spDh2o_log1, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# A few (untrue) outliers.
# Plot by species over h2o (Better visualization!):
visreg(lmer_rdmc_spDh2o_log1, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Larix and Picea show the greatest difference between water treatments!

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_rdmc_spDh2o_log1 <- emmeans::emmeans(lmer_rdmc_spDh2o_log1,
                                              specs = list(pairwise ~ h2o | species | Soil_depth),
                                              type = "response", adjust = "tukey")

emmip(lmer_rdmc_spDh2o_log1, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = rdmc_dry_g_fresh_g_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ dry ~ matter ~ content)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))

### ----------------------------------------------- Root volume ----------------------------------------------- ###

### rootvol_aff_class_cm3 ###
lmer_vol_spDh2o <- lme4::lmer(rootvol_aff_class_cm3 ~ species * Soil_depth * h2o
                              + (1|bloc_pairs)
                              + (1|bloc_pairs : species)
                              + (1|bloc_pairs : Soil_depth)
                              + (1|bloc_pairs : h2o)
                              + (1|block_plot),
                              data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_vol_spDh2o) # Strong funnel shape!

lmer_vol_spDh2o_sqrt1 <- lme4::lmer(rootvol_aff_class_cm3_sqrt1 ~ species * Soil_depth * h2o
                                    + (1|bloc_pairs)
                                    + (1|bloc_pairs : species)
                                    + (1|bloc_pairs : Soil_depth)
                                    + (1|bloc_pairs : h2o)
                                    + (1|block_plot),
                                    data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_vol_spDh2o_sqrt1) # Take log transformation, it looks better!

# Use "lmer_vol_spDh2o_log1" for the first paper!
lmer_vol_spDh2o_log1 <- lme4::lmer(rootvol_aff_class_cm3_log1 ~ species * Soil_depth * h2o
                                   + (1|bloc_pairs)
                                   + (1|bloc_pairs : species)
                                   + (1|bloc_pairs : Soil_depth)
                                   + (1|bloc_pairs : h2o)
                                   + (1|block_plot),
                                   data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_vol_spDh2o_log1) # Looks good!

summary(as_lmerModLmerTest(lmer_vol_spDh2o_log1)) # Yes, produces results!
anova(as_lmerModLmerTest(lmer_vol_spDh2o_log1), type = 3)
qqnorm(resid(lmer_vol_spDh2o_log1))
qqline(resid(lmer_vol_spDh2o_log1)) # Wow, very good!
hist(resid(lmer_vol_spDh2o_log1), breaks = 100)
car::vif(lmer_vol_spDh2o_log1) # Actualluy a good Vif!
r.squaredGLMM(lmer_vol_spDh2o_log1) # R2m = 0.3670386, R2c = 0.6112955 Bloc seems to explain a lot of variation!

fixef(lmer_vol_spDh2o_log1)
ranef(lmer_vol_spDh2o_log1)
coef(lmer_vol_spDh2o_log1)
VarCorr(lmer_vol_spDh2o_log1)
plot(lmer_vol_spDh2o_log1, form = resid(.) ~ fitted(.) | species)
plot(lmer_vol_spDh2o_log1, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_vol_spDh2o_log1, form = resid(.) ~ fitted(.) | h2o)
plot(lmer_vol_spDh2o_log1, form = rootvol_aff_class_cm3_log1 ~ resid(.))
plot(ranef(lmer_vol_spDh2o_log1))
lattice::qqmath(ranef(lmer_vol_spDh2o_log1)) # First CI band leaves zero!

### Model visualization ###
# Two-way interaction (species * h2o):
plot(allEffects(lmer_vol_spDh2o_log1), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_vol_spDh2o_log1))
# Summary statistics:
plot_model(lmer_vol_spDh2o_log1, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_vol_spDh2o_log1, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_vol_spDh2o_log1, type = "eff", terms = c("species"),
           show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_vol_spDh2o_log1, type = "eff", terms = c("h2o"),
           show.data = T) + theme_bw()
# Overall effect of water treatment per species:
emmip(lmer_vol_spDh2o_log1, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)

# Plot by h2o over species (Not too useful!):
visreg(lmer_vol_spDh2o_log1, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Plot by species over h2o (Better visualization!):
visreg(lmer_vol_spDh2o_log1, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Quercus and Acer show the greatest difference between water treatments!

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_vol_spDh2o_log1 <- emmeans::emmeans(lmer_vol_spDh2o_log1,
                                             specs = list(pairwise ~ h2o | species | Soil_depth),
                                             type = "response", adjust = "tukey")

emmip(lmer_vol_spDh2o_log1, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = rootvol_aff_class_cm3_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ volume)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))


#### ------------------------------ Corrected LMMs (main variable) ----------------------------------- ####

# Daniel Schoenig advice [12.19.2023]: Hier die 3 Optionen (random intercept; random slope on main effects:
# random slope on main effects and interactions), für die beiden Modelle (Daniel's suggestions throw errors):

# Old model:
lmer_rbd_spDh2o <- lme4::lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                              + (1|bloc_pairs)
                              + (1|bloc_pairs : species)
                              + (1|bloc_pairs : Soil_depth)
                              + (1|bloc_pairs : h2o)
                              + (1|block_plot),
                              data = absorptive, REML = TRUE)

# Daniel's suggestion:
lmer_rbd_spDh2o_randInter <- lme4::lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                                        + (1 | block + block_plot),
                                        data = absorptive, REML = TRUE)
# Error: Invalid grouping factor specification, block + block_plot
# In addition: Warning message:
# In Ops.factor(block, block_plot) : ‘+’ not meaningful for factors

# Model with species as a fixed effect and random effects for block and block_plot
lmer_rbd_spDh2o_randInter <- lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                                  + (1 | block) + (1 | block:block_plot),
                                  data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
lmer_rbd_spDh2o_randInter <- lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                                  + (1 | block) + (1 | block_plot:block),
                                  data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')

# ChatGPT 4.0 correction: Assuming block_plot is nested within block
lmer_rbd_spDh2o_randInter <- lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                                  + (1 | block/block_plot),
                                  data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')

# Assuming block and block_plot are crossed random effects
lmer_rbd_spDh2o_randInter <- lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                                  + (1 | block) + (1 | block_plot),
                                  data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')

# Remove potentially repeated random terms (species and block_plot):
lmer_rbd_spDh2o_randBlock <- lme4::lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                                        + (1 | block),
                                        data = absorptive, REML = TRUE)
# This model fits without issues!!!
# GPT-4o [16.05.2024]:
# Perhaps Soil_depth does not introduce enough variability to warrant a random effect of it's own?
# Too many levels can lead to singular fits.
# If you need to explicitly account for the nesting of Soil_depth within block_plot but cannot model
# it as a random effect due to the error, ensure Soil_depth is properly accounted for as a fixed effect.
# This approach captures its influence without overwhelming the model with excessive complexity.
# However, the above model does not account for the nesting of soil depths within plots (block_plot)
# Model with nested random effects for block and block_plot
lmer_rbd_spDh2o_randInter <- lme4::lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                                        + (1 | block) + (1 | block_plot),
                                        data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')

# Daniels suggestion:
lmer_rbd_spDh2o_randInter <- lme4::lmer(ntips_length_n_cm ~ species * Soil_depth * h2o
                                        + (1 + species + Soil_depth + h2o | block + block_plot),
                                        data = absorptive, REML = TRUE)
# Daniels suggestion
lmer_rbd_spDh2o_randInter <- lme4::lmer(ntips_length_n_cm_log1 ~ species * Soil_depth * h2o
                                        + (1 + species * Soil_depth * h2o | block_plot + block_plot),
                                        data = absorptive, REML = TRUE)

summary(lmer_rbd_spDh2o_randBlock)
plot(lmer_rbd_spDh2o_randBlock) # Needs a log!
lmer_rbd_spDh2o_randBlock_log <- lme4::lmer(ntips_length_n_cm_log1 ~ species * Soil_depth * h2o
                                            + (1 | block),
                                            data = absorptive, REML = TRUE)
plot(lmer_rbd_spDh2o_randBlock_log) # Better!
summary(lmer_rbd_spDh2o_randBlock_log)
anova(as_lmerModLmerTest(lmer_rbd_spDh2o_randBlock_log), type = 3)
qqnorm(resid(lmer_rbd_spDh2o_randBlock_log))
qqline(resid(lmer_rbd_spDh2o_randBlock_log)) # That residuals leave fitted line at both ends is normal!
hist(resid(lmer_rbd_spDh2o_randBlock_log), breaks = 100) # Looks normal, except for a few outliers!
r.squaredGLMM(lmer_rbd_spDh2o_randBlock_log) # R2m = 0.4791645, R2c = 0.4847929
# (Almost the same fit as the more complex model given it';s random effect structure)
plot(lmer_rbd_spDh2o_randBlock_log, form = resid(.) ~ fitted(.) | species)
plot(lmer_rbd_spDh2o_randBlock_log, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_rbd_spDh2o_randBlock_log, form = resid(.) ~ fitted(.) | h2o)
plot(ranef(lmer_rbd_spDh2o_randBlock_log))
lattice::qqmath(ranef(lmer_rbd_spDh2o_randBlock_log)) # Looks good as long as it's overlapping zero!?

### Model visualization ###
# Two-way interaction (species * h2o):
plot(allEffects(lmer_rbd_spDh2o_randBlock_log), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_rbd_spDh2o_randBlock_log))
# Summary statistics:
plot_model(lmer_rbd_spDh2o_randBlock_log, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_rbd_spDh2o_randBlock_log, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Betula, Quercus, Pinus, and (Picea) show the strongest response!
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_rbd_spDh2o_randBlock_log, type = "eff", terms = c("species"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_rbd_spDh2o_randBlock_log, type = "eff", terms = c("h2o"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect of water treatment per species:
emmip(lmer_rbd_spDh2o_randBlock_log, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)
# Betula and Quercus show strongest response!
# Plot by h2o over species (Not too useful!):
visreg(lmer_rbd_spDh2o_randBlock_log, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Both h2o = high and low follow the same pattern!
# Plot by species over h2o (Better visualization!):
visreg(lmer_rbd_spDh2o_randBlock_log, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_rbd_spDh2o_randBlock_log <- emmeans::emmeans(lmer_rbd_spDh2o_randBlock_log,
                                                      specs = list(pairwise ~ h2o | species | Soil_depth),
                                                      type = "response", adjust = "tukey")
# Results are almost exactly the same as with the more complex model!
# Below plot is overplotted!
plot(emm_lmer_rbd_spDh2o_randBlock_log, comparisons = TRUE)

emmip(lmer_rbd_spDh2o_randBlock_log, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = ntips_length_n_cm_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ branching ~ density)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))

## After reading Popovic et al. (2024) and Smith et al. (2024), I came to the
# conclusion that the LMMs should reflect the experimental design, but
# should not involve too many random interactions. The block pairs are arbitrary
# and not reproducible --> Delete!!!

# Assuming block and block_plot are crossed random effects. I include them as
# crossed random effects because this is essentially what Stephane recommended
# before. Additionally, Alain Paquette and Stephane told me that plot is NOT
# nested in block.
# Example for nesting from: https://online.stat.psu.edu/stat503/lesson/14/14.1
# When factor B is nested in levels of factor A, the levels of the nested factor
# don't have exactly the same meaning under each level of the main factor,
# in this case factor A. In a nested design, the levels of factor (B) are not
# identical to each other at different levels of factor (A), although they might
# have the same labels. For example, if A is school and B is teacher, teacher 1
# will differ between the schools. This has to be kept in mind when trying to
# determine if the design is crossed or nested. To be crossed, the same teacher
# needs to teach at all the schools. --> In my study, the same tree species
# are present in all blocks. ==> I have a crossed, random design!

# # Use "lmer_rbdlog1_spDh2o_BlockPlot" in the second paper!
lmer_rbdlog1_spDh2o_BlockPlot <- lme4::lmer(ntips_length_n_cm_log1 ~ species * Soil_depth * h2o
                                            + (1 | block) + (1 | block_plot),
                                            data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
plot(lmer_rbdlog1_spDh2o_BlockPlot) # Variance increases slightly with the mean, but I would say it's not too bad!
summary(lmer_rbdlog1_spDh2o_BlockPlot)
anova(as_lmerModLmerTest(lmer_rbdlog1_spDh2o_BlockPlot), type = 3)
qqnorm(resid(lmer_rbdlog1_spDh2o_BlockPlot))
qqline(resid(lmer_rbdlog1_spDh2o_BlockPlot)) # That residuals leave fitted line at both ends is normal!
hist(resid(lmer_rbdlog1_spDh2o_BlockPlot), breaks = 100) # Looks normal, except for a few outliers!
r.squaredGLMM(lmer_rbdlog1_spDh2o_BlockPlot) # R2m = 0.4791645, R2c = 0.4847929
# (Exactly the same fit as the more complex models)
plot(lmer_rbdlog1_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | species)
plot(lmer_rbdlog1_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_rbdlog1_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | h2o)
plot(ranef(lmer_rbdlog1_spDh2o_BlockPlot))
lattice::qqmath(ranef(lmer_rbdlog1_spDh2o_BlockPlot)) # Looks good as long as it's overlapping zero!?

### Model visualization ###
# According to Popovic et al. (2024), model graphs should also be published
# Two-way interaction (species * h2o):
plot(allEffects(lmer_rbdlog1_spDh2o_BlockPlot), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_rbdlog1_spDh2o_BlockPlot))
# Summary statistics:
plot_model(lmer_rbdlog1_spDh2o_BlockPlot, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_rbdlog1_spDh2o_BlockPlot, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Betula, Quercus, Pinus, and (Picea) show the strongest response!
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_rbdlog1_spDh2o_BlockPlot, type = "eff", terms = c("species"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_rbdlog1_spDh2o_BlockPlot, type = "eff", terms = c("h2o"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect of water treatment per species
# This was the type of graph that Popovic et al. (2024) included. Since I am
# showing factors and not time on the x-axis, there should be NO lines.
# Why are the CIs all the same in my graph (balanced model?)?:
emmip(lmer_rbdlog1_spDh2o_BlockPlot, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)
# NOTE: Results may be misleading due to involvement in interactions
# Same species show strongest response!
# Plot by h2o over species (Not too useful!):
visreg(lmer_rbdlog1_spDh2o_BlockPlot, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Both h2o = high and low follow the same pattern!
# Plot by species over h2o (Better visualization!):
visreg(lmer_rbdlog1_spDh2o_BlockPlot, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_rbdlog1_spDh2o_BlockPlot <- emmeans::emmeans(lmer_rbdlog1_spDh2o_BlockPlot,
                                                      specs = list(pairwise ~ h2o | species | Soil_depth),
                                                      type = "response", adjust = "tukey")
# Results are almost exactly the same as with the more complex model!
# Below plot is overplotted!
plot(emm_lmer_rbdlog1_spDh2o_BlockPlot, comparisons = TRUE)

emmip(lmer_rbdlog1_spDh2o_BlockPlot, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = ntips_length_n_cm_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ branching ~ density)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))


### log(x + 1) transformed RLD_cm_cm3 ###
# Use "lmer_rldlog1_spDh2o_BlockPlot" in the second paper!
lmer_rldlog1_spDh2o_BlockPlot <- lme4::lmer(RLD_cm_cm3_log1 ~ species * Soil_depth * h2o
                                            + (1 | block) + (1 | block_plot),
                                            data = absorptive, REML = TRUE)
# Fits even WITHOUT isSingular warning!
plot(lmer_rldlog1_spDh2o_BlockPlot) # Looks okay!
summary(lmer_rldlog1_spDh2o_BlockPlot)
anova(as_lmerModLmerTest(lmer_rldlog1_spDh2o_BlockPlot), type = 3)
qqnorm(resid(lmer_rldlog1_spDh2o_BlockPlot))
qqline(resid(lmer_rldlog1_spDh2o_BlockPlot)) # That residuals leave fitted line at both ends is normal!
hist(resid(lmer_rldlog1_spDh2o_BlockPlot), breaks = 100) # Looks normal, except for a few outliers!
r.squaredGLMM(lmer_rldlog1_spDh2o_BlockPlot) # R2m = 0.6915637 (previously = 0.4791645), R2c =  0.7563447 (previously 0.4847929)
# The simpler model shows a much better fit (very good fit)
plot(lmer_rldlog1_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | species)
plot(lmer_rldlog1_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_rldlog1_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | h2o)
plot(ranef(lmer_rldlog1_spDh2o_BlockPlot))
lattice::qqmath(ranef(lmer_rldlog1_spDh2o_BlockPlot)) # One CI leaves the zero line.

### Model visualization ###
# According to Popovic et al. (2024), model graphs should also be published
# Two-way interaction (species * h2o):
plot(allEffects(lmer_rldlog1_spDh2o_BlockPlot), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_rldlog1_spDh2o_BlockPlot))
# Summary statistics:
plot_model(lmer_rldlog1_spDh2o_BlockPlot, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_rldlog1_spDh2o_BlockPlot, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Betula and Pinus show the strongest response!
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_rldlog1_spDh2o_BlockPlot, type = "eff", terms = c("species"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_rldlog1_spDh2o_BlockPlot, type = "eff", terms = c("h2o"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect of water treatment per species
emmip(lmer_rldlog1_spDh2o_BlockPlot, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)
# NOTE: Results may be misleading due to involvement in interactions
# Plot by h2o over species (Not too useful!):
visreg(lmer_rldlog1_spDh2o_BlockPlot, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Both h2o = high and low follow the same pattern!
# Plot by species over h2o (Better visualization!):
visreg(lmer_rldlog1_spDh2o_BlockPlot, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_rldlog1_spDh2o_BlockPlot <- emmeans::emmeans(lmer_rldlog1_spDh2o_BlockPlot,
                                                      specs = list(pairwise ~ h2o | species | Soil_depth),
                                                      type = "response", adjust = "tukey")
# Results are almost exactly the same as with the more complex model!
# However, some of the P-values changed slightly, so I need to adjust graphs in paper!
# Below plot is overplotted!
plot(emm_lmer_rldlog1_spDh2o_BlockPlot, comparisons = TRUE)

emmip(lmer_rldlog1_spDh2o_BlockPlot, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = RLD_cm_cm3_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ length ~ density)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))


### log(x + 1) transformed SRL_m_g1 ###
# Use "lmer_srllog1_spDh2o_BlockPlot" in the second paper!
lmer_srllog1_spDh2o_BlockPlot <- lme4::lmer(SRL_m_g_log1 ~ species * Soil_depth * h2o
                                            + (1 | block) + (1 | block_plot),
                                            data = absorptive, REML = TRUE)
# boundary (singular) fit: see help('isSingular')
# The singularity could be caused by the zero variance of block_plot
# I introduce two random intercepts into my model (block_plot and block)
# Random intercept block shows some variance that gets captured
plot(lmer_srllog1_spDh2o_BlockPlot) # No funnel shape, no modes. Stephane said that a few outliers are not a problem!
summary(lmer_srllog1_spDh2o_BlockPlot)
anova(as_lmerModLmerTest(lmer_srllog1_spDh2o_BlockPlot), type = 3)
qqnorm(resid(lmer_srllog1_spDh2o_BlockPlot))
qqline(resid(lmer_srllog1_spDh2o_BlockPlot)) # That residuals leave fitted line at both ends is normal!
hist(resid(lmer_srllog1_spDh2o_BlockPlot), breaks = 100) # Looks normal, except for a few outliers!
r.squaredGLMM(lmer_srllog1_spDh2o_BlockPlot) # R2m = 0.5543202, R2c = 0.5765576
plot(lmer_srllog1_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | species)
plot(lmer_srllog1_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_srllog1_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | h2o)
plot(ranef(lmer_srllog1_spDh2o_BlockPlot))
lattice::qqmath(ranef(lmer_srllog1_spDh2o_BlockPlot)) # One CI leaves the zero line.

### Model visualization ###
# According to Popovic et al. (2024), model graphs should also be published
# Two-way interaction (species * h2o):
plot(allEffects(lmer_srllog1_spDh2o_BlockPlot), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_srllog1_spDh2o_BlockPlot))
# Summary statistics:
plot_model(lmer_srllog1_spDh2o_BlockPlot, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_srllog1_spDh2o_BlockPlot, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Quercus and Pinus show the strongest response!
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_srllog1_spDh2o_BlockPlot, type = "eff", terms = c("species"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_srllog1_spDh2o_BlockPlot, type = "eff", terms = c("h2o"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect of water treatment per species
emmip(lmer_srllog1_spDh2o_BlockPlot, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)
# NOTE: Results may be misleading due to involvement in interactions
# Plot by h2o over species (Not too useful!):
visreg(lmer_srllog1_spDh2o_BlockPlot, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Both h2o = high and low follow the same pattern!
# Plot by species over h2o (Better visualization!):
visreg(lmer_srllog1_spDh2o_BlockPlot, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_srllog1_spDh2o_BlockPlot <- emmeans::emmeans(lmer_srllog1_spDh2o_BlockPlot,
                                                      specs = list(pairwise ~ h2o | species | Soil_depth),
                                                      type = "response", adjust = "tukey")
# Results are almost exactly the same as with the more complex model!
# Below plot is overplotted!
plot(emm_lmer_srllog1_spDh2o_BlockPlot, comparisons = TRUE)

emmip(lmer_srllog1_spDh2o_BlockPlot, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = SRL_m_g_log1), data = absorptive,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Specific ~ root ~ length)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))


### sqrt() transformed RTD_m_g ###
# Use "lmer_rtdsqrt_spDh2o_BlockPlot" in the second paper!
lmer_rtdsqrt_spDh2o_BlockPlot <- lme4::lmer(RTD_g_cm3_sqrt ~ species * Soil_depth * h2o
                                            + (1 | block) + (1 | block_plot),
                                            data = absorptive_rtd_outRem, REML = TRUE)
# Wow, model fits without problems!
plot(lmer_rtdsqrt_spDh2o_BlockPlot) # Except for four outliers, it looks okay!
summary(lmer_rtdsqrt_spDh2o_BlockPlot)
anova(as_lmerModLmerTest(lmer_rtdsqrt_spDh2o_BlockPlot), type = 3)
qqnorm(resid(lmer_rtdsqrt_spDh2o_BlockPlot))
qqline(resid(lmer_rtdsqrt_spDh2o_BlockPlot)) # That residuals leave fitted line at both ends is normal!
hist(resid(lmer_rtdsqrt_spDh2o_BlockPlot), breaks = 100) # Looks normal, except for a few outliers!
r.squaredGLMM(lmer_rtdsqrt_spDh2o_BlockPlot) # R2m = 0.1973886, R2c = 0.4825925
plot(lmer_rtdsqrt_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | species)
# Error in `[[<-.data.frame`(`*tmp*`, j, value = c(2L, 2L, 2L, 2L, 2L, 4L,  : 
# replacement has 240 rows, data has 239
# This error comes from the fact that I removed one extreme outlier.
plot(lmer_rtdsqrt_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | Soil_depth)
plot(lmer_rtdsqrt_spDh2o_BlockPlot, form = resid(.) ~ fitted(.) | h2o)
plot(ranef(lmer_rtdsqrt_spDh2o_BlockPlot))
lattice::qqmath(ranef(lmer_rtdsqrt_spDh2o_BlockPlot)) # Four CIs leave the zero line.

### Model visualization ###
# According to Popovic et al. (2024), model graphs should also be published
# Two-way interaction (species * h2o):
plot(allEffects(lmer_rtdsqrt_spDh2o_BlockPlot), multiline = T, rug = F, ci.style = "line",
     show.data = T) + theme_bw()
# Water treatment (High vs. Low) effect per species:
plot(predictorEffect(predictor = "h2o", mod = lmer_rtdsqrt_spDh2o_BlockPlot))
# Summary statistics:
plot_model(lmer_rtdsqrt_spDh2o_BlockPlot, show.values = T, sort.est = NULL,
           vline.color = "grey", value.offset = -0.3)
# Two-way interaction divided in facets H2O = High vs. Low:
plot_model(lmer_rtdsqrt_spDh2o_BlockPlot, type = "eff", terms = c("species", "h2o"),
           show.data = F, dodge = 0.5, colors = c("blue", "red")) + theme_bw()
# Only Quercus shows somewhat of a response!
# Overall effect for species (not too useful, but good to see data points):
plot_model(lmer_rtdsqrt_spDh2o_BlockPlot, type = "eff", terms = c("species"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect for h2o (not too useful, but good to see data points):
plot_model(lmer_rtdsqrt_spDh2o_BlockPlot, type = "eff", terms = c("h2o"),
           jitter = 1, show.data = T) + theme_bw()
# Overall effect of water treatment per species
emmip(lmer_rtdsqrt_spDh2o_BlockPlot, species ~ h2o , cov.reduce = range, CIs = T, dodge = 0.5)
# NOTE: Results may be misleading due to involvement in interactions
# Something is definitely happening with Quercus RTD!
# Plot by h2o over species (Not too useful!):
visreg(lmer_rtdsqrt_spDh2o_BlockPlot, 
       xvar = "species",
       by = "h2o",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)
# Both h2o = high and low follow the same pattern!
# Plot by species over h2o (Better visualization!):
visreg(lmer_rtdsqrt_spDh2o_BlockPlot, 
       xvar = "h2o",
       by = "species",
       breaks = 3,
       type = "conditional",
       data = absorptive,
       alpha = 0.05,
       nn = 101,
       jitter = FALSE,
       plot = TRUE)

### emmeans (Tukey’s Honest Significant Difference Test) ###

emm_lmer_rtdsqrt_spDh2o_BlockPlot <- emmeans::emmeans(lmer_rtdsqrt_spDh2o_BlockPlot,
                                                      specs = list(pairwise ~ h2o | species | Soil_depth),
                                                      type = "response", adjust = "tukey")
# Results are almost exactly the same as with the more complex model!
# However, some of the P-values changed slightly, so I need to adjust graphs in paper!
# Below plot is overplotted!
plot(emm_lmer_rtdsqrt_spDh2o_BlockPlot, comparisons = TRUE)

emmip(lmer_rtdsqrt_spDh2o_BlockPlot, ~ h2o | species | Soil_depth, CIs = TRUE, type = "response") +
  geom_point(aes(x = h2o, y = RTD_g_cm3_sqrt), data = absorptive_rtd_outRem,
             pch = 2, color = "blue") +
  labs(x = expression(H[2]*O),
       y = expression(Root ~ tissue ~ density)) +
  labs(colour = expression(H[2]*O)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"))
# Warning message:
# Removed 1 row containing missing values or values outside the scale range (`geom_point()`).
# This warning occurs because I removed one outlier (data.frame = absorptive_rtd_outRem).

#### --------------------------- Corrected LMMs (supplementary variables) --------------------------- ####






#### ---------------------- Percentages increase, decrease, or difference ---------------------- ####

### RBD (n cm-1) ###
# New percentages after laptop was stolen (stayed mostly the same):
mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Betula papyrifera"
                                        & absorptive$h2o == "High" &
                                          absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 1.1075 n cm-1
mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Betula papyrifera"
                                        & absorptive$h2o == "Low" &
                                          absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 2.4125 n cm-1
(2.4125 - 1.1075) / 1.1075 * 100 # 117.833 % Increase from High to Low (Don't round!).

mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Betula papyrifera"
                                        & absorptive$h2o == "High" &
                                          absorptive$Soil_depth == "5-10")], na.rm = TRUE) # 0.65 n cm-1
mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Betula papyrifera"
                                        & absorptive$h2o == "Low" &
                                          absorptive$Soil_depth == "5-10")], na.rm = TRUE) # 1.3975 n cm-1
(1.3975 - 0.65) / 0.65 * 100 # 115 % Increase from High to Low (Don't round!).

mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Quercus rubra"
                                        & absorptive$h2o == "High" &
                                          absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 0.5625 n cm-1
mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Quercus rubra"
                                        & absorptive$h2o == "Low" &
                                          absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 1.89 n cm-1
(1.89 - 0.5625) / 0.5625 * 100 # 236 % Increase from High to Low (Don't round!).

mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Larix laricina"
                                        & absorptive$h2o == "High" &
                                          absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 3.695 n cm-1
mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Larix laricina"
                                        & absorptive$h2o == "Low" &
                                          absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 1.34 n cm-1
(1.34 - 3.695) / 3.695 * 100 # -63.73478 % Decrease from High to Low!

mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Pinus strobus"
                                        & absorptive$h2o == "High" &
                                          absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 0.1925 n cm-1
mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Pinus strobus"
                                        & absorptive$h2o == "Low" &
                                          absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 0.87 n cm-1
(0.87 - 0.1925) / 0.1925 * 100 # 351.9481 % Increase from High to Low (Don't round!).

mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Pinus strobus"
                                        & absorptive$h2o == "High" &
                                          absorptive$Soil_depth == "5-10")], na.rm = TRUE) # 0.54 n cm-1
mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Pinus strobus"
                                        & absorptive$h2o == "Low" &
                                          absorptive$Soil_depth == "5-10")], na.rm = TRUE) # 1.585 n cm-1
(1.585 - 0.54) / 0.54 * 100 # 193.5185 % Increase from High to Low (Don't round!).

mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Picea glauca"
                                        & absorptive$h2o == "High" &
                                          absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 0.2775 n cm-1
mean(absorptive$ntips_length_n_cm[which(absorptive$species == "Picea glauca"
                                        & absorptive$h2o == "Low" &
                                          absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 1.55 n cm-1
(1.55 - 0.2775) / 0.2775 * 100 # 458.5586 % Increase from High to Low (Don't round!).

### Root length density (cm cm-3) ###
# Only the Species x Soil depth interaction was significant! --> Olivier Villemaire-Cote:
# If double interactions are significant, don't interpret single fixed effects!
mean(absorptive$RLD_cm_cm3[which(absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 39.20611 cm cm-3
mean(absorptive$RLD_cm_cm3[which(absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 4.674809 cm cm-3
(4.674809 - 39.20611) / 39.20611 * 100 # -88.07633 % Decrease from 0-5 cm to 20-30 cm soil depth!

# [20.11.2023] Christian comment: Clear trend of lower RLD from high- to low-water in the 0-5 cm soil depth:
mean(absorptive$RLD_cm_cm3[which(absorptive$Soil_depth == "0-5" & absorptive$h2o == "High")], na.rm = TRUE) # 48.8465 cm cm-3
mean(absorptive$RLD_cm_cm3[which(absorptive$Soil_depth == "0-5" & absorptive$h2o == "Low")], na.rm = TRUE) # 29.56572 cm cm-3
(29.56572 - 48.8465) / 48.8465 * 100 # -39.47218 % Decrease from High to low water in the 0-5 cm soil depth!

mean(absorptive$RLD_cm_cm3[which(absorptive$species == "Larix laricina"
                                 & absorptive$h2o == "High" &
                                   absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 44.94575 cm cm-3
mean(absorptive$RLD_cm_cm3[which(absorptive$species == "Larix laricina"
                                 & absorptive$h2o == "Low" &
                                   absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 21.05451 cm cm-3
(21.05451 - 44.94575) / 44.94575 * 100 # -53.15573 % Decrease from High to Low!

mean(absorptive$RLD_cm_cm3[which(absorptive$species == "Larix laricina"
                                 & absorptive$h2o == "High" &
                                   absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 1.398042 cm cm-3
mean(absorptive$RLD_cm_cm3[which(absorptive$species == "Larix laricina"
                                 & absorptive$h2o == "Low" &
                                   absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 3.789398 cm cm-3
(3.789398 - 1.398042) / 1.398042 * 100 # 171.0504 % Increase from High to Low (Don't round!).

mean(absorptive$RLD_cm_cm3[which(absorptive$species == "Pinus strobus"
                                 & absorptive$h2o == "High" &
                                   absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 4.716085 cm cm-3
mean(absorptive$RLD_cm_cm3[which(absorptive$species == "Pinus strobus"
                                 & absorptive$h2o == "Low" &
                                   absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 1.558565 cm cm-3
(1.558565 - 4.716085) / 4.716085 * 100 # -66.95214 % Decrease from High to Low!

mean(absorptive$RLD_cm_cm3[which(absorptive$species == "Picea glauca"
                                 & absorptive$h2o == "High" &
                                   absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 60.76597 cm cm-3
mean(absorptive$RLD_cm_cm3[which(absorptive$species == "Picea glauca"
                                 & absorptive$h2o == "Low" &
                                   absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 28.13935 cm cm-3
(28.13935 - 60.76597) / 60.76597 * 100 # -53.69226 % Decrease from High to Low!

### Specific root length (m g-1) ###
mean(absorptive$srl_m_g[which(absorptive$h2o == "High")], na.rm = TRUE) # 62.14167 m g-1
mean(absorptive$srl_m_g[which(absorptive$h2o == "Low")], na.rm = TRUE) # 54.61 m g-1
(54.61 - 62.14167) / 62.14167 * 100 # -12.12016 % Decrease from High to Low!

mean(absorptive$srl_m_g[which(absorptive$species == "Larix laricina"
                              & absorptive$h2o == "High" &
                                absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 42.65 m g-1
mean(absorptive$srl_m_g[which(absorptive$species == "Larix laricina"
                                 & absorptive$h2o == "Low" &
                                   absorptive$Soil_depth == "0-5")], na.rm = TRUE) # 23.75 m g-1
(23.75 - 42.65) / 42.65 * 100 # -44.31419 % Decrease from High to Low!

mean(absorptive$srl_m_g[which(absorptive$species == "Larix laricina"
                              & absorptive$h2o == "High" &
                                absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 26.15 m g-1
mean(absorptive$srl_m_g[which(absorptive$species == "Larix laricina"
                              & absorptive$h2o == "Low" &
                                absorptive$Soil_depth == "20-30")], na.rm = TRUE) # 32.65 m g-1
(32.65 - 26.15) / 26.15 * 100 # 24.8566 % Increase from High to Low (Don't round!).

### Root tissue density (g cm-3) ###
mean(absorptive_rtd_outRem$rtd_g_cm3[which(absorptive_rtd_outRem$Soil_depth == "0-5")],
     na.rm = TRUE) # 0.3277872 g cm-3
mean(absorptive_rtd_outRem$rtd_g_cm3[which(absorptive_rtd_outRem$Soil_depth == "20-30")],
     na.rm = TRUE) # 0.1697917 g cm-3
(0.1697917 - 0.3277872) / 0.3277872 * 100 # -48.20063 % Decrease from 0-5 cm to 20-30 cm soil depth!

mean(absorptive_rtd_outRem$rtd_g_cm3[which(absorptive_rtd_outRem$species == "Larix laricina"
                                           & absorptive_rtd_outRem$h2o == "High" &
                                             absorptive_rtd_outRem$Soil_depth == "0-5")], na.rm = TRUE) # 0.2335 g cm-3
mean(absorptive_rtd_outRem$rtd_g_cm3[which(absorptive_rtd_outRem$species == "Larix laricina"
                                           & absorptive_rtd_outRem$h2o == "Low" &
                                             absorptive_rtd_outRem$Soil_depth == "0-5")], na.rm = TRUE) # 0.10125 g cm-3
(0.10125 - 0.2335) / 0.2335 * 100 # -56.63812 % Decrease from High to Low!

mean(absorptive_rtd_outRem$rtd_g_cm3[which(absorptive_rtd_outRem$species == "Quercus rubra"
                                           & absorptive_rtd_outRem$h2o == "High" &
                                             absorptive_rtd_outRem$Soil_depth == "5-10")], na.rm = TRUE)
# 0.115 g cm-3
mean(absorptive_rtd_outRem$rtd_g_cm3[which(absorptive_rtd_outRem$species == "Quercus rubra"
                                           & absorptive_rtd_outRem$h2o == "Low" &
                                             absorptive_rtd_outRem$Soil_depth == "5-10")], na.rm = TRUE)
# 0.6615 g cm-3
(0.6615 - 0.115) / 0.115 * 100 # 475.2174 % Increase from High to Low.







