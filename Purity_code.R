## setup environment ----
setwd("E:/species")
load('Purity/puritysession.RData')
## install the latest release of biomod2
#devtools::install_github('biomodhub/biomod2')

## load the required packages
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(randomForest)
library(raster)
library(rgdal)
library(usdm)
library(sf)
library(tidyverse)
library(rgeos)
library(scales)
library(fasterize)
library(rJava)
library(dismo)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
# Obtaining species data from GBIF
set.seed(100)
mysp <- read.csv("Purity/Purity/occurrence.csv")
glimpse(mysp)
myspecies <- mysp%>%
  dplyr::select(species, Longitude , Latitude)

unique(myspecies$species)

#Loading predictor variables
lst <- list.files(path="data_kenya/",pattern='asc$',all.files = TRUE, full.names = T) 
preds<-stack(lst)

plot(preds[[4]])
###################Choose background points ##################
length(which(!is.na(values(subset(preds,1)))))

#####Select the variables to use###############################################
poly<-readOGR("data/Kenya/Kenya.shp")

spg<- myspecies%>%dplyr::select(Longitude,Latitude)
head(spg)
spg$vulture<- 1
spg<- spg%>%drop_na()
nrow(spg)
coordinates(spg)<-c('Longitude', 'Latitude')
class(spg)
head(spg)

crs(spg)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

spg<-spg[poly,] 

plot(spg,add=T)

ex<-raster ::extract (preds,spg)
head(ex)

v<-usdm::vifstep(ex)
preds<-usdm::exclude(preds, v)

spg<-as.data.frame(spg)
#write.csv(spg, "Carol/final.csv")

preds <-
  stack(
    c(
      Aspect  = "data_kenya/aspect_output.asc",
      Elevation  = "data_kenya/elev_output.asc",
      HII= 'data_kenya/hpi_output.asc',
      NDVI = 'data_kenya/ndvi_output.asc',
      Slope = 'data_kenya/slope_output.asc',
      bio_13 = 'data_kenya/wc2.1_30s_bio_13.asc',
      bio_14 = 'data_kenya/wc2.1_30s_bio_14.asc',
      bio_18 = 'data_kenya/wc2.1_30s_bio_18.asc',
      bio_19 = 'data_kenya/wc2.1_30s_bio_19.asc',
      bio_3 = 'data_kenya/wc2.1_30s_bio_3.asc',
      bio_7 = 'data_kenya/wc2.1_30s_bio_7.asc'
    )
  )


## format the data ----
rp_data <- 
  BIOMOD_FormatingData(
    resp.var = spg['vulture'],
    resp.xy = spg[, c('Longitude', 'Latitude')],
    expl.var = preds,
    resp.name = "vulture",
    PA.nb.rep = 1,
    PA.nb.absences = 10000,
    PA.strategy = 'random'
  )

## formatted object summary
rp_data

## plot of selected pseudo-absences
plot(rp_data)

BIOMOD_ModelingOptions()

## define individual models options ---- 
rp_opt <- 
  BIOMOD_ModelingOptions(
    MAXENT = list(path_to_maxent.jar = 'E:/species/maxent/maxent.jar', maximumiterations = 1000)
  )

## run the individual models ----

rp_models <- BIOMOD_Modeling(bm.format = rp_data,
                                    modeling.id = 'AllModels',
                                    models = c('MAXENT'),
                                    bm.options = rp_opt,
                                    nb.rep = 2,
                                    data.split.perc = 70,
                                    metric.eval = c('TSS','ROC'),
                                    var.import = 3,
                                    do.full.models = FALSE,
                                    seed.val = 42)

#myLoadedModels <- BIOMOD_LoadModels("Raptors/Raptors.demo1ensemble.models.out", models = c("GLM", "GBM", "RF", "GAM", "ANN", "MARS", "MAXENT.Phillips"))

## asses individual models quality ----

myGLMs <- BIOMOD_LoadModels(bm.out = rp_models, algo = c('MAXENT'))

mods <- get_built_models(myBiomodModelOut, run = 'RUN1')

myRespPlot2D <- 
  bm_PlotResponseCurves(
    bm.out  = rp_models,
    models.chosen = myGLMs,
    xed.var = 'mean')
    
 
y<-myRespPlot2D$tab

df.wide <- pivot_wider(y, names_from = expl.name, values_from = c(expl.val, pred.val)) 

###########################################for each variable\\\\\\\\\\\\\\\\\\
Aspect<- df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_Aspect, pred.val_Aspect, species)%>%
  group_by(id)%>%
  summarise(expl.val_Aspect = mean(expl.val_Aspect), pred.val_Aspect = mean(pred.val_Aspect))
  
 glimpse(Aspect)

p1<-ggplot(Aspect, aes(x = expl.val_Aspect, y = pred.val_Aspect), color= pred.name) + 
  geom_line() +
  xlab("Aspect (º)") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))

Elevation<-df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_Elevation, pred.val_Elevation, species)%>%
  group_by(id)%>%
  summarise(expl.val_Elevation = mean(expl.val_Elevation), pred.val_Elevation = mean(pred.val_Elevation))
  
glimpse(Elevation)

p2<-ggplot(Elevation, aes(x = expl.val_Elevation, y = pred.val_Elevation), color= pred.name) + 
  geom_line() +
  xlab("Elevation (m)") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))

HII<- df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_HII, pred.val_HII, species)%>%
  group_by(id)%>%
  summarise(expl.val_HII = mean(expl.val_HII), pred.val_HII = mean(pred.val_HII))


p3<-ggplot(HII, aes(x = expl.val_HII, y = pred.val_HII)) + 
  geom_line() +
  xlab("Human Influence Index") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))


NDVI<- df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_NDVI, pred.val_NDVI, species)%>%
  group_by(id)%>%
  summarise(expl.val_NDVI = mean(expl.val_NDVI), pred.val_NDVI = mean(pred.val_NDVI))


p4<-ggplot(NDVI, aes(x = expl.val_NDVI, y = pred.val_NDVI)) + 
  geom_line() +
  xlab("Normalized Difference Vegetation Index") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))


bio_13<- df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_bio_13, pred.val_bio_13, species)%>%
  group_by(id)%>%
  summarise(expl.val_bio_13 = mean(expl.val_bio_13), pred.val_bio_13 = mean(pred.val_bio_13))


p5<-ggplot(bio_13, aes(expl.val_bio_13,pred.val_bio_13)) + 
  geom_line() +
  xlab("Precipitation of Wettest Month (mm)") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))


bio_14<- df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_bio_14, pred.val_bio_14, species)%>%
  group_by(id)%>%
  summarise(expl.val_bio_14 = mean(expl.val_bio_14), pred.val_bio_14 = mean(pred.val_bio_14))

p6<-ggplot(bio_14, aes(expl.val_bio_14,pred.val_bio_14)) + 
  geom_line() +
  xlab("Precipitation of Driest Month (Mm)") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))

bio_18<- df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_bio_18, pred.val_bio_18, species)%>%
  group_by(id)%>%
  summarise(expl.val_bio_18 = mean(expl.val_bio_18), pred.val_bio_18 = mean(pred.val_bio_18))

p7<-ggplot(bio_18, aes(expl.val_bio_18,pred.val_bio_18)) + 
  geom_line() +
  xlab("Precipitation of Warmest Quarter (mm)") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))

bio_19<- df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_bio_19, pred.val_bio_19, species)%>%
  group_by(id)%>%
  summarise(expl.val_bio_19 = mean(expl.val_bio_19), pred.val_bio_19 = mean(pred.val_bio_19))

p8<-ggplot(bio_19, aes(expl.val_bio_19,pred.val_bio_19)) + 
  geom_line() +
  xlab("Precipitation of coldest Quarter (mm) ") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))


bio_3<- df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_bio_3, pred.val_bio_3, species)%>%
  group_by(id)%>%
  summarise(expl.val_bio_3 = mean(expl.val_bio_3), pred.val_bio_3 = mean(pred.val_bio_3))

p9<-ggplot(bio_3, aes(expl.val_bio_3,pred.val_bio_3)) + 
  geom_line() +
  xlab("Isothermality (%)") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))

bio_7<-df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_bio_7, pred.val_bio_7, species)%>%
  group_by(id)%>%
  summarise(expl.val_bio_7 = mean(expl.val_bio_7), pred.val_bio_7 = mean(pred.val_bio_7))

p10<-ggplot(bio_7, aes(expl.val_bio_7,pred.val_bio_7)) + 
  geom_line() +
  xlab("Temperature Annual range (ºC)") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))

slope<- df.wide %>%
  mutate(species = word(pred.name, sep = "\\_"))%>%
  dplyr::select(id, expl.val_Slope, pred.val_Slope, species)%>%
  group_by(id)%>%
  summarise(expl.val_Slope = mean(expl.val_Slope), pred.val_Slope = mean(pred.val_Slope))

p11<-ggplot(slope, aes(expl.val_Slope,pred.val_Slope)) + 
  geom_line() +
  xlab("Slope (º)") + ylab("Probability of presence") + theme(plot.title = element_text(hjust = 0.5))



cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, labels = "auto")

glimpse(df.wide)


########################################################################################

## get models evaluation scores
rp_models_scores <- get_evaluations(rp_models)


#write.csv(evaluation_table, "Carol/evaluation_table.csv")

## ProLau_models_scores is a 5 dimension array containing the scores of the models
dim(rp_models_scores)
dimnames(rp_models_scores)


# Represent mean evaluation scores
bm_PlotEvalMean(bm.out = rp_models)

## check variable importance
(rp_models_var_import <- get_variables_importance(rp_models))
varibleimp<-tibble::rownames_to_column(as.data.frame(rp_models_var_import), "Variable")

varibleimp1 <- varibleimp %>%
  rowwise()%>%
  group_by(expl.var)%>%
  summarise(m = mean(var.imp), sd = sd(var.imp))%>%
  mutate( se=sd/sqrt(3))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, 3-1))

# Standard Error
ggplot(varibleimp1) +
  geom_bar( aes(x=expl.var, y=m), stat="identity", fill="grey", alpha=0.5) +
  geom_errorbar( aes(x=expl.var, ymin=m-se, ymax=m+se), width=0.4, colour="black", alpha=0.9, linewidth=1.5) +
  ggtitle("Variable Importance")+ xlab("Variables") + ylab("Mean")


## individual models response plots
#rp_glm <- BIOMOD_LoadModels(rp_models, models='GLM')
#rp_gbm <- BIOMOD_LoadModels(rp_models, models='GBM')
#rp_rf <- BIOMOD_LoadModels(rp_models, models='RF')
#rp_gam <- BIOMOD_LoadModels(rp_models, models='GAM')
# rp_max <- BIOMOD_LoadModels(rp_models, models='MAXENT.Phillips')

#glm_eval_strip <- 
#biomod2::response.plot2(
#   models  = rp_glm,
#   Data = get_formal_data(rp_models,'expl.var'), 
#   show.variables= get_formal_data(rp_models,'expl.var.names'),
#   do.bivariate = FALSE,
#   fixed.var.metric = 'median',
#   legend = FALSE,
#   display_title = FALSE,
#   data_species = get_formal_data(rp_models,'resp.var')
# )

# gam_eval_strip <- 
# biomod2::response.plot2(
#  models  = rp_max,
# Data = get_formal_data(rp_models,'expl.var'), 
# show.variables= get_formal_data(rp_models,'expl.var.names'),
# do.bivariate = FALSE,
# fixed.var.metric = 'median',
# legend = FALSE,
#  display_title = FALSE,
#  data_species = get_formal_data(rp_models,'resp.var')
#  )

# Model ensemble models
rp_ensemble_models  <- BIOMOD_EnsembleModeling(bm.mod = rp_models,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c('EMmean', 'EMca'),
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.7),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3,
                                      seed.val = 42)
rp_ensemble_models 


## asses ensemble models quality ----
(rp_ensemble_models_scores <- get_evaluations(rp_ensemble_models))

# Project single models
file.proj <- paste0("vulture", "vulture/proj_current/", "vulture", ".Current.projection.out")
if (file.exists(file.proj)) {
  myBiomodProj <- get(load(file.proj))
} else {
  myBiomodProj <- BIOMOD_Projection(bm.mod = rp_models,
                                    proj.name = 'Current',
                                    new.env = preds,
                                    models.chosen = 'all')
}
myBiomodProj
plot(myBiomodProj)

rp_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(bm.em = rp_ensemble_models , 
                                                              bm.proj = myBiomodProj,
                                                              models.chosen = 'all',
                                                              metric.binary = 'all',
                                                              metric.filter = 'all')


##############################################################################################################################3


## future projections

## load 2040 bioclim variables
bioclim_ZA_2040_BC45 <-
  stack(
    c(
      Aspect  = "data_kenya/aspect_output.asc",
      Elevation  = "data_kenya/elev_output.asc",
      HII= 'data_kenya/hpi_output.asc',
      NDVI = 'data_kenya/ndvi_output.asc',
      Slope = 'data_kenya/slope_output.asc',
      bio_13 = 'Future/ssp245/2140/HadGEM3.GC31.LL_ssp245_2021.2040.13.asc',
      bio_14 = 'Future/ssp245/2140/HadGEM3.GC31.LL_ssp245_2021.2040.14.asc',
      bio_18 = 'Future/ssp245/2140/HadGEM3.GC31.LL_ssp245_2021.2040.18.asc',
      bio_19 = 'Future/ssp245/2140/HadGEM3.GC31.LL_ssp245_2021.2040.19.asc',
      bio_3 = 'Future/ssp245/2140/HadGEM3.GC31.LL_ssp245_2021.2040.3.asc',
      bio_7 = 'Future/ssp245/2140/HadGEM3.GC31.LL_ssp245_2021.2040.7.asc'
    )
  )


rp_models_proj_2040_BC45 <-BIOMOD_Projection(bm.mod =rp_models,
                  proj.name = 'Future204045',
                  new.env = bioclim_ZA_2040_BC45,
                  models.chosen = 'all',
                  build.clamping.mask = TRUE)


rp_ensemble_models_proj_2040_BC45<- BIOMOD_EnsembleForecasting(bm.em = rp_ensemble_models , 
                                                              bm.proj = rp_models_proj_2040_BC45,
                                                              models.chosen = 'all',
                                                              metric.binary = 'all',
                                                              metric.filter = 'all')

## load 2040 BC85bioclim variables
bioclim_ZA_2040_BC85 <-
  stack(
    c(
      Aspect  = "data_kenya/aspect_output.asc",
      Elevation  = "data_kenya/elev_output.asc",
      HII= 'data_kenya/hpi_output.asc',
      NDVI = 'data_kenya/ndvi_output.asc',
      Slope = 'data_kenya/slope_output.asc',
      bio_13 = 'Future/ssp285/2140/HadGEM3.GC31.LL_ssp585_2021.2040.13.asc',
      bio_14 = 'Future/ssp285/2140/HadGEM3.GC31.LL_ssp585_2021.2040.14.asc',
      bio_18 = 'Future/ssp285/2140/HadGEM3.GC31.LL_ssp585_2021.2040.18.asc',
      bio_19 = 'Future/ssp285/2140/HadGEM3.GC31.LL_ssp585_2021.2040.19.asc',
      bio_3 = 'Future/ssp285/2140/HadGEM3.GC31.LL_ssp585_2021.2040.3.asc',
      bio_7 = 'Future/ssp285/2140/HadGEM3.GC31.LL_ssp585_2021.2040.7.asc'
    )
  )

rp_models_proj_2040_BC85 <-BIOMOD_Projection(bm.mod =rp_models,
                                             proj.name = 'Future204085',
                                             new.env = bioclim_ZA_2040_BC85,
                                             models.chosen = 'all',
                                             build.clamping.mask = TRUE)


rp_ensemble_models_proj_2040_BC85<- BIOMOD_EnsembleForecasting(bm.em = rp_ensemble_models , 
                                                               bm.proj = rp_models_proj_2040_BC85,
                                                               models.chosen = 'all',
                                                               metric.binary = 'all',
                                                               metric.filter = 'all')


#########################################################################################################################

save.image(file='Purity/puritysession.RData')

