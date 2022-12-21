####ALLISON DENNERT, Cross-watershed salmon lead traits analyses
####"MANUSCRIPT TITLE" by A. Dennert, E. Elle, J. Reynolds

#This R code analyzes/visualizes 1) isotopes, 2) % nitrogen 3) leaf mass per area,
#4) leaf area, and 5) leaf greenness in salmonberry (RUSP), false lily-of-the-valley
#(MADI), false azalea (MEFE), blueberry species (VASP), and foamflower (TITR).
#These data were collected across 14 streams of varying chum and pink salmon 
#spawning densities.

####0. PACKAGE LOADING####
library(dplyr)
library(lme4)
library(visreg)
library(ggplot2)
library(tidyr)
library(car)
library(MASS)
library(glmmTMB)
library(DHARMa)
library(PNWColors)

####00. READ, CLEAN, TRANSFORM DATA####
data <- read.csv("Raw Data/final_data.csv", strip.white = TRUE)
str(data)
data$species <- as.factor(data$species)
data$stream <- as.factor(data$stream)
data$quadrat.id <- as.factor(data$quadrat.id)

#calculate "northness" and "eastness" from aspect data, as data are circular;
#sin() and cos() functions operate with radians, so we must convert aspect from
#degrees to radians and then apply the trigonometric transformations
data <- data %>% mutate(northness = cos((aspect*(pi/180)))) %>% 
  mutate(eastness = sin((aspect*(pi/180))))

#calculate relative soil moisture based on soil moisture reference plot
data <- data %>% mutate(rel.soil.moisture = avg.soil.moisture - ref.soil.moisture)

#calculate %N and C:N ratios; must convert ug to mg
data <- data %>% mutate(percent.N = ((0.001*total.n.ug)/sample.weight.mg)*100) %>% 
  mutate(C.N.ratio = total.c.ug/total.n.ug)

#standardize (center and scale by standard deviation) all continuous predictors;
#the scale function centers and scales by default, and the c() on each scale()
#function removes additional attributes created by the scale function that may
#cause issues with other packages
data$salmon.density.scaled <- c(scale(data$salmon.density))
data$chum.density.scaled <- c(scale(data$chum.density))
data$pink.density.scaled <- c(scale(data$pink.density))
data$dist.upstream.scaled <- c(scale(data$dist.upstream))
data$dist.from.stream.scaled <- c(scale(data$dist.from.stream))
data$slope.scaled <- c(scale(data$slope))
data$northness.scaled <- c(scale(data$northness))
data$eastness.scaled <- c(scale(data$eastness))
data$avg.canopy.cover.scaled <- c(scale(data$avg.canopy.cover))
data$rel.soil.moisture.scaled <- c(scale(data$rel.soil.moisture))

#remove 5 streams that are a part of another study
data <- data %>% dplyr::filter(stream %in% c("Beales","Bullock Main","Clatse","Fancy",
"Fannie", "Farm Bay","Goatbushu","Hooknose","Jane","Kill","Kunsoot","Lee",
"Quartcha","Sagar"))

####01. NITROGEN ISOTOPE MODELLING####
isotope_mod <- glmmTMB(n.15 ~ species * log(salmon.density.scaled + 1)
                         + dist.upstream.scaled +
                         dist.from.stream.scaled #+ northness.scaled + eastness.scaled +
                        + rel.soil.moisture.scaled + avg.canopy.cover.scaled 
                        + slope +
                   + (1|stream/quadrat.id),
                   family = gaussian(), #ziformula = ~0,
                   #dispformula = ~dist.upstream.scaled,
                   data = data, na.action = "na.omit")

pal <- pnw_palette("Cascades", 5)
ggplot(data, aes(fill = species, y = n.15,  x = salmon.density.scaled)) + 
  geom_point() +
  geom_smooth() +
  scale_fill_manual(values = pal) +
  theme_classic(30) +
  theme(axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", size = 1, fill = NA),
        #panel.spacing.x = unit(0.3,"line"),
        #panel.spacing.y = unit(0.3,"line"),
        strip.background = element_rect(colour = "white",size = 1, fill = NA),
        line = element_line(size = 0.25),
        strip.text.y.right = element_text(angle = 0)) +
  labs(x = "Salmon Density", y = "nitrogen-15")

#check variance inflation factors for collinearity among predictors
performance::check_collinearity(isotope_mod)

#check residuals and test model assumptions ising the DHARMa package
res <- simulateResiduals(isotope_mod, n = 1000)
plot(res)

#Check Individual Predictors
data2 <- tidyr::drop_na(data)
plotResiduals(res$scaledResiduals, data2$salmon.density.scaled)
plotResiduals(res$scaledResiduals, data2$species)
plotResiduals(res$scaledResiduals, data2$dist.upstream.scaled)
plotResiduals(res$scaledResiduals, data2$dist.from.stream.scaled)
plotResiduals(res$scaledResiduals, data2$northness.scaled)
plotResiduals(res$scaledResiduals, data2$eastness.scaled)
plotResiduals(res$scaledResiduals, data2$avg.canopy.cover.scaled)
plotResiduals(res$scaledResiduals, data2$rel.soil.moisture.scaled)


# Goodness of fit test
testUniformity(simulationOutput = res) 
# Deviation issue most likely reflective of large sample size (DHARMa vignette),
# dispersion looks good, line very straight given sample size

#Check for outliers
testOutliers(simulationOutput = res) ## outliers likely due to sample size
#will not remove any simply to improve fit

# Check Dispersion
testDispersion(res) # No over/underdispersion 
#(Small p-value (<0.05) indicates dispersion problem)

#once model assumptions and fit have been checked, assess the model outputs
summary(isotope_mod)

####02. %N and C:N MODELLING####

percent_n_mod <- glmmTMB(percent.N ~ salmon.density.scaled + 1 * species + dist.upstream.scaled +
                         dist.from.stream.scaled + northness.scaled + eastness.scaled +
                         avg.canopy.cover.scaled + rel.soil.moisture.scaled
                       + (1|stream/quadrat.id),
                       family = gaussian(), #ziformula = ~0,
                       #dispformula = ~1,
                       data = data, na.action = "na.omit")
res <- simulateResiduals(percent_n_mod)
plot(res)
summary(percent_n_mod)

c_n_mod <- glmmTMB(C.N.ratio ~ log(salmon.density) * species + dist.upstream.scaled +
                           dist.from.stream.scaled + northness.scaled + eastness.scaled +
                           avg.canopy.cover.scaled + rel.soil.moisture.scaled
                         + (1|stream/quadrat.id),
                         family = gaussian(), #ziformula = ~0,
                         #dispformula = ~1,
                         data = data, na.action = "na.omit")
summary(c_n_mod)
res <- simulateResiduals(c_n_mod)
plot(res)

####03. LEAF MASS PER AREA MODELLING####

mass_mod <- glmmTMB(punch.weight.mg ~ log(salmon.density) * species + dist.upstream.scaled +
                         dist.from.stream.scaled + northness.scaled + eastness.scaled +
                         avg.canopy.cover.scaled + rel.soil.moisture.scaled
                       + (1|stream/quadrat.id),
                       family = gaussian(), #ziformula = ~0,
                       #dispformula = ~1,
                       data = data, na.action = "na.omit")
summary(mass_mod)
res <- simulateResiduals(mass_mod)
plot(res)

####04. LEAF AREA MODELLING####

area_mod <- glmmTMB(leaf.area ~ log(salmon.density) * species + dist.upstream.scaled +
                         dist.from.stream.scaled + northness.scaled + eastness.scaled +
                         avg.canopy.cover.scaled + rel.soil.moisture.scaled
                       + (1|stream/quadrat.id),
                       family = gaussian(), #ziformula = ~0,
                       #dispformula = ~1,
                       data = data, na.action = "na.omit")
summary(area_mod)
res <- simulateResiduals(area_mod)
plot(res)

####05. LEAF GREENNESS MODELLING####

greenness_mod <- glmmTMB(percent.green ~ log(salmon.density) * species + dist.upstream.scaled +
                         dist.from.stream.scaled + northness.scaled + eastness.scaled +
                         avg.canopy.cover.scaled + rel.soil.moisture.scaled
                       + (1|stream/quadrat.id),
                       family = gaussian(), #ziformula = ~0,
                       #dispformula = ~1,
                       data = data, na.action = "na.omit")
summary(greenness_mod)
res <- simulateResiduals(greenness_mod)
plot(res)
