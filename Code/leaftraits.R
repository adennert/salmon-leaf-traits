####ALLISON DENNERT, Cross-watershed salmon lead traits analyses
####"MANUSCRIPT TITLE" by A. Dennert, E. Elle, J. Reynolds

#This R code analayzes and visualizes that relationships between spawning 
#salmon and:
#1) nitrogen isotopes, 
#2) % nitrogen 
#3) leaf mass per area,
#4) leaf area, and 
#5) leaf greenness 
#in salmonberry (RUSP), false lily-of-the-valley (MADI), false azalea (MEFE), 
#blueberry species (VASP), and foamflower (TITR).
#These data were collected across 14 streams of varying chum and pink salmon 
#spawning densities, along with a variety of environmental co-variates.

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
library(ggeffects)

####00. READ, CLEAN, TRANSFORM DATA####
set.seed(1)
data <- read.csv("Raw Data/final_data.csv", strip.white = TRUE)
str(data)
data$species <- as.factor(data$species)
data$stream <- as.factor(data$stream)
data$quadrat.id <- as.factor(data$quadrat.id)

levels(data$species)
data$species = factor(data$species,levels(data$species)[c(1,3,2,5,4)])

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

#Hocking et al. 2011 found that a non-linear relationship between salmon density
#and predictors had the best model fit; thus, salmon density will be 
#transformed using log(x + 1); this transformation must be done before scaling,
#as scaling produces negative numbers which log transforming cannot handle
data <- data %>% mutate(log.salmon.density = log(salmon.density + 1))

#standardize (center and scale by standard deviation) all continuous predictors;
#the scale function centers and scales by default, and the c() on each scale()
#function removes additional attributes created by the scale function that may
#cause issues with other packages
data$salmon.density.scaled <- c(scale(data$salmon.density))
data$log.salmon.density.scaled <- c(scale(data$log.salmon.density))
data$chum.density.scaled <- c(scale(data$chum.density))
data$pink.density.scaled <- c(scale(data$pink.density))
data$dist.upstream.scaled <- c(scale(data$dist.upstream))
data$dist.from.stream.scaled <- c(scale(data$dist.from.stream))
data$slope.scaled <- c(scale(data$slope))
data$northness.scaled <- c(scale(data$northness))
data$eastness.scaled <- c(scale(data$eastness))
data$avg.canopy.cover.scaled <- c(scale(data$avg.canopy.cover))
data$rel.soil.moisture.scaled <- c(scale(data$rel.soil.moisture))
data$percent.N.scaled <- c(scale(data$percent.N))

#remove 5 streams that are a part of another study; keep 14 streams
data <- data %>% dplyr::filter(stream %in% c("Beales", "Bullock Main", "Clatse",
                                             "Fancy", "Fannie", "Farm Bay",
                                             "Goatbushu", "Hooknose", "Jane",
                                             "Kill", "Kunsoot", "Lee", "Quartcha",
                                             "Sagar")) %>% 
#remove 5th species that is not of interest and has limited data
  filter(species != "TITR")

####01. NITROGEN ISOTOPE MODELLING####
hist(data$n.15, breaks = 100) #use normal distribution in lme4

isotope_mod <- lme4::lmer(n.15 ~ species * log.salmon.density.scaled +
                         dist.upstream.scaled +
                         dist.from.stream.scaled + 
                         rel.soil.moisture.scaled + 
                         avg.canopy.cover.scaled +
                         slope.scaled
                       + (1|stream/quadrat.id),
                       data = data, na.action = "na.omit")

#check variance inflation factors for collinearity among predictors
performance::check_collinearity(isotope_mod)

#check residuals and test model assumptions using the DHARMa package
res <- simulateResiduals(isotope_mod, n = 1000)
plot(res)

#Check Individual Predictors
data_na <- data %>% tidyr::drop_na(n.15) 
data_na$species <- as.character(data_na$species) #changed to enable visualization
plotResiduals(res$scaledResiduals, data_na$log.salmon.density.scaled)
plotResiduals(res$scaledResiduals, data_na$species)
plotResiduals(res$scaledResiduals, data_na$dist.upstream.scaled)
plotResiduals(res$scaledResiduals, data_na$dist.from.stream.scaled)
plotResiduals(res$scaledResiduals, data_na$rel.soil.moisture.scaled)
plotResiduals(res$scaledResiduals, data_na$avg.canopy.cover.scaled)
plotResiduals(res$scaledResiduals, data_na$slope.scaled)
#potential for heteroscedasiticy

#Check for outliers
testOutliers(simulationOutput = res) #limited number of outliers; 
#likely due to sample size (F. Hartig, DHARMa vignette); will not remove any 
#simply to improve fit; outlier test not significant

#Check Dispersion
testDispersion(res) # No over/underdispersion 

#once model assumptions and fit have been checked, assess the model outputs
summary(isotope_mod)

#check for evidence of heteroscedasticity
performance::check_heteroscedasticity(isotope_mod)
#heteroscedasticity present

#calculate heteroscedasticity-robust sandwich estimators of standard errors
sand <- clubSandwich::vcovCR(isotope_mod, type = "CR2")
clubSandwich::coef_test(isotope_mod, vcov = sand, test = "z")
#From clubSandwich vignette: "CR2" is the "bias-reduced linearization" 
#adjustment proposed by Bell and McCaffrey (2002) and further developed in 
#Pustejovsky and Tipton (2017). The adjustment is chosen so that the
#variance-covariance estimator is exactly unbiased under a user-specified 
#working model

#calculate the conditional R^2, which accounts for both fixed and random effects
performance::r2(isotope_mod)
#Conditional R^2 = 0.714

#calculate species-level relationships between salmon and response
emmeans::emtrends(isotope_mod, ~species, 
                  var = "log.salmon.density.scaled")
emmeans::emmip(isotope_mod, species ~ log.salmon.density.scaled, cov.reduce = range)

####02. %N MODELLING####
hist(data$percent.N, breaks = 100) #use normal distribution in lme

percent_mod <- lme4::lmer(percent.N ~ species * log.salmon.density.scaled +
                            dist.upstream.scaled +
                            dist.from.stream.scaled + 
                            rel.soil.moisture.scaled + 
                            avg.canopy.cover.scaled +
                            slope.scaled
                          + (1|stream/quadrat.id),
                          data = data, na.action = "na.omit")

#check variance inflation factors for collinearity among predictors
performance::check_collinearity(percent_mod)

#check residuals and test model assumptions using the DHARMa package
res <- simulateResiduals(percent_mod, n = 1000)
plot(res)

#Check Individual Predictors
data_na <- data %>% tidyr::drop_na(percent.N) 
data_na$species <- as.character(data_na$species)
plotResiduals(res$scaledResiduals, data_na$log.salmon.density.scaled)
plotResiduals(res$scaledResiduals, data_na$species)
plotResiduals(res$scaledResiduals, data_na$dist.upstream.scaled)
plotResiduals(res$scaledResiduals, data_na$dist.from.stream.scaled)
plotResiduals(res$scaledResiduals, data_na$rel.soil.moisture.scaled)
plotResiduals(res$scaledResiduals, data_na$avg.canopy.cover.scaled)
plotResiduals(res$scaledResiduals, data_na$slope.scaled)
#potential for heteroscedasiticy

#Check for outliers
testOutliers(simulationOutput = res) #limited number of outliers; 
#likely due to sample size (F. Hartig, DHARMa vignette); will not remove any 
#simply to improve fit; outlier test not significant

#Check Dispersion
testDispersion(res) # No over/underdispersion 

#once model assumptions and fit have been checked, assess the model outputs
summary(percent_mod)

#check for evidence of heteroscedasticity
performance::check_heteroscedasticity(percent_mod)
#heteroscedasticity not present

#will calculate conservative standard errors despite homoscedasticity;
#calculate heteroscedasticity-robust sanwich estimators of standard errors
sand <- clubSandwich::vcovCR(percent_mod, type = "CR2")
clubSandwich::coef_test(percent_mod, vcov = sand, test = "z")
#From clubSandwich vignette: "CR2" is the "bias-reduced linearization" 
#adjustment proposed by Bell and McCaffrey (2002) and further developed in 
#Pustejovsky and Tipton (2017). The adjustment is chosen so that the
#variance-covariance estimator is exactly unbiased under a user-specified 
#working model

#calculate the conditional R^2, which accounts for both fixed and random effects
performance::r2(percent_mod)
#Conditional R^2 = 0.478

#calculate species-level relationships between salmon and response
emmeans::emtrends(percent_mod, ~species, 
                  var = "log.salmon.density.scaled")
emmeans::emmip(percent_mod, species ~ log.salmon.density.scaled, cov.reduce = range)

####03. LEAF MASS PER AREA MODELLING####
hist(data$punch.weight.mg, breaks = 100) #use gamma distribution in glmmTMB

#must use glmmTMB has lme4 models won't converge
mass_mod <- glmmTMB(punch.weight.mg ~ species * log.salmon.density.scaled +
                         dist.upstream.scaled +
                         dist.from.stream.scaled + 
                         rel.soil.moisture.scaled + 
                         avg.canopy.cover.scaled +
                         slope.scaled
                       + (1|stream/quadrat.id),
                       dispformula = ~ log.salmon.density.scaled +
                         species + avg.canopy.cover.scaled,
                       family = Gamma(link = "log"),
                       data = data, na.action = "na.omit")

#check variance inflation factors for collinearity among predictors
performance::check_collinearity(mass_mod)

#check residuals and test model assumptions using the DHARMa package
res <- simulateResiduals(mass_mod, n = 1000)
plot(res)

#Check Individual Predictors
data_na <- data %>% tidyr::drop_na(punch.weight.mg) 
data_na$species <- as.character(data_na$species)
plotResiduals(res$scaledResiduals, data_na$log.salmon.density.scaled)
plotResiduals(res$scaledResiduals, data_na$species)
plotResiduals(res$scaledResiduals, data_na$dist.upstream.scaled)
plotResiduals(res$scaledResiduals, data_na$dist.from.stream.scaled)
plotResiduals(res$scaledResiduals, data_na$rel.soil.moisture.scaled)
plotResiduals(res$scaledResiduals, data_na$avg.canopy.cover.scaled)
plotResiduals(res$scaledResiduals, data_na$slope.scaled)
#potential for heteroscedasticity at only one predictor

#Check for outliers
testOutliers(simulationOutput = res) #limited number of outliers; 
#likely due to sample size (F. Hartig, DHARMa vignette); will not remove any 
#simply to improve fit; outlier test not significant

#Check Dispersion
testDispersion(res) # No over/underdispersion 

#once model assumptions and fit have been checked, assess the model outputs 
summary(mass_mod)

#check for evidence of heteroscedasticity
performance::check_heteroscedasticity(mass_mod)
#cannot check with this function on non-gaussian models

#cannot calculate heteroscedasticity-robust sandwich estimators of standard errors
#with glmmTMB; thus, interpretation of coefficients with confidence intervals
#near zero must be conservative

#calculate the conditional R^2, which accounts for both fixed and random effects
performance::r2(mass_mod)
#after removing model's dispersion parameters, conditional R^2 = 0.387

#calculate species-level relationships between salmon and response
emmeans::emtrends(mass_mod, ~species, 
                  var = "log.salmon.density.scaled")
emmeans::emmip(mass_mod, species ~ log.salmon.density.scaled, cov.reduce = range)

####04. LEAF AREA MODELLING####
hist(data$leaf.area, breaks = 100) #use gamma distribution in glmmTMB

area_mod <- glmmTMB(leaf.area ~ species * log.salmon.density.scaled +
                      dist.upstream.scaled +
                      dist.from.stream.scaled + 
                      #northness.scaled + 
                      #eastness.scaled +
                      rel.soil.moisture.scaled + 
                      avg.canopy.cover.scaled +
                      slope.scaled
                    + (1|stream/quadrat.id),
                    dispformula = ~ log.salmon.density.scaled +
                    species,
                    family = Gamma(link = "log"),
                    data = data, na.action = "na.omit")

#check variance inflation factors for collinearity among predictors
performance::check_collinearity(area_mod)

#check residuals and test model assumptions using the DHARMa package
res <- simulateResiduals(area_mod, n = 1000)
plot(res)

#Check Individual Predictors
data_na <- data %>% tidyr::drop_na(leaf.area) 
data_na$species <- as.character(data_na$species)
plotResiduals(res$scaledResiduals, data_na$log.salmon.density.scaled)
plotResiduals(res$scaledResiduals, data_na$species)
plotResiduals(res$scaledResiduals, data_na$dist.upstream.scaled)
plotResiduals(res$scaledResiduals, data_na$dist.from.stream.scaled)
plotResiduals(res$scaledResiduals, data_na$rel.soil.moisture.scaled)
plotResiduals(res$scaledResiduals, data_na$avg.canopy.cover.scaled)
plotResiduals(res$scaledResiduals, data_na$slope.scaled)
#potential for heteroscedasticity at only one predictor

#Check for outliers
testOutliers(simulationOutput = res) #no outliers

#Check Dispersion
testDispersion(res) # No over/underdispersion 

#once model assumptions and fit have been checked, assess the model outputs
summary(area_mod)

#check for evidence of heteroscedasticity
performance::check_heteroscedasticity(area_mod)
#cannot check with this function on non-gaussian models

#cannot calculate heteroscedasticity-robust sandwich estimators of standard errors
#with glmmTMB; thus, interpretation of coefficients with confidence intervals
#near zero must be conservative

#calculate the conditional R^2, which accounts for both fixed and random effects
performance::r2(area_mod)
#after removing model's dispersion parameters, conditional R^2 = 0.695

#calculate species-level relationships between salmon and response
emmeans::emtrends(area_mod, ~species, 
                  var = "log.salmon.density.scaled")
emmeans::emmip(area_mod, species ~ log.salmon.density.scaled, cov.reduce = range)

####05. LEAF GREENNESS MODELLING####
hist(data$percent.green, breaks = 100) #use normal distribution in lme4

green_mod <- lme4::lmer(percent.green ~ species * log.salmon.density.scaled +
                       dist.upstream.scaled +
                       dist.from.stream.scaled + 
                       rel.soil.moisture.scaled + 
                       avg.canopy.cover.scaled +
                       slope.scaled
                     + (1|stream/quadrat.id),
                     data = data, na.action = "na.omit")

#check variance inflation factors for collinearity among predictors
performance::check_collinearity(green_mod)

#check residuals and test model assumptions using the DHARMa package
res <- simulateResiduals(green_mod, n = 1000)
plot(res)

#Check Individual Predictors
data_na <- data %>% tidyr::drop_na(percent.green) 
data_na$species <- as.character(data_na$species)
plotResiduals(res$scaledResiduals, data_na$log.salmon.density.scaled)
plotResiduals(res$scaledResiduals, data_na$species)
plotResiduals(res$scaledResiduals, data_na$dist.upstream.scaled)
plotResiduals(res$scaledResiduals, data_na$dist.from.stream.scaled)
plotResiduals(res$scaledResiduals, data_na$rel.soil.moisture.scaled)
plotResiduals(res$scaledResiduals, data_na$avg.canopy.cover.scaled)
plotResiduals(res$scaledResiduals, data_na$slope.scaled)
#potential for heteroscedasticity

#Check for outliers
testOutliers(simulationOutput = res) #no outliers present

#Check Dispersion
testDispersion(res) # No over/underdispersion 

#once model assumptions and fit have been checked, assess the model outputs
summary(green_mod)

#check for evidence of heteroscedasticity
performance::check_heteroscedasticity(green_mod)

#heteroscedasticity present
#calculate heteroscedasticity-robust sandwich estimators of standard errors
sand <- clubSandwich::vcovCR(green_mod, type = "CR2")
clubSandwich::coef_test(green_mod, vcov = sand, test = "z")
#From clubSandwich vignette: "CR2" is the "bias-reduced linearization" 
#adjustment proposed by Bell and McCaffrey (2002) and further developed in 
#Pustejovsky and Tipton (2017). The adjustment is chosen so that the
#variance-covariance estimator is exactly unbiased under a user-specified 
#working model

#calculate the conditional R^2, which account for both fixed and random effects
performance::r2(green_mod)
#conditional R^2 = 0.595

#calculate species-level relationships between salmon and response
emmeans::emtrends(green_mod, ~species, 
                  var = "log.salmon.density.scaled")
emmeans::emmip(green_mod, species ~ log.salmon.density.scaled, cov.reduce = range)

####06. FIGURES####

#Nitrogen-15 Figure----

#create model predictions using heteroscedasticity-consistent standard errors
predict_n.15 <- ggpredict(isotope_mod, vcov.fun = "vcovHC", vcov.type = "CR2",
                          terms = c("log.salmon.density.scaled[n=100]",
                                              "species")) %>% 
  #rename columns to column names in original data
  rename(log.salmon.density.scaled = x,
         species = group) %>% 
  mutate(log.salmon.density = 
           log.salmon.density.scaled*sd(data$log.salmon.density)+
           mean(data$log.salmon.density)) %>% 
  #limit predictions to the observed values rather than extrapolating
  filter((species == 'RUSP' & 
            (log.salmon.density.scaled >= 
               min(filter(data, species == 'RUSP')$log.salmon.density.scaled) & 
               log.salmon.density.scaled <= 
            max(filter(data, species == 'RUSP')$log.salmon.density.scaled))) |
         (species == 'MADI' & 
            (log.salmon.density.scaled >= 
               min(filter(data, species == 'MADI')$log.salmon.density.scaled) & 
               log.salmon.density.scaled <= 
            max(filter(data, species == 'MADI')$log.salmon.density.scaled))) |
         (species == 'VASP' & 
            (log.salmon.density.scaled >= 
               min(filter(data, species == 'VASP')$log.salmon.density.scaled) & 
               log.salmon.density.scaled <= 
            max(filter(data, species == 'VASP')$log.salmon.density.scaled))) |
         (species == 'MEFE' & 
            (log.salmon.density.scaled >= 
                min(filter(data, species == 'MEFE')$log.salmon.density.scaled) & 
                log.salmon.density.scaled <= 
            max(filter(data, species == 'MEFE')$log.salmon.density.scaled))))

#create a custom colour palette
  pal <- pnw_palette("Cascades", 5)
  
#plot the raw data with model predictions on top
  
ggplot(data, aes(x = log.salmon.density, y = n.15)) +
geom_point(aes(colour = species)) +
geom_ribbon(data = predict_n.15,
            aes(y = predicted, ymin = conf.low, ymax = conf.high,
                fill = species), 
            alpha = 0.2) +
geom_line(data = predict_n.15,
          aes(y = predicted, colour = species)) +
scale_color_manual(values = pal, labels = c("False Lily","Salmonberry",
                                            "False Azalea", "Blueberry"),
                   name = "") +
scale_fill_manual(values = pal, labels = c("False Lily","Salmonberry",
                                           "False Azalea", "Blueberry"),
                  name = "") +
labs(x = "log(Salmon Density)", y = expression(paste(delta^{15} ~ N))) +
theme_classic(30) +
facet_wrap(~species, scales = "fixed") +
theme(strip.text = element_blank(), axis.line = element_line(size = 0.25),
      legend.text = element_text(size=30)) +
annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 2)+
annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 2)

ggsave("Figures/isotopes.png",  height=9, width=16)
ggsave("Figures/isotopes.pdf",  height=9, width=16)

#Percent Nitrogen Figure----
  
#create model predictions using heteroscedasticity-consistent standard errors
predict_percent <- ggpredict(percent_mod, vcov.fun = "vcovHC", vcov.type = "CR2",
                             terms = c("log.salmon.density.scaled[n=100]",
                                                   "species")) %>% 
  #rename columns to column names in original data
  rename(log.salmon.density.scaled = x,
         species = group) %>% 
  mutate(log.salmon.density = 
           log.salmon.density.scaled*sd(data$log.salmon.density)+
           mean(data$log.salmon.density)) %>% 
  #limit predictions to the observed values rather than extrapolating
  filter((species == 'RUSP' & 
            (log.salmon.density.scaled >= 
               min(filter(data, species == 'RUSP')$log.salmon.density.scaled) & 
               log.salmon.density.scaled <= 
               max(filter(data, species == 'RUSP')$log.salmon.density.scaled))) |
           (species == 'MADI' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'MADI')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'MADI')$log.salmon.density.scaled))) |
           (species == 'VASP' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'VASP')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'VASP')$log.salmon.density.scaled))) |
           (species == 'MEFE' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'MEFE')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'MEFE')$log.salmon.density.scaled))))

#create a custom colour palette
pal <- pnw_palette("Cascades", 5)

ggplot(data, aes(x = log.salmon.density, y = percent.N)) +
  geom_point(aes(colour = species)) +
  geom_ribbon(data = predict_percent,
              aes(y = predicted, ymin = conf.low, ymax = conf.high,
                  fill = species), 
              alpha = 0.2) +
  geom_line(data = predict_percent,
            aes(y = predicted, colour = species)) +
  scale_color_manual(values = pal, labels = c("False Lily","Salmonberry",
                                              "False Azalea", "Blueberry"),
                     name = "") +
  scale_fill_manual(values = pal, labels = c("False Lily","Salmonberry",
                                             "False Azalea", "Blueberry"),
                    name = "") +
  labs(x = "log(Salmon Density)", y = "% Nitrogen") +
  theme_classic(30) +
  facet_wrap(~species, scales = "fixed") +
  theme(strip.text = element_blank(), axis.line = element_line(size = 0.25),
        legend.text = element_text(size=30)) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 2)+
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 2)

ggsave("Figures/percent_nitrogen.png",  height=9, width=16)
ggsave("Figures/percent_nitrogen.pdf",  height=9, width=16)

#Leaf Mass Figure----

#create model predictions using heteroscedasticity consistent standard errors
predict_mass <- ggpredict(mass_mod, vcov.fun = "vcovHC", vcov.type = "CR2",
                             terms = c("log.salmon.density.scaled[n=100]",
                                       "species")) %>% 
  #rename columns to column names in original data
  rename(log.salmon.density.scaled = x,
         species = group) %>% 
  mutate(log.salmon.density = 
           log.salmon.density.scaled*sd(data$log.salmon.density)+
           mean(data$log.salmon.density)) %>% 
  #limit predictions to the observed values rather than extrapolating
  filter((species == 'RUSP' & 
            (log.salmon.density.scaled >= 
               min(filter(data, species == 'RUSP')$log.salmon.density.scaled) & 
               log.salmon.density.scaled <= 
               max(filter(data, species == 'RUSP')$log.salmon.density.scaled))) |
           (species == 'MADI' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'MADI')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'MADI')$log.salmon.density.scaled))) |
           (species == 'VASP' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'VASP')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'VASP')$log.salmon.density.scaled))) |
           (species == 'MEFE' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'MEFE')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'MEFE')$log.salmon.density.scaled))))

#create a custom colour palette
pal <- pnw_palette("Cascades", 5)

#plot the raw data with model predictions on top
ggplot(data, aes(x = log.salmon.density, y = punch.weight.mg)) +
  geom_point(aes(colour = species)) +
  geom_ribbon(data = predict_mass,
              aes(y = predicted, ymin = conf.low, ymax = conf.high,
                  fill = species), 
              alpha = 0.2) +
  geom_line(data = predict_mass,
            aes(y = predicted, colour = species)) +
  scale_color_manual(values = pal, labels = c("False Lily","Salmonberry",
                                              "False Azalea", "Blueberry"),
                     name = "") +
  scale_fill_manual(values = pal, labels = c("False Lily","Salmonberry",
                                             "False Azalea", "Blueberry"),
                    name = "") +
  labs(x = "log(Salmon Density)", y = "Leaf Mass") +
  theme_classic(30) +
  facet_wrap(~species, scales = "fixed") +
  theme(strip.text = element_blank(), axis.line = element_line(size = 0.25),
        legend.text = element_text(size=30)) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 2)+
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 2)

ggsave("Figures/leaf_mass.png",  height=9, width=16)
ggsave("Figures/leaf_mass.pdf",  height=9, width=16)

#Leaf Area Figure----

#create model predictions using heteroscedasticity consistent standard errors
predict_area <- ggpredict(area_mod, vcov.fun = "vcovHC", vcov.type = "CR2",
                             terms = c("log.salmon.density.scaled[n=100]",
                                       "species")) %>% 
  #rename columns to column names in original data
  rename(log.salmon.density.scaled = x,
         species = group) %>% 
  mutate(log.salmon.density = 
           log.salmon.density.scaled*sd(data$log.salmon.density)+
           mean(data$log.salmon.density)) %>% 
  #limit predictions to the observed values rather than extrapolating
  filter((species == 'RUSP' & 
            (log.salmon.density.scaled >= 
               min(filter(data, species == 'RUSP')$log.salmon.density.scaled) & 
               log.salmon.density.scaled <= 
               max(filter(data, species == 'RUSP')$log.salmon.density.scaled))) |
           (species == 'MADI' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'MADI')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'MADI')$log.salmon.density.scaled))) |
           (species == 'VASP' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'VASP')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'VASP')$log.salmon.density.scaled))) |
           (species == 'MEFE' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'MEFE')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'MEFE')$log.salmon.density.scaled))))

#create a custom colour palette
pal <- pnw_palette("Cascades", 5)

#plot the raw data with model predictions on top
ggplot(data, aes(x = log.salmon.density, y = leaf.area)) +
  geom_point(aes(colour = species)) +
  geom_ribbon(data = predict_area,
              aes(y = predicted, ymin = conf.low, ymax = conf.high,
                  fill = species), 
              alpha = 0.2) +
  geom_line(data = predict_area,
            aes(y = predicted, colour = species)) +
  scale_color_manual(values = pal, labels = c("False Lily","Salmonberry",
                                              "False Azalea", "Blueberry"),
                     name = "") +
  scale_fill_manual(values = pal, labels = c("False Lily","Salmonberry",
                                             "False Azalea", "Blueberry"),
                    name = "") +
  labs(x = "log(Salmon Density)", y = "Leaf Area") +
  theme_classic(30) +
  facet_wrap(~species, scales = "free_y") +
  theme(strip.text = element_blank(), axis.line = element_line(size = 1),
        legend.text = element_text(size=30)) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 2)+
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 2)

ggsave("Figures/leaf_area.png",  height=9, width=16)
ggsave("Figures/leaf_area.pdf",  height=9, width=16)

#Leaf Greenness Figure----

#create model predictions using heteroscedasticity consistent standard errors
predict_green <- ggpredict(green_mod, vcov.fun = "vcovHC", vcov.type = "CR2",
                             terms = c("log.salmon.density.scaled[n=100]",
                                       "species")) %>% 
  #rename columns to column names in original data
  rename(log.salmon.density.scaled = x,
         species = group) %>% 
  mutate(log.salmon.density = 
           log.salmon.density.scaled*sd(data$log.salmon.density)+
           mean(data$log.salmon.density)) %>% 
  #limit predictions to the observed values rather than extrapolating
  filter((species == 'RUSP' & 
            (log.salmon.density.scaled >= 
               min(filter(data, species == 'RUSP')$log.salmon.density.scaled) & 
               log.salmon.density.scaled <= 
               max(filter(data, species == 'RUSP')$log.salmon.density.scaled))) |
           (species == 'MADI' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'MADI')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'MADI')$log.salmon.density.scaled))) |
           (species == 'VASP' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'VASP')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'VASP')$log.salmon.density.scaled))) |
           (species == 'MEFE' & 
              (log.salmon.density.scaled >= 
                 min(filter(data, species == 'MEFE')$log.salmon.density.scaled) & 
                 log.salmon.density.scaled <= 
                 max(filter(data, species == 'MEFE')$log.salmon.density.scaled))))

#create a custom colour palette
pal <- pnw_palette("Cascades", 5)

#plot the raw data with model predictions on top
ggplot(data, aes(x = log.salmon.density, y = percent.green)) +
  geom_point(aes(colour = species)) +
  geom_ribbon(data = predict_green,
              aes(y = predicted, ymin = conf.low, ymax = conf.high,
                  fill = species), 
              alpha = 0.2) +
  geom_line(data = predict_green,
            aes(y = predicted, colour = species)) +
  scale_color_manual(values = pal, labels = c("False Lily","Salmonberry",
                                              "False Azalea", "Blueberry"),
                     name = "") +
  scale_fill_manual(values = pal, labels = c("False Lily","Salmonberry",
                                             "False Azalea", "Blueberry"),
                    name = "") +
  labs(x = "log(Salmon Density)", y = "Proportion Green") +
  theme_classic(30) +
  facet_wrap(~species, scales = "fixed") +
  theme(strip.text = element_blank(), axis.line = element_line(size = 0.25),
        legend.text = element_text(size=30)) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 2)+
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 2)

ggsave("Figures/prop_green.png",  height=9, width=16)
ggsave("Figures/prop_green.pdf",  height=9, width=16)
