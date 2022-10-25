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

####00. READ, CLEAN, TRANSFORM DATA####
data <- read.csv("Raw Data/final_data.csv", strip.white = TRUE)

####01. NITROGEN ISOTOPE MODELLING####

####02. %N and C:N MODELLING####

####03. LEAF MASS PER AREA MODELLING####

####04. LEAF AREA MODELLING####

####05. LEAF GREENNESS MODELLING####