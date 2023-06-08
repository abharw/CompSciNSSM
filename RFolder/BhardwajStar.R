# Author: Arav Bhardwaj
# Date: 09-08-22
# StarPlot.R
# Analyzes star data

# ENVIRONMENT SETUP
rm(list=ls())
setwd("C:\Users\abhardwaj24\Desktop\RFolder")

# CODE

# Read in CSV
stars <- read.csv("http://chemistry.ncssm.edu/data/HIP_star.csv")

# Find names in dataset and save to a .txt file
sink("starnames.txt")
names(stars)
sink()

# Plot RA and DE
jpeg("starchart.jpg")
plot(stars$DE ~ stars$RA, main = "Plot of all stars with right ascension (RA) and declination (DE)",
    xlab = "RA", ylab = "DE", pch = "*")
dev.off()

# Plot bar plots in 2 x 2 matrix and save to .jpg file
jpeg("boxplots.jpg")
par(mfrow = c(2, 2))
boxplot(stars$Vmag, main = "Vmag")
boxplot(stars$Plx, main = "Plx")
boxplot(stars$pmRA, main = "pmRA")
boxplot(stars$pmDE, main = "pmDE")
# Reset plotting matrix
par(mfrow = c(1, 1))
dev.off()

# Filter stars with Vmag between 8 and 9 and save plot to .jpg file
jpeg("VmagsvsBV.jpg")
stars.filter <- (stars$Vmag >= 8 & stars$Vmag <= 9)
plot(stars$B.V[stars.filter] ~ stars$Vmag[stars.filter], main = "Plot of Filtered values and B.V",
     xlab = "Vmag[vmagfiltered]", ylab = "B.V[vmagfiltered]", pch="*")
dev.off()

# Calculate log (base 10) of Plx
logPlx <- log10(stars$Plx)
logL <- (15-stars$Vmag-5*logPlx)/(2.5)
# Save plot to .jpg
jpeg("BVvslogL.jpg")
# Plot logL vs. B.V. (filtered)
plot(logL[stars.filter] ~ stars$B.V[stars.filter], main = "Plot of B.V. and logL (luminosity)",
     xlab = "B.V[vmagfiltered]", ylab = "logL[vmagFiltered]", pch = "*")
# Regress the data
lmfit <- lm(logL[stars.filter] ~ stars$B.V[stars.filter])
# Plot line on top of data
abline(lmfit)
dev.off()

# EOF