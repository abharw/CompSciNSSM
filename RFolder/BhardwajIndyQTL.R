# R script for analyzing Shimomura  QTL data
# Arav Bhardwaj
# BhardwajIndyQTL.R script
# October 15, 2022

# Clean things up
rm(list=ls())
setwd("/Users/abhardwaj24/Desktop/RFolder")

# Install and load QTL package
install.packages("qtl")
library(qtl)

# Load in data
cross <- read.cross("csv", file = "shimomura.csv", genotypes = c("A", "H", "B"),
                    na.strings = "-",alleles = c("A", "B"))

# Get summary statistics
summary(cross)

# Jittermap
jittermap(cross)

# HISTOGRAMS #
# Phase angle of entrainment
phase <- cross$pheno$phase
hist(phase, xlab = "Phase", main = "Phase Angle of Entrainment")
# Amplitude of circadian rhythmicity
amp <- cross$pheno$amp
hist(amp, xlab = "Amplitude", main = "Amplitude of Circadian Rhythmicity")
# Circadian activity level
act <- cross$pheno$act
hist(act, xlab = "Activity Level", main = "Circadian Activity Level")

# Genetic map and missing data
plot.map(cross)
plotMissing(cross)

# QQ PLOTS #
# Phase angle of entrainment
qqnorm(phase, main = "Phase Angle of Entrainment")
qqline(phase)

# Amplitude of circadian rhythmicity
qqnorm(amp, main = "Amplitude of Circadian Rhythmicity")
qqline(amp)

# Circadian activity level
qqnorm(act, main = "Circadian Activity Level")
qqline(act)

# CROSS CALCS #
cross <- calc.genoprob(cross, step = 2.0, off.end = 0, error.prob = 1.0e-4,
                       map.function = "haldane", stepwidth = "fixed")

cross <- sim.geno(cross, step = 2.0, off.end = 0, error.prob = 1.0e-4,
                  map.function = "haldane", stepwidth = "fixed")


# MAIN SCANS #
# Phase angle of entrainment
cross.scanPhase <- scanone(cross, pheno.col = 7, model = "normal", method = "em")
cross.scanPhase.perm <- scanone(cross, pheno.col = 7, model = "normal", method = "em", n.perm = 100)
# Amplitude of circadian rhythmicity
cross.scanAmp <- scanone(cross, pheno.col = 8, model = "normal", method = "em")
cross.scanAmp.perm <- scanone(cross, pheno.col = 8, model = "normal", method = "em", n.perm = 100)
# Circadian activity level
cross.scanAct <- scanone(cross, pheno.col = 9, model = "normal", method = "em")
cross.scanAct.perm <- scanone(cross, pheno.col = 9, model = "normal", method = "em", n.perm = 100)

summary(cross.scanPhase)
summary(cross.scanAmp)
summary(cross.scanAct)

# PLOT MAIN SCANS #
# Phase angle of entrainment
plot(cross.scanPhase, main="Phase Angle of Entrainment")
threshPh <- summary(cross.scanPhase.perm, alpha = c(0.37, 0.10, 0.05))
abline(h = threshPh[1], col = "blue")
abline(h = threshPh[2], col = "red")
abline(h = threshPh[3], col = "green")
# Amplitude of circadian rhythmicity
plot(cross.scanAmp, main="Amplitude of Circadian Rhythmicity")
threshAmp <- summary(cross.scanAmp.perm, alpha = c(0.37, 0.10, 0.05))
abline(h = threshAmp[1], col = "blue")
abline(h = threshAmp[2], col = "red")
abline(h = threshAmp[3], col = "green")
# Circadian activity level
plot(cross.scanAct, main="Circadian Activity Level")
threshHW <- summary(cross.scanAct.perm, alpha = c(0.37, 0.10, 0.05))
abline(h = threshHW[1], col = "blue")
abline(h = threshHW[2], col = "red")
abline(h = threshHW[3], col = "green")

# SUMMARY #
# Phase
summary(cross.scanPhase, perm = cross.scanPhase.perm, alpha = 0.05)
# Amp
summary(cross.scanAmp, perm = cross.scanAmp.perm, alpha = 0.05)
# Act
summary(cross.scanAct, perm = cross.scanAct.perm, alpha = 0.05)

# EFFECT PLOTS #
# Blood pressure
firstPhase <- find.marker(cross, chr = 7, pos = 34.7)
effectplot(cross, pheno.col = 7, mname1 = firstPhase)
secondPhase <- find.marker(cross, chr = 7, pos = 54.8)
effectplot(cross, pheno.col = 7, mname1 = secondPhase)
