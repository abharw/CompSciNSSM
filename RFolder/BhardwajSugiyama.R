# R script for analyzing Sugiyama QTL data
# Arav Bhardwaj
# BhardwajSugiyama.R script
# October 15, 2022

# Clean things up
rm(list=ls())
setwd("/Users/abhardwaj24/Desktop/RFolder")

# Install and load QTL package
install.packages("qtl")
library(qtl)

# Load in data
cross <- read.cross("csv", file = "sugiyamashort.csv", genotypes = c("C", "H", "B"),
                    na.strings = "-",alleles = c("C", "B"))

# Get summary statistics
summary(cross)

# HISTOGRAMS #
# Blood pressure
bp <- cross$pheno$BP_final
hist(bp, xlab = "BP", main = "Blood Pressure")
# Heart rate
hr <- cross$pheno$HR_final
hist(hr, xlab = "Heart Rate", main = "Heart Rate")
# Heart Weight
hw <- cross$pheno$heart_wt
hist(hw, xlab = "Heart Weight", main = "Heart Weight")

# Genetic map and missing data
plot.map(cross)
plotMissing(cross)

# QQ PLOTS #
# Blood pressure
qqnorm(bp, main = "Blood Pressure")
qqline(bp)

# Heart rate
qqnorm(hr, main = "Heart Rate")
qqline(hr)

# Heart weight
qqnorm(hw, main = "Heart Weight")
qqline(hw)

# CROSS CALCS #
cross <- calc.genoprob(cross, step = 2.0, off.end = 0, error.prob = 1.0e-4,
                       map.function = "haldane", stepwidth = "fixed")

cross <- sim.geno(cross, step = 2.0, off.end = 0, error.prob = 1.0e-4,
                  map.function = "haldane", stepwidth = "fixed")

# MAIN SCANS #
# Blood pressure
cross.scanBP <- scanone(cross, pheno.col = 3, model = "normal", method = "em")
cross.scanBP.perm <- scanone(cross, pheno.col = 3, model = "normal", method = "em", n.perm = 100)
# Heart rate
cross.scanHR <- scanone(cross, pheno.col = 4, model = "normal", method = "em")
cross.scanHR.perm <- scanone(cross, pheno.col = 4, model = "normal", method = "em", n.perm = 100)
# Heart weight
cross.scanHW <- scanone(cross, pheno.col = 6, model = "normal", method = "em")
cross.scanHW.perm <- scanone(cross, pheno.col = 6, model = "normal", method = "em", n.perm = 100)

# PLOT MAIN SCANS #
# Blood pressure
plot(cross.scanBP, main="Blood Pressure")
threshBP <- summary(cross.scanBP.perm, alpha = c(0.37, 0.10, 0.05))
abline(h = threshBP[1], col = "blue")
abline(h = threshBP[2], col = "red")
abline(h = threshBP[3], col = "green")
# Heart rate
plot(cross.scanHR, main="Heart Rate")
threshHR <- summary(cross.scanHR.perm, alpha = c(0.37, 0.10, 0.05))
abline(h = threshHR[1], col = "blue")
abline(h = threshHR[2], col = "red")
abline(h = threshHR[3], col = "green")
# Heart weight
plot(cross.scanHW, main="Heart Weight")
threshHW <- summary(cross.scanHW.perm, alpha = c(0.37, 0.10, 0.05))
abline(h = threshHW[1], col = "blue")
abline(h = threshHW[2], col = "red")
abline(h = threshHW[3], col = "green")

# SUMMARY #
# Blood pressure
summary(cross.scanBP, perm = cross.scanBP.perm, alpha = 0.05)
# Heart rate
summary(cross.scanHR, perm = cross.scanHR.perm, alpha = 0.05)
# Heart weight
summary(cross.scanHW, perm = cross.scanHW.perm, alpha = 0.05)

# EFFECT PLOTS #
# Blood pressure
firstBP <- find.marker(cross, chr = 3, pos = 48.7)
effectplot(cross, pheno.col = 3, mname1 = firstBP)
secondBP <- find.marker(cross, chr = 4, pos = 12.0)
effectplot(cross, pheno.col = 3, mname1 = secondBP)
# Heart rate
firstHR <- find.marker(cross, chr = 4, pos = 59.8)
effectplot(cross, pheno.col = 4, mname1 = firstHR)
# Heart weight
firstHW <- find.marker(cross, chr = 6, pos = 52.2)
effectplot(cross, pheno.col = 6, mname1 = firstHW)

# EOF