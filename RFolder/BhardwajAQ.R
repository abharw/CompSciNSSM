# Author: Arav Bhardwaj
# Date: 09-08-22
# AirQuality.R
# Analyzes airquality data

# SETUP
rm(list=ls())
setwd("C:\Users\abhardwaj24\Desktop\RFolder")

# CODE

# Read in CSV
airquality <- read.csv("airquality.csv", header = T)

# Print and save names to .txt file
sink("names.txt")
names(airquality)
sink()

# Print and save summary statistics to a .tx file
sink("summary.txt")
summary(airquality)
sink()

# Perform linear regression on wind and temperature (x-values) and
# ozone (y-value) and save to a .txt file
multi_fit <- lm(airquality$Ozone ~ airquality$Wind + airquality$Temp)
sink("multireg.txt")
summary(multi_fit)
sink()

# Perform and plot regression on temperature (x-value) and ozone (y-value)
# and save to a .jpg file
linear_fit <- lm(airquality$Ozone ~ airquality$Temp)
jpeg("regressionplot.jpg")
plot(airquality$Ozone ~ airquality$Temp, main="Plot of Temp vs. Ozone",
     ylab="Ozone",
     xlab="Temp")
abline(linear_fit)
dev.off()

# Plot four boxplots in a 2 x 2 matrix and save to a .jpg file
jpeg("boxplots.jpg")
# Create a 2 x 2 plotting matrix
par(mfrow = c(2, 2))
boxplot(airquality$Wind, main="Wind", horizontal = T)
boxplot(airquality$Temp, main="Temp", horizontal = T)
boxplot(airquality$Solar.R, main="Solar Radiation", horizontal = T)
boxplot(airquality$Ozone, main="Ozone", horizontal = T)
# Reset plotting matrix
par(mfrow = c(1, 1))
dev.off()

# Calculate and save BIC scores to a .txt file
# Includes analyses for each BIC score
sink("bicscores.txt")
bic0 <- BIC(lm(airquality$Ozone ~ 1)) # 1148.801
# This score is a baseline for evaluating the subsequent scores
bic1 <- BIC(lm(airquality$Ozone ~ airquality$Temp))  # 1075.967
sprintf("BIC Temp -> Ozone: %f", bic1)
# Ozone is affected by temperature, as the difference between
# this BIC score and the baseline is much greater than 10
bic2 <- BIC(lm(airquality$Ozone ~ airquality$Solar.R))  # 1091.84
sprintf("BIC Solar.R -> Ozone: %f", bic2)
# There is also a relationship between ozone and solar radiation, as there is a significant
# difference between this BIC score and the baseline. However, the difference is not as great
# as the previous BIC score, indicating that there is a weaker relationship between ozone and
# solar radiation as compared to ozone and temperature
BIC(lm(airquality$Solar.R ~ 1)) # 1737.428
# This score is a bseline for evaluating the next score
bic3 <- BIC(lm(airquality$Solar.R ~ airquality$Ozone))  # 1315.552
sprintf("BIC Ozone -> Solar.R: %f", bic3)
# Solar radiation is greatly affected by ozone, as the difference
# between this BIC score and the baseline is quite large
bic4 <- BIC(lm(airquality$Temp ~ airquality$Solar.R)) + BIC(lm(airquality$Ozone ~ airquality$Temp)) # 2141.483
sprintf("BIC Solar.R -> Temp -> Ozone: %f", bic4)
# This BIC score is much higher than the previous scores pertaining to ozone,
# suggesting that looking at solar radiation and temperature together is more significant.
bic5 <- BIC(lm(airquality$Temp ~ airquality$Solar.R)) + BIC(lm(airquality$Ozone ~ airquality$Solar.R)) + BIC(lm(airquality$Ozone ~ airquality$Temp))
sprintf("BIC Ozone <- Solar.R -> Temp -> Ozone: %f", bic5)
# This score tells us that ozone is greatly affected by solar radiation and
# temperature (which is also affected by solar radiation)
sink()

# EOF