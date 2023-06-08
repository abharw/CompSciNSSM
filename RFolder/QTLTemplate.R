# R script for analyzing Sugiyama QTL data
# Arav Bhardwaj
# hyper.qtl script
# October 15, 2022
#
# clean things up
rm(list=ls())
setwd("/Users/abhardwaj24/Desktop/RFolder")
#
# load the QTL library
# NOTE!  I first had to INSTALL the library using: install.packages("qtl")
install.packages("qtl")
# Now I can use the package qtl
library(qtl)

# Load the data!
cross <- read.cross("csv", file = "sugiyamashort.csv", genotypes = c("C", "H", "B"),
                    na.strings = "-",alleles = c("C", "B"))
# Sometimes the genetic markers are too close. Jittermap will move them apart slightly so my results are better.
jittermap(cross)
#  A summary of the cross gives me some basic data.  Nice!
summary(cross)
# I need to see what phenotypes are in the dataset.  names does that for me
names(cross$pheno)
# take a look at my data, make sure it's pretty clean
#  I should NOT see any really big red spots, especially in the bottom right corner under the diagonal
cross <- est.rf(cross)
plotRF(cross)
# It's nice to see my genetic map -- all of the horizontal lines are genetic markers that have been inserted.
plot.map(cross)
#  It's often the case that I have missing data -- plot.missing shows me where it is
plotMissing(cross)
# a histogram of my BP phenotype -- if I get a "normal" distribution (a bell-shaped curve), then I can be pretty confident of the data
hist(cross$pheno$BP_final)
#  I'm lazy....I hate to type "cross$pheno$bp every time I need to study the BP, so I give it a short name.
bp <- cross$pheno$BP_final
#  Another diagnostic....if my data is relatively clean, I should get a nice 45 degree diagonal line
qqnorm(bp)
qqline(bp)
# Now I'm going to generate a mainscan.  First, I calculate what the scan SHOULD look like, so I'm going to calculate a genetic probability map.
cross <- calc.genoprob(cross, step = 2.0, off.end = 0, error.prob = 1.0e-4,
                       map.function = "haldane", stepwidth = "fixed")
# Run a simulated geno probability calculations
cross <- sim.geno(cross, step = 2.0, off.end = 0, error.prob = 1.0e-4,
                       map.function = "haldane", stepwidth = "fixed")
# Perform the mainscan for the QTL
cross.scanBP <- scanone(cross, pheno.col = 1, model = "normal", method = "em")
#  I'm only going to run this for 100 "permulations" -- typically you do 500-1000, but that takes a LONG time.
cross.scanBP.perm <- scanone(cross, pheno.col = 1, model = "normal", method = "em", n.perm = 100)
# plot the mainscan
plot(cross.scanBP, main="Mainscan plot of BP")
#  I'm putting threshold likes at 63% confidence, 90% confidence, and 95% confidence.
thresh <- summary(cross.scanBP.perm, alpha = c(0.37, 0.10, 0.05))
abline(h = thresh[1], col = "blue")
abline(h = thresh[2], col = "red")
abline(h = thresh[3], col = "green")
#  I'd like to see a text-based output of my scan
summary(cross.scanBP, perm = cross.scanBP.perm, alpha = 0.05)

# do an effect plot
# once you see an effect plot, you'll understand what it does!


# second effect plot


#
#  I'm going to do a confidence interval plot, try to zoom in on the scan



# All done
#EOF