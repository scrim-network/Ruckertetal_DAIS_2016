##################################################################################################
#
#  -file = "LHS_plots_C.R"   Code written September 2014 edited March 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in the LHS workspace and creates supplementary Fig 1 and 2
#
#   -NOTE: The graphs will be saved as tif files.
#
##==============================================================================
## Copyright 2016 Kelsey Ruckert
## This file is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This file is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this file.  If not, see <http://www.gnu.org/licenses/>.
###################################################################################################

# Load in the saved workspace from LHS precalibration
load("Scratch/DAIS_precalibration_LHS_relative_2.RData")

# Install and open packages:
#install.packages('ash')
library(ash)
#install.packages('fields')
library(fields)
#install.packages('RColorBrewer')
library(RColorBrewer)
#install.packages('plotrix')
library(plotrix)

# Create color palettes
mypalette <- brewer.pal(9,"YlGnBu")
boxpalette <- brewer.pal(11,"Spectral")
expert_assessment <- brewer.pal(5,"Pastel2")

# Source plotting functions
source("Scripts/put_fig_letter.r")
source("Scripts/plot_rangefn.R")

#------------------------------------- Transparent Color Function -------------------------------------------
makeTransparent<- function(somecolor, alpha=100){
  someColor = someColor
  newColor<-col2rgb(someColor)
  apply(newColor,2 ,
        function(curcoldata)
        {rgb(red=curcoldata[1],
             green=curcoldata[2],
             blue=curcoldata[3], alpha=alpha,
             maxColorValue=255)})
}

#------------------------------------- Expert Assessments -------------------------------------------
# Set up previous study ranges for the year 2100.
# the previous ranges are in cm so divide by 100 to get meters
# Non-model ranges
pfeffer = c(12.8, 14.6, 61.9)/100 #Low2, low1, and high estimate in Pfeffer et al. 2008
Bamber = c(-2, 14, 83)/100 # 5%, Median, 95% estimates in Bamber and Aspinall 2013
IPCC_AR5 = c(-6, 4, 12)/100 # 5%, Median, 95% estimates in IPCC AR5

# Based on models with projections using RCP 8.5
Little = c(-8, 2.4, 13.3)/100 # 5%, Median, 95% estimates in Little et al. 2013
Kopp = c(-11, 4, 33)/100 # 5%, Median, 95% estimates in Kopp et al. 2014
Ritz = c(2, 11.9, 29.6)/100 # 5%, Median, 95% estimates in Ritz et al. 2015
Golledge = c(0.1, 0.39) # 'low' and 'high' simulations in Golledge et al. 2015
present = 2

model_assessment_colors = c("#000000", "#004949", "#009292", "#49E9BD")
nonmodel_assessment_colors = brewer.pal(5,"Blues")[3:5]

#-------------------------- Set widths and heights ------------------------------
inches_to_dpi = function(inch){ inch * 300 }

text_column_width   = 5.2
minimum_width       = 2.63
full_page_width     = 7.5
full_page_height    = 8.75
single_panel_height = 4

#-------------------------- Estimate Root Mean Square Error ------------------------------

lhs_mean <- rep(NA,enddate)
lhs_median <- rep(NA,enddate)
for(i in 1:enddate){
  lhs_mean[i] = mean(dais.pre.cali[,i])
    lhs_median[i] = median(dais.pre.cali[,i])
}

res_mean <-
res_median <-
res_optim <- rep(NA,length(obs.years[1:3]))
for(i in 1:length(obs.years[1:3])){
    res_mean[i] <- median(windows[i,]) - lhs_mean[obs.years[i]]
    res_median[i] <- median(windows[i,]) - lhs_median[obs.years[i]]
    #res_optim[i] <- median(windows[i,]) - (best.project[obs.years[i]]-mean(best.project[SL.1961_1990]))
}

res_mean.inst <- median(windows[4,]) - (lhs_mean[obs.years[4]] - lhs_mean[239992])
res_median.inst <- median(windows[4,]) - (lhs_median[obs.years[4]] - lhs_median[239992])
#res_optim[i] <- median(windows[i,]) - (best.project[obs.years[i]]-mean(best.project[SL.1961_1990]))

res_mean = c(res_mean, res_mean.inst)
res_median = c(res_median, res_median.inst)

# Estimate and return the root means square error
rmse_mean = sqrt(mean(res_mean^2))
rmse_median = sqrt(mean(res_median^2))
#rmse_optim = sqrt(mean(res_optim^2))

print(paste('mean rmse = ', rmse_mean))
print(paste('median rmse = ', rmse_median))
#print(paste('optim command without noise rmse = ', rmse_optim))

###################################### SUPPLEMENTARY FIGURES ############################################
#------------------------------------- Supplementary Figure 2 -------------------------------------------

#pdf(file="Figures/SuppFigures/suppFig4_dais_LHS.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
png(file="Scratch/S2_Fig_inst2_2.tif", family="Helvetica", width=text_column_width,
height=single_panel_height*2, units="in",pointsize=12, res=300)
par(mfrow=c(2,1),mgp=c(1.5,.5,0), mar=c(3.5,4,1,2)) # set figure dimensions

# Panel a) Hindcasts during the Last interglacial
plot(un.sflig$pdf, main="",lwd=3, col="gray91", xlab="Projected AIS volume loss", sub="during Last interglacial [SLE m]",
ylab="Probability Density",xlim=c(-10,25), ylim=c(-0.45,0.6), yaxt="n")
#ylab="Probability Density",xlim=c(min(un.sflig$pdf$x),20), ylim=c(-0.45,0.6), yaxt="n")

lines(LIG.sflig$pdf, col=mypalette[3], lwd=2)
lines(LGM.sflig$pdf, col=mypalette[5], lwd=2)
lines(MH.sflig$pdf, col=mypalette[1], lwd=2)
lines(present.sflig$pdf, col=mypalette[7], lwd=2)
lines(all.sflig$pdf, col=mypalette[9], lwd=2)

# Add box and whisker plots
probs = c(0.05, 0.95)
#where = c(-0.03, -0.08, -0.16, -0.24, -0.32, -0.4)
add.hor.box(all.prob_proj[,1], probs, width.size = 0.1, where.at = -0.05, tick.length = 0.01, line.width = 2, color = mypalette[9])
#points(dais.pre.cali[sur.all,obs.years[1]], -0.03, pch=15, col = mypalette[9], cex=0.75)
add.hor.box(present.prob_proj[,1], probs, width.size = 0.1, where.at = -0.13, tick.length = 0.01, line.width = 2, color = mypalette[7])
add.hor.box(MHprob_proj[,1], probs, width.size = 0.1, where.at = -0.21, tick.length = 0.01, line.width = 2, color = mypalette[1])
add.hor.box(LGMprob_proj[,1], probs, width.size = 0.1, where.at = -0.29, tick.length = 0.01, line.width = 2, color = mypalette[5])
add.hor.box(LIGprob_proj[,1], probs, width.size = 0.1, where.at = -0.37, tick.length = 0.01, line.width = 2, color = mypalette[3])
add.hor.box(unprob_proj[,1], probs, width.size = 0.1, where.at = -0.45, tick.length = 0.01, line.width = 2, color = "gray91")

abline(h=0.4, lty=2)
abline(v=0, lty=2, col="gray")
text(14, 0.4+0.1, cex=0.75, "Last Interglacial period\nData Constraint")
text(14, 0.4-0.1, cex=0.75, "Model Inversion\nThis study")

# Add range of observational constraint during the Last interglacial
plotrange(windows[1,2], (windows[1,1] + windows[1,2])/2, windows[1,1], year=F, height=0.55, color="black")

axis(3, labels = FALSE)
put.fig.letter(label="a.", location="topleft", font=2)
#=======
# Panel b) Projections in the year 2100
par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 1, 2))
plot(un.sf2100$pdf, main="",lwd=3, col="gray91", xlab="Projected AIS volume loss", sub= "in 2100 [SLE m]",
ylab="Probability Density",
xlim=c(-0.25,4.75), ylim=c(-5.5,14.25), yaxt="n") #ylim=c(-18,78)
#xlim=c(-3,1.75), ylim=c(-3,10.5))#, yaxt="n") #ylim=c(-18,78)

lines(LIG.sf2100$pdf, col=mypalette[3], lwd=2)
lines(LGM.sf2100$pdf, col=mypalette[5], lwd=2)
lines(MH.sf2100$pdf, col=mypalette[1], lwd=2)
lines(present.sf2100$pdf, col=mypalette[7], lwd=2)
lines(all.sf2100$pdf, col=mypalette[9], lwd=2)

# Add box and whisker plots
probs = c(0.05, 0.95)
#where = c(-0.2, -0.5, -1, -1.5, -2, -2.5)
add.hor.box(all.prob_proj[,6], probs, width.size = 1.2, where.at = -0.5, tick.length = 0.1, line.width = 2, color = mypalette[9])
#points(dais.pre.cali[sur.all,240100], -0.2, pch=15, col= mypalette[9], cex=0.75)
add.hor.box(present.prob_proj[,6], probs, width.size = 1.2, where.at = -1.5, tick.length = 0.1, line.width = 2, color = mypalette[7])
add.hor.box(MHprob_proj[,6], probs, width.size = 1.2, where.at = -2.5, tick.length = 0.1, line.width = 2, color = mypalette[1])
add.hor.box(LGMprob_proj[,6], probs, width.size = 1.2, where.at = -3.5, tick.length = 0.1, line.width = 2, color = mypalette[5])
add.hor.box(LIGprob_proj[,6], probs, width.size = 1.2, where.at = -4.5, tick.length = 0.1, line.width = 2, color = mypalette[3])
add.hor.box(unprob_proj[,6], probs, width.size = 1.2, where.at = -5.5, tick.length = 0.1, line.width = 2, color = "gray91")

abline(h=7, lty=2)
abline(v=0, lty=2, col="gray")
text(1.25, 7+0.5, cex=0.65, "Expert Assessments")
text(1.25, 7-1, cex=0.65, "Model Inversion\nThis study")

place.where = c(8.2, 9.2, 10.2, 11.2, 12.6, 13.6, 14.6)
#place.where = c(8.2,9.7,11.2,12.7,14.2)
width = 0.75  # Width of bars

rect(-2.3, 11.5, 4.3, 16.0, col = "gray93", border=NA)

# Add ranges from model and non-model expert assessments
# Model

polygon(x = c(Little[1], Little[3], Little[3], Little[1]),
y = c((place.where[1]), (place.where[1]),  (place.where[1] - width), (place.where[1]-width)),
border=NA, col = model_assessment_colors[1])
points(Little[2], (place.where[1]-0.2), col="black", pch="|", cex=0.8)

polygon(x = c(Kopp[1], Kopp[3], Kopp[3], Kopp[1]),
y = c((place.where[2]), (place.where[2]),  (place.where[2] - width), (place.where[2]-width)),
border=NA, col = model_assessment_colors[2])
points(Kopp[2], (place.where[2]-0.2), col="black", pch="|", cex=0.8)

polygon(x = c(Ritz[1], Ritz[3], Ritz[3], Ritz[1]),
y = c((place.where[3]), (place.where[3]),  (place.where[3] - width), (place.where[3]-width)),
border=NA, col = model_assessment_colors[3])
points(Ritz[2], (place.where[3]-0.2), col="black", pch="|", cex=0.8)

polygon(x = c(Golledge[1], Golledge[2], Golledge[2], Golledge[1]),
y = c((place.where[4]), (place.where[4]),  (place.where[4] - width), (place.where[4]-width)),
border=NA, col = model_assessment_colors[4])

# Non-model

polygon(x = c(IPCC_AR5[1], IPCC_AR5[3], IPCC_AR5[3], IPCC_AR5[1]),
y = c((place.where[5]), (place.where[5]),  (place.where[5] - width), (place.where[5]-width)),
border=NA, col = nonmodel_assessment_colors[1])
points(IPCC_AR5[2], (place.where[5]-0.2), col="black", pch="|", cex=0.8)

polygon(x = c(pfeffer[1], pfeffer[3], pfeffer[3], pfeffer[1]),
y = c((place.where[6]), (place.where[6]),  (place.where[6] - width), (place.where[6]-width)),
border=NA, col = nonmodel_assessment_colors[2])
points(pfeffer[2], (place.where[6]-0.2), col="black", pch="|", cex=0.8)

polygon(x = c(Bamber[1], Bamber[3], Bamber[3], Bamber[1]),
y = c((place.where[7]), (place.where[7]),  (place.where[7] - width), (place.where[7]-width)),
border=NA, col = nonmodel_assessment_colors[3])
points(Bamber[2], (place.where[7]-0.2), col="black", pch="|", cex=0.8)

axis(3, labels = FALSE)
put.fig.letter(label="b.", location="topleft", font=2)
box()

legend.names = c("Median & 90% C.I. (Little et al. 2013)","Median & 90% C.I. (Kopp et al. 2014)",
"Median & 90% C.I. (Ritz et al. 2015)", "'Low' & 'High' estimate (Golledge et al. 2015)",
"Median & 90% C.I. (Church et al. 2013)", "'Low2', 'Low1', &'High' estimate\n(Pfeffer et al. 2008)",
"Median & 90% C.I. (Bamber & Aspinall 2013)", paste("Unconstrained L.H. fits (", sample_length, ")", sep=""),
paste("LIG constraint fits (", length(surLIG), ")", sep=""),
paste("LGM constraint fits (", length(surLGM), ")", sep=""),
paste("MH constraint fits (", length(surMH), ")", sep=""),
paste("Instrumental constraint fits (", length(sur9311trend), ")", sep=""),
#paste("All L.H. constrained fits (0)", sep=""),# length(sur.all), ")", sep=""),
paste("All L.H. constrained fits (", length(sur.all), ")", sep=""))

legend("topright", legend.names, cex=0.575, pch=15, #bty="n",
col=c(model_assessment_colors, nonmodel_assessment_colors,"gray91", mypalette[3], mypalette[5], mypalette[1], mypalette[7],mypalette[9]))#,
#y.intersp=c(0.5,0.6,0.6,0.6, 0.6, 0.6,0.6,0.6,0.6,0.6,0.6))

dev.off()

#------------------------------------- Supplementary Figure 1 -------------------------------------------
# # LHS hindcasts & projection
# width=1920, height=1080
png(file="Scratch/S1_Fig_inst2.tif", family="Helvetica", width=text_column_width, height=single_panel_height*2, units="in",pointsize=12, res=300)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))
# jpeg(file="nSuppFig3_dais_mcmcLHS.jpeg", family="Helvetica", width=1590, height=1920, units="px", pointsize=40)
# par(mfrow=c(3,2), mgp=c(1.5,.5,0),mar=c(4, 4, 3, 2))

# Last interglacial 240 kyr Bp - 2010 AD
plot(date[1:240010], AIS_melt-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", lwd=1,
     xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]", ylim=c(-40,20), xaxt="n")
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
#Uncomment if there are hindcasts passing through the instrumental period
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
abline(h=0, lty=2, col="black")
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
lines(date, lhs_mean, col="black", lwd=1)
lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="gray", lwd=1.5)
lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="gray", lwd=1.5)
lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="gray", lwd=1.5)
ticks=c(-200000,-150000,-100000,-50000,0)
axis(side=1, at=ticks, labels=expression(-200,-150,-100,-50,0))
put.fig.letter(label="a.", location="topleft", font=2)

axis(1, at = date[kyrbp_25]:date[240010], labels = FALSE, col=boxpalette[1], lwd=2)

# Last glacial maximum 25 kyr Bp - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
plot(date[kyrbp_25:240010], AIS_melt[kyrbp_25:240010]-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", 
     lwd=1,xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]", ylim=c(-40,20), xaxt="n") # -30, 7
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
lines(date, lhs_mean, col="black", lwd=1)
# abline(h=0, lty=2, col="gray")
lines(c(last.glacialmax[2],last.glacialmax[2]), c(windows[2,1],windows[2,2]), col="gray", lwd=1.5)
lines(last.glacialmax, c(windows[2,1],windows[2,1],windows[2,1]), col="gray", lwd=1.5)
lines(last.glacialmax, c(windows[2,2],windows[2,2],windows[2,2]), col="gray", lwd=1.5)
ticks=c(-25000,-20000,-15000,-10000,-5000,0)
axis(side=1, at=ticks, labels=expression(-25,-20,-15,-10,-5,0))
put.fig.letter(label="b.", location="topleft", font=2)

box(col = boxpalette[1])
axis(1, at = date[kyrbp_6]:date[240010], labels = FALSE, col=boxpalette[2], lwd=2)

# Mid-Holocene6 kyr Bp - 2010 AD 
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

plot(date[kyrbp_6:240010], AIS_melt[kyrbp_6:240010]-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", lwd=1,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", ylim=c(-15,15)) #-8, 4
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
lines(date, lhs_mean, col="black", lwd=1)
lines(c(holocene[2],holocene[2]), c(windows[3,1],windows[3,2]), col="gray", lwd=1.5)
lines(holocene, c(windows[3,1],windows[3,1],windows[3,1]), col="gray", lwd=1.5)
lines(holocene, c(windows[3,2],windows[3,2],windows[3,2]), col="gray", lwd=1.5)
put.fig.letter(label="c.", location="topleft", font=2)

box(col = boxpalette[2])
axis(1, at = date[AD_1880]:date[240010], labels = FALSE, col=boxpalette[3], lwd=2)

# Present day 1880 - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))

plot(date[AD_1880:240010], Project_melt[AD_1880:240010]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=1,
     xlab="Date [yr AD]", ylab="AIS Volume loss [SLE m]", ylim=c(-0.6,0.6), xaxt="n") #(-0.10,0.05)
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
lines(date, lhs_mean, col="black", lwd=1)
#err_positive = windows[4,2]
#err_negative = windows[4,1]
#present = 2
#arrows(present, err_positive, present, err_negative, length=0, lwd=1.5, col="black")
#points(present,median(windows[4,]), col="black", pch=8, cex=0.75)
ticks=c(-120,-100,-50,0)
axis(side=1, at=ticks, labels=expression(1880,1900,1950,2000))
put.fig.letter(label="d.", location="topleft", font=2)

box(col = boxpalette[3])
axis(1, at = date[240000]:date[240020], labels = FALSE, col=boxpalette[4], lwd=2)

# Future 2010 - 2300 AD
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

plot(date[240000:enddate], Project_melt[240000:enddate]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=1,
     xlab="Date [yr AD]", ylab="AIS Volume loss [SLE m]", xlim=c(-10,290),ylim=c(-0.4,20), xaxt="n") #ylim=c(-0.5,5)
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
lines(date, lhs_mean, col="black", lwd=1)
#err_positive = windows[4,2]
#err_negative = windows[4,1]
#arrows(present, err_positive, present, err_negative, length=0, lwd=1.5, col="black")
ticks=c(0,50,100,150,200,250,300)
axis(side=1, at=ticks, labels=expression(2000,2050,2100,2150,2200,2250,2300))
put.fig.letter(label="e.", location="topleft", font=2)

box(col = boxpalette[4])

#plot.new()
plot(1, 1, col="white", xlab="", ylab="", xlim=c(-10,290),ylim=c(-10,20), xaxt="n", yaxt="n")
legend("topleft", expression(paste("±2 ", sigma, " Data constraint")), pch="I", bty="n", col="gray", cex=0.85)
legend("top", c("", "Optim best-fit", "Mean", "No Constraints", "Last interglacial constraint", "Last glacial maximum constraint","Mid-Holocene constraint",
                   "Instrumental period (1992-2011)
constraint", "All constraints"),lty=c(NA,1,1,NA,NA,NA,NA,NA,NA,NA), pch=c(NA,NA,NA,15,15,15,15,15,15,15),
       col=c("white", "coral3", "black", "gray91", mypalette[3], mypalette[5], mypalette[1], mypalette[7], mypalette[9]), bty="n", cex=0.85)
dev.off()

###################################### END ############################################

