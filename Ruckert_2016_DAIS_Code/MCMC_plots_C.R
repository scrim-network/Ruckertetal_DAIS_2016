##################################################################################################
#
#  -file = "MCMC_plots.R"   Code written September 2014 edited March 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in the MCMC workspace
#
#   -NOTE: The graphs will be saved as tif files in the current working directory.
#
###################################################################################################

load("Scratch/Workspace/DAIS_MCMC_R_C_calibration_relative_8e5.RData")
#load("Scratch/Workspace/DAIS_MCMC_R_C_calibration_realtive.RData")
#load("Scratch/Workspace/DAIS_MCMC_R_C_calibration.RData") # Load in the saved workspace from MCMC calibration
#load("Scratch/Workspace/DAIS_MCMC_Matlabcalibration_1234_C.RData")

#install.packages('ash')
library(ash)
#install.packages('fields')
library(fields)
#install.packages('RColorBrewer')
library(RColorBrewer)
#install.packages('plotrix')
library(plotrix)

mypalette <- brewer.pal(9,"YlGnBu")
boxpalette <- brewer.pal(11,"Spectral")
expert_assessment <- brewer.pal(5,"Pastel2")

source("Scripts/put_fig_letter.r")
source("Scripts/plot_rangefn.R")
source('Scripts/plot_PdfCdfSf.R')

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

someColor = c("red", "goldenrod2", "blue", "darkorange2")
studies.color <- makeTransparent(someColor, 100)
someColor = c("goldenrod2")
this.study.color <- makeTransparent(someColor, 150)
this.study.color_95 <- makeTransparent(someColor, 100)
this.study.color_99 <- makeTransparent(someColor, 75)

someColor = "darkorange2"
nonoise.color_68 <- makeTransparent(someColor, 150)
nonoise.color_90 <- makeTransparent(someColor, 100)
nonoise.color_95 <- makeTransparent(someColor, 75)

#-------------------------- Set widths and heights ------------------------------
inches_to_dpi = function(inch){ inch * 300 }

text_column_width   = 5.2
minimum_width       = 2.63
full_page_width     = 7.5
full_page_height    = 8.75
single_panel_height = 4

#------------------------------------- Expert Assessments -------------------------------------------
#Set up previous study ranges for the year 2100.
#the previous ranges are in cm so divide by 100 to get meters
# Non-model ranges
pfeffer = c(12.8, 14.6, 61.9)/100 #Low2, low1, and high estimate in Pfeffer et al. 2008
Bamber = c(-2, 14, 83)/100 # 5%, Median, 95% estimates in Bamber and Aspinall 2013
# IPCC_AR5 = c(-15, 4, 23)/100 # 5%, Median, 95% estimates in IPCC AR5
IPCC_AR5 = c(-6, 4, 12)/100 # 5%, Median, 95% estimates in IPCC AR5

# Based on models with projections using RCP 8.5
Little = c(-8, 2.4, 13.3)/100 # 5%, Median, 95% estimates in Little et al. 2013
Kopp = c(-11, 4, 33)/100 # 5%, Median, 95% estimates in Kopp et al. 2014
Ritz = c(2, 11.9, 29.6)/100 # 5%, Median, 95% estimates in Ritz et al. 2015
Golledge = c(0.1, 0.39) # 'low' and 'high' simulations in Golledge et al. 2015
present = 2

model_assessment_colors = c("#000000", "#004949", "#009292", "#49E9BD")
nonmodel_assessment_colors = brewer.pal(5,"Blues")[3:5]

#------------------------------------- Calculate 90% and 99% credible intervals -------------------------------------------
mcmc_5 <-
  mcmc_95 <-
  mcmc_2point5 <-
  mcmc_975 <-
#mcmc_point5 <-
  mcmc_std <-
  mcmc_mean<-
#mcmc_995 <-
  mcmc_5_NN <-
  mcmc_95_NN <-
  mcmc_2point5_NN <-
  mcmc_975_NN <-
#mcmc_point5_NN <-
  mcmc_std_NN <-
  mcmc_mean_NN<-rep(NA,enddate)
#mcmc_995_NN <- # 240300 years
for(i in 1:enddate){
  mcmc_5[i] = quantile(proj.mcmc.1961_1990[,i],0.05) # 90%
  mcmc_95[i] = quantile(proj.mcmc.1961_1990[,i],0.95)
  
  mcmc_2point5[i] = quantile(proj.mcmc.1961_1990[,i],0.025) # 95%
  mcmc_975[i] = quantile(proj.mcmc.1961_1990[,i],0.975)
  
  #mcmc_point5[i] = quantile(proj.mcmc.1961_1990[,i],0.005) # 99%
  #mcmc_995[i] = quantile(proj.mcmc.1961_1990[,i],0.995)
  
  mcmc_std[i] = sd(proj.mcmc.1961_1990[,i]) # 1 standard deviation
  mcmc_mean[i] = mean(proj.mcmc.1961_1990[,i]) #mean
  
  mcmc_5_NN[i] = quantile(proj.mcmc.anomaly[,i],0.05) # 90%
  mcmc_95_NN[i] = quantile(proj.mcmc.anomaly[,i],0.95)
  
  mcmc_2point5_NN[i] = quantile(proj.mcmc.anomaly[,i],0.025) # 95%
  mcmc_975_NN[i] = quantile(proj.mcmc.anomaly[,i],0.975)
  
  #mcmc_point5_NN[i] = quantile(proj.mcmc.anomaly[,i],0.005) # 99%
  #mcmc_995_NN[i] = quantile(proj.mcmc.anomaly[,i],0.995)
  
  mcmc_std_NN[i] = sd(proj.mcmc.anomaly[,i]) # 1 standard deviation
  mcmc_mean_NN[i] = mean(proj.mcmc.anomaly[,i]) #mean
}

x_mcmc_90=c(mcmc_5, rev(mcmc_95)); y_mcmc_90=c(date, rev(date))
x_mcmc_95=c(mcmc_2point5, rev(mcmc_975)); y_mcmc_95=c(date, rev(date))
#x_mcmc_99=c(mcmc_point5, rev(mcmc_995)); y_mcmc_99=c(date, rev(date))

x_mcmc_90_NN=c(mcmc_5_NN, rev(mcmc_95_NN)); y_mcmc_90_NN=c(date, rev(date))
x_mcmc_95_NN=c(mcmc_2point5_NN, rev(mcmc_975_NN)); y_mcmc_95_NN=c(date, rev(date))
#x_mcmc_99_NN=c(mcmc_point5_NN, rev(mcmc_995_NN)); y_mcmc_99_NN=c(date, rev(date))

mcmc_plusSTD = mcmc_mean + mcmc_std
mcmc_minusSTD = mcmc_mean - mcmc_std
plus_minusSTD = c(mcmc_minusSTD, rev(mcmc_plusSTD))

mcmc_plusSTD_NN = mcmc_mean_NN + mcmc_std_NN
mcmc_minusSTD_NN = mcmc_mean_NN - mcmc_std_NN
plus_minusSTD_NN = c(mcmc_minusSTD_NN, rev(mcmc_plusSTD_NN))

y_STD_2100=c(date[1:240100], rev(date[1:240100]))
plus_minusSTD_2100 = c(mcmc_minusSTD[1:240100], rev(mcmc_plusSTD[1:240100]))
plus_minusSTD_2100_NN = c(mcmc_minusSTD_NN[1:240100], rev(mcmc_plusSTD_NN[1:240100]))

#-------------------------- Estimate Root Mean Square Error ------------------------------

mcmc_median <- rep(NA,enddate)
for(i in 1:enddate){
    mcmc_median[i] = median(proj.mcmc.1961_1990[,i])
}

res_mean <-
res_median <- rep(NA,length(obs.years[1:3]))
for(i in 1:length(obs.years)){
    res_mean[i] <- median(windows[i,]) - mcmc_mean[obs.years[i]]
    res_median[i] <- median(windows[i,]) - mcmc_median[obs.years[i]]
}

res_mean.inst <- median(windows[4,]) - (mcmc_mean[obs.years[4]] - mcmc_mean[239992])
res_median.inst <- median(windows[4,]) - (mcmc_median[obs.years[4]] - mcmc_median[239992])

res_mean = c(res_mean, res_mean.inst)
res_median = c(res_median, res_median.inst)

# Estimate and return the root means square error
rmse_mean = sqrt(mean(res_mean^2))
rmse_median = sqrt(mean(res_median^2))

print(paste('mean rmse = ', rmse_mean))
print(paste('median rmse = ', rmse_median))

print("2100:")
print(paste('mean 2100 = ', mcmc_mean[240100]))
print(paste('median 2100 = ', mcmc_median[240100]))
print(paste('5% = ', mcmc_5[240100], ' and 95% =' , mcmc_95[240100]))
print(paste('mean - sigma = ', plus_minusSTD[240100], ' mean + sigma =' , mcmc_plusSTD[240100]))

print("LIG:")
print(paste('mean LIG = ', mcmc_mean[obs.years[1]]))
print(paste('median LIG = ', mcmc_median[obs.years[1]]))
print(paste('5% = ', mcmc_5[obs.years[1]], ' and 95% =' , mcmc_95[obs.years[1]]))
print(paste('2.5% = ', mcmc_2point5[obs.years[1]], ' and 97.5% =' , mcmc_975[obs.years[1]]))
print(paste('mean - sigma = ', plus_minusSTD[obs.years[1]], ' mean + sigma =' , mcmc_plusSTD[obs.years[1]]))

#------------------------------------- Create mean estimate projection and output R, H, & F ----------------------
source("Scripts/iceflux.mult_func_outRHF_C.R")
#Project the best fit from the fitted parameters:
standards = c(Tf, rho_w, rho_i, rho_m, Toc_0, Rad0, Volo)
RHF_outputs = iceflux_RHF(mean.parameters, project.forcings, standards)
mean_RHF = RHF_outputs$SLE[90000:140000]-mean(RHF_outputs$SLE[SL.1961_1990])
#mean_RHF = mean_RHF + rnorm(length(mean_RHF), mean=0, sd=bias.mean)
mean_RHF = mean_RHF + rnorm(length(mean_RHF), mean=0, sd=mean(paleo.bias))
  
###################################### MAIN FIGURES ############################################
#------------------------------------- Figure 1 -------------------------------------------
#pdf(file="Figures/SuppFigures/Ruckertetal_daisRHF_Sfig2a.pdf", family="Helvetica", height=5.4, width=6.7,pointsize=11)
png(file="Scratch/Figures/Fig1_C_instPaleo.tif", family="Helvetica", width=6,
    height=6, units="in",pointsize=12, res=300)
par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1)) # set figure dimensions
# Last interglacial 240 kyr Bp - 2010 AD
plot(date[90000:140000], mcmc_mean[90000:140000], #RHF_outputs$SLE[90000:140000]-mean(RHF_outputs$SLE[SL.1961_1990]),
     typ="l", col="black", lwd=2, xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]",
     ylim=c(-1,6), xaxt="n")
abline(h=0, lty=2, col="black")

lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="black", lwd=2)
lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="black", lwd=2)
lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="black", lwd=2)
points(last.interglacial[2],median(windows[1,]), col="black", pch=8, cex=0.585)
ticks=c(-150000,-140000,-130000,-120000,-110000,-100000)
axis(side=1, at=ticks, labels=expression(-150, -140, -130,-120,-110,-100))
put.fig.letter(label="a.", location="topleft", font=2)

plot(date[90000:140000], RHF_outputs$Rad[90000:140000], col="black", lwd=2,typ="l",
     xlab="Date [Kyr BP]", ylab="Radius [m]", xaxt="n")
ticks=c(-150000,-140000,-130000,-120000,-110000,-100000)
axis(side=1, at=ticks, labels=expression(-150, -140, -130,-120,-110,-100))
put.fig.letter(label="b.", location="topleft", font=2)

plot(date[90000:140000], RHF_outputs$WatDepth[90000:140000], col="black", lwd=2,typ="l",
     xlab="Date [Kyr BP]", ylab="Water depth at the grounding line 
(H) [m]", xaxt="n")
ticks=c(-150000,-140000,-130000,-120000,-110000,-100000)
axis(side=1, at=ticks, labels=expression(-150, -140, -130,-120,-110,-100))
put.fig.letter(label="c.", location="topleft", font=2)

plot(date[90000:140000], RHF_outputs$IceFlux[90000:140000], col="black", lwd=2,typ="l",
     xlab="Date [Kyr BP]", ylab="Ice flux at the grounding line 
(F) [m/yr]", xaxt="n")
ticks=c(-150000,-140000,-130000,-120000,-110000,-100000)
axis(side=1, at=ticks, labels=expression(-150, -140, -130,-120,-110,-100))
put.fig.letter(label="d.", location="topleft", font=2)
dev.off()

#------------------------------------- Figure 2 -------------------------------------------
# MCMC hindcasts & projection
png(file="Scratch/Figures/Fig2_C_instPaleo.tif", family="Helvetica", width=text_column_width,
    height=single_panel_height*2, units="in",pointsize=12, res=300)
# pdf(file="Figures/Ruckertetal_dais90MCMC_Fig1.pdf", family="Helvetica", width=6.7, height=8.1, pointsize=12)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

#pdf(file="NEWnFig1_dais_mcmcLHS.pdf", family="Helvetica",height=2.7, width=6.7,pointsize=11)
#par(mfrow=c(1,2),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1)) # set figure dimensions
# Last interglacial 240 kyr Bp - 2010 AD
plot(date[1], proj.mcmc.1961_1990[1,1], type="l", col="powderblue", lwd=2,
     xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(date[1], date[235000]), ylim=c(-25,10), xaxt="n")

polygon(y_mcmc_95, x_mcmc_95, col=this.study.color_99, border=NA)
polygon(y_mcmc_90, x_mcmc_90, col=this.study.color_95, border=NA)
polygon(y_mcmc_90, plus_minusSTD, col=this.study.color, border=NA)

abline(h=0, lty=3, col="black")

lines(date, mcmc_mean, col="brown")

ticks=c(-200000,-150000,-100000,-50000,0)
axis(side=1, at=ticks, labels=expression(-200,-150,-100,-50,0))
put.fig.letter(label="a.", location="topleft", font=2)

axis(1, at = date[90000]:date[240010], labels = FALSE, col=boxpalette[1], lwd=2)

# Last interglacial 150 kyr Bp - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))

plot(date[90000], proj.mcmc.1961_1990[1,90000], type="l", col="powderblue", 
     lwd=2,xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]",
     xlim=c(date[90000], date[235000]), ylim=c(-5,6), xaxt="n")

polygon(y_mcmc_95, x_mcmc_95, col=this.study.color_99, border=NA)
polygon(y_mcmc_90, x_mcmc_90, col=this.study.color_95, border=NA)
polygon(y_mcmc_90, plus_minusSTD, col=this.study.color, border=NA)

lines(date, mcmc_mean, col="brown")
lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="black", lwd=2)
lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="black", lwd=2)
lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="black", lwd=2)
points(last.interglacial[2],median(windows[1,]), col="black", pch=8, cex=0.585)
ticks=c(-150000,-125000,-100000, -75000, -50000, -25000, 0)
axis(side=1, at=ticks, labels=expression(-150, -125, -100,-75,-50,-25, 0))
put.fig.letter(label="b.", location="topleft", font=2)

box(col = boxpalette[1])
axis(1, at = date[kyrbp_25]:date[240010], labels = FALSE, col=boxpalette[1], lwd=2)

# Last glacial maximum 25 kyr Bp - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

#pdf(file="NEWnFig25kyr.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
#par(mfrow=c(1,1),mgp=c(1.5,.5,0)) # set figure dimensions
plot(date[kyrbp_25], proj.mcmc.1961_1990[1,kyrbp_25], type="l", col="powderblue", 
     lwd=2,xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]",
     xlim=c(date[kyrbp_25], date[238000]), ylim=c(-21,0.5), xaxt="n")

polygon(y_mcmc_95, x_mcmc_95, col=this.study.color_99, border=NA)
polygon(y_mcmc_90, x_mcmc_90, col=this.study.color_95, border=NA)
polygon(y_mcmc_90, plus_minusSTD, col=this.study.color, border=NA)

lines(date, mcmc_mean, col="brown")
lines(c(last.glacialmax[2],last.glacialmax[2]), c(windows[2,1],windows[2,2]), col="black", lwd=2)
lines(last.glacialmax, c(windows[2,1],windows[2,1],windows[2,1]), col="black", lwd=2)
lines(last.glacialmax, c(windows[2,2],windows[2,2],windows[2,2]), col="black", lwd=2)
points(last.glacialmax[2],median(windows[2,]), col="black", pch=8, cex=0.585)
ticks=c(-25000,-20000,-15000,-10000,-5000,0)
axis(side=1, at=ticks, labels=expression(-25,-20,-15,-10,-5,0))
put.fig.letter(label="c.", location="topleft", font=2)

box(col = boxpalette[1])
axis(1, at = date[kyrbp_6]:date[240010], labels = FALSE, col=boxpalette[2], lwd=2)

# Mid-Holocene6 kyr Bp - 2010 AD 
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
#par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

#pdf(file="NEWnFig6kyr.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
#par(mfrow=c(1,1),mgp=c(1.5,.5,0)) # set figure dimensions
plot(date[kyrbp_6], proj.mcmc.1961_1990[1,kyrbp_6], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(date[kyrbp_6], date[240010]),ylim=c(-6,3)) #ylim=c(-5,0.5))

polygon(y_mcmc_95, x_mcmc_95, col=this.study.color_99, border=NA)
polygon(y_mcmc_90, x_mcmc_90, col=this.study.color_95, border=NA)
polygon(y_mcmc_90, plus_minusSTD, col=this.study.color, border=NA)

lines(date, mcmc_mean, col="brown")
lines(c(holocene[2],holocene[2]), c(windows[3,1],windows[3,2]), col="black", lwd=2)
lines(holocene, c(windows[3,1],windows[3,1],windows[3,1]), col="black", lwd=2)
lines(holocene, c(windows[3,2],windows[3,2],windows[3,2]), col="black", lwd=2)
points(holocene[2],median(windows[3,]), col="black", pch=8, cex=0.585)
put.fig.letter(label="d.", location="topleft", font=2)

box(col = boxpalette[2])
axis(1, at = date[AD_1880]:date[240010], labels = FALSE, col=boxpalette[3], lwd=2)

# Present day 1880 - 2010 AD
#par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

#pdf(file="NEWnFig1880.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
#par(mfrow=c(1,1),mgp=c(1.5,.5,0)) # set figure dimensions
plot(date[AD_1880], proj.mcmc.1961_1990[1,AD_1880], type="l", col="powderblue", lwd=2,
     xlab="Date [yr AD]", ylab="AIS Volume loss [SLE m]",
     xlim=c(date[AD_1880], date[240010]), xaxt="n", ylim=c(-0.1,0.1)) #ylim=c(-2,3))

polygon(y_mcmc_95, x_mcmc_95, col=this.study.color_99, border=NA)
polygon(y_mcmc_90, x_mcmc_90, col=this.study.color_95, border=NA)
polygon(y_mcmc_90, plus_minusSTD, col=this.study.color, border=NA)

lines(date,proj.mcmc.1961_1990[585,] , col="lightgray", lwd=0.95)
lines(date,proj.mcmc.1961_1990[1225,] , col="gray", lwd=0.95)

lines(date, mcmc_mean, col="brown")
#err_positive = windows[4,2]
#err_negative = windows[4,1]
#present = 2
#arrows(present, err_positive, present, err_negative, length=0, lwd=2, col="black")
#points(present,median(windows[4,]), col="black", pch=8, cex=0.585)
ticks=c(-120,-100,-50,0)
axis(side=1, at=ticks, labels=expression(1880,1900,1950,2000))
put.fig.letter(label="e.", location="topleft", font=2)

box(col = boxpalette[3])
axis(1, at = date[240000]:date[240020], labels = FALSE, col=boxpalette[4], lwd=2)

# Future 2010 - 2300 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
plot(date[240000], proj.mcmc.1961_1990[1,240000], type="l", col="powderblue", lwd=2,
     xlab="Date [yr AD]", ylab="AIS Volume loss [SLE m]",
     xlim=c(-10,285), xaxt="n", ylim=c(-0.1,7)) #ylim=c(-2,7))

polygon(y_mcmc_95, x_mcmc_95, col=this.study.color_99, border=NA)
polygon(y_mcmc_90, x_mcmc_90, col=this.study.color_95, border=NA)
polygon(y_mcmc_90, plus_minusSTD, col=this.study.color, border=NA)

lines(date,proj.mcmc.1961_1990[585,] , col="lightgray", lwd=0.95)
lines(date,proj.mcmc.1961_1990[1225,] , col="gray", lwd=0.95)

lines(date, mcmc_mean, col="brown")
#err_positive = windows[4,2]
#err_negative = windows[4,1]
#arrows(present, err_positive, present, err_negative, length=0, lwd=2, col="black")
ticks=c(0,50,100,150,200,250,300)
axis(side=1, at=ticks, labels=expression(2000,2050,2100,2150,2200,2250,2300))
put.fig.letter(label="f.", location="topleft", font=2)

box(col = boxpalette[4])

legend("topleft", expression(paste("±2 ", sigma, " Data constraint")), pch="I", bty="n", col="black")
legend("topleft", c("", "Mean hindcast & projection", "95% credible interval", "90% credible interval",
expression(paste("±1 ", sigma, " credible interval"))),
lty=c(NA,1,NA,NA,NA), pch=c(NA,NA,15,15,15), bty="n",
       col=c("white", "brown", this.study.color_99, this.study.color_95, this.study.color))

dev.off()

#------------------------------------- Figure 3 -------------------------------------------
# Function to find SLE values in certain years 'fn.prob.proj'
# Calculate the pdf, cdf, and sf of AIS melt estimates in:
NN_mcmc.prob_proj <- fn.prob.proj(proj.mcmc.anomaly, year.pcs, subset_length, un.constr=T)
NN_LIG.mcmc <- plot.sf(NN_mcmc.prob_proj[,1], make.plot=F) # 120,000 BP (Last interglacial)
NN_2100.mcmc <- plot.sf(NN_mcmc.prob_proj[,6], make.plot=F) # 2100

#pdf(file="Figures/Ruckertetal_daisLIG2100pdf_Fig2aa.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
png(file="Scratch/Figures/Fig3_C_instPaleo.tif", family="Helvetica", width=text_column_width,
    height=single_panel_height*2, units="in",pointsize=12, res=300)
par(mfrow=c(2,1),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1)) # set figure dimensions

plot(LIG.sf.mcmc$pdf, main="",lwd=3, col=this.study.color, xlab="Hindcasted AIS Volume loss", sub="during Last interglacial [SLE m]",
ylab="Probability Density", ylim=c(-0.07,0.46), yaxt="n",xlim=c(min(LIG.sf.mcmc$pdf$x),17))

lines(LIG.sf.mcmc$pdf, col=this.study.color, lwd=2)
#lines(NN_LIG.mcmc$pdf, col=nonoise.color_68, lwd=2)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,1], probs, width.size = 0.07, where.at = -0.05, tick.length = 0.01, line.width = 2, color = this.study.color)
#add.hor.box(NN_mcmc.prob_proj[,1], probs, width.size = 0.07, where.at = -0.12, tick.length = 0.01, line.width = 2, color = nonoise.color_68)

# abline(h=1, lty=2)
abline(h=0.4, lty=2)
text(11, 0.4+0.06, cex=0.8, "Last Interglacial period
Data Constraint")
text(11, 0.4-0.06, cex=0.8, "Model Inversion
This study")

plotrange(windows[1,2], (windows[1,1] + windows[1,2])/2, windows[1,1], year=F, height=0.45, color="black")
put.fig.letter(label="a.", location="topleft", font=2)

#par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 1, 2))
plot(sf.2100.mcmc$pdf, main="", col=this.study.color,lwd=2,
     xlab="Projected AIS Volume loss", sub="in 2100 [SLE m]", 
     ylab="Probability Density", 
xlim=c(-2.2,1.5),
ylim=c(-0.4,3.5), yaxt="n") #ylim=c(-7,53), yaxt="n") #2.4

lines(sf.2100.mcmc$pdf, col=this.study.color, lwd=2)
#lines(NN_2100.mcmc$pdf, col=nonoise.color_68, lwd=2)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,6], probs, width.size = 0.5, where.at = -.25, tick.length = 0.06, line.width = 2, color = this.study.color)
#add.hor.box(NN_mcmc.prob_proj[,6], probs, width.size = 1, where.at = -1.3, tick.length = 0.2, line.width = 2, color = nonoise.color_68)

# abline(h=30, lty=2)
# text(0.6, 30+3, cex=0.585, "Expert Assessments")
# text(0.6, 30-5, cex=0.585, "Model Inversion
#      This study")
# plotrange(Little[1], Little[2], Little[3], year=F, height=34, color="pink")
# plotrange(IPCC_AR5[1], IPCC_AR5[2], IPCC_AR5[3], year=F, height=38, color="purple")
# plotrange(Kopp[1], Kopp[2], Kopp[3], year=F, height=42, color="orange")
# plotrange(pfeffer[1], pfeffer[2], pfeffer[3], year=F, height=46, color="gray")
# plotrange(Bamber[1], Bamber[2], Bamber[3], year=F, height=50, color="red")

# abline(h=5, lty=2)
abline(h=1.9, lty=2)
text(1, 1.9+0.25, cex=0.8, "Expert Assessments")
text(1, 1.9-0.3, cex=0.8, "Model Inversion
This study")
#plotrange(Little[1], Little[2], Little[3], year=F, height=1.35, color="pink")
#plotrange(IPCC_AR5[1], IPCC_AR5[2], IPCC_AR5[3], year=F, height=1.58, color="purple")
#plotrange(Kopp[1], Kopp[2], Kopp[3], year=F, height=1.81, color="orange")
#plotrange(pfeffer[1], pfeffer[2], pfeffer[3], year=F, height=2.04, color="gray")
#plotrange(Bamber[1], Bamber[2], Bamber[3], year=F, height=2.27, color="red")

place.where = c(2.36, 2.65, 2.88, 3.11, 3.34)
width = 0.2  # Width of bars

polygon(x = c(Little[1], Little[3], Little[3], Little[1]),
y = c((place.where[1]), (place.where[1]),  (place.where[1] - width), (place.where[1]-width)),
border=NA, col = expert_assessment[1])
points(Little[2], (place.where[1]-0.1), col="black", pch="|")

polygon(x = c(IPCC_AR5[1], IPCC_AR5[3], IPCC_AR5[3], IPCC_AR5[1]),
y = c((place.where[2]), (place.where[2]),  (place.where[2] - width), (place.where[2]-width)),
border=NA, col = expert_assessment[2])
points(IPCC_AR5[2], (place.where[2]-0.1), col="black", pch="|")

polygon(x = c(Kopp[1], Kopp[3], Kopp[3], Kopp[1]),
y = c((place.where[3]), (place.where[3]),  (place.where[3] - width), (place.where[3]-width)),
border=NA, col = expert_assessment[3])
points(Kopp[2], (place.where[3]-0.1), col="black", pch="|")

polygon(x = c(pfeffer[1], pfeffer[3], pfeffer[3], pfeffer[1]),
y = c((place.where[4]), (place.where[4]),  (place.where[4] - width), (place.where[4]-width)),
border=NA, col = expert_assessment[4])
points(pfeffer[2], (place.where[4]-0.1), col="black", pch="|")

polygon(x = c(Bamber[1], Bamber[3], Bamber[3], Bamber[1]),
y = c((place.where[5]), (place.where[5]),  (place.where[5] - width), (place.where[5]-width)),
border=NA, col = expert_assessment[5])
points(Bamber[2], (place.where[5]-0.1), col="black", pch="|")

put.fig.letter(label="b.", location="topleft", font=2)

#plot.new()
legend("topleft", c("This study", "Median & 90% C.I. (Little et al. 2013)","Median & 90% C.I. (Church et al. 2013)",
                 "Median & 90% C.I. (Kopp et al. 2014)", "Mean & 90% C.I. (Pfeffer et al. 2008)","Median & 90% C.I. 
(Bamber & Aspinall 2013)"), pch=15, bty="n", cex = 0.75,
       col=c(this.study.color,expert_assessment), y.intersp=c(0.4,0.5,0.5,0.5, 0.5, 0.6))
dev.off()

#------------------------------------- Figure 4 -------------------------------------------
round(mcmc_mean[240100],2)
round(mcmc_std[240100],2)

png(file="Scratch/Figures/Fig4_C_instPaleo.tif", family="Helvetica", units="in", width=text_column_width, height=single_panel_height, pointsize=11, res=300)
par(mfrow=c(1,1), mgp=c(1.5,.5,0),mar=c(4, 4, 2, 1))

plot(date[240000], proj.mcmc.1961_1990[1,240000], type="l", col="powderblue", lwd=2,
     xlab="Year", ylab="AIS Volume loss [SLE m]", 
     xlim=c(-15,106), ylim=c(-0.01,1.4), xaxt="n")#ylim=c(-0.1,1.35), xaxt="n")

# polygon(y_mcmc_90, plus_minusSTD, col=studies.color[2], border=NA)
polygon(y_STD_2100, plus_minusSTD_2100, col=studies.color[2], border=NA)
#polygon(y_STD_2100, plus_minusSTD_2100_NN, col=studies.color[4], border=NA)
#lines(date[1:240100], mcmc_mean_NN[1:240100], col="darkorange2")
#text(85, -0.25, paste(round(mcmc_mean_NN[240100],2), "±", round(mcmc_std_NN[240100],2), sep = ""), col="darkorange2", cex=0.9)

lines(date[1:240100], mcmc_mean[1:240100], col="goldenrod2")
abline(v=100, lty=2)
text(85, 0.2, paste(round(mcmc_mean[240100],2), "±", round(mcmc_std[240100],2), sep = ""), col="goldenrod2", cex=0.9)
# text(85, 0.2, "0.06 ± 0.05", col="goldenrod2", cex=0.75)

# DeConto and Pollard 2016
highPliocene = c((1.05-0.3), 1.05, (1.05+0.3))
lowPliocene = c((0.64-0.49), 0.64, (0.64+0.49))

place.where = c(104, 108)
width = 2.5  # Width of bars

polygon(x = c((place.where[1]), (place.where[1]),  (place.where[1] - width), (place.where[1]-width)),
        y = c(highPliocene[1], highPliocene[3], highPliocene[3], highPliocene[1]), 
        border=NA, col = studies.color[1])
points((place.where[1]-1), highPliocene[2], col="red", pch="_")
text(82.5, 1.1, "1.05 ± 0.30", col="red", cex=0.9)

polygon(x = c((place.where[2]), (place.where[2]),  (place.where[2] - width), (place.where[2]-width)),
        y = c(lowPliocene[1], lowPliocene[3], lowPliocene[3], lowPliocene[1]), 
        border=NA, col = studies.color[3])
points((place.where[2]-1), lowPliocene[2], col="blue", pch="_")
text(82.5, 0.6, "0.64 ± 0.49", col="blue", cex=0.9)

abline(h=0, lty=2, col="gray")

ticks=c(0,20,40,60,80,100)
axis(side=1, at=ticks, labels=expression(2000,2020,2040,2060,2080,2100))

# legend.names = c(expression(paste("±1", sigma, " This study (neglecting cliff 
# instability")), expression(paste("±1", sigma, " DeConto & Pollard 2016 
# (10-20m Pliocene target)")), 
#                  expression(paste("±1", sigma, " DeConto & Pollard 2016 
# (5-15m Pliocene target)")))
legend("topleft", legend=c(expression(paste("±1", sigma)), expression(paste("±1", sigma)),
expression(paste("±1", sigma))), pch=15, col=c(studies.color[2], studies.color[1], studies.color[3]),
bty="n", cex=0.75)#, y.intersp=c(0.5,1.6,1.9))
#(8, 1.42)
legend(-12, 1.45, legend=c("This study neglecting cliff instability",
"DeConto & Pollard 2016 (10-20m Pliocene target)", "DeConto & Pollard 2016 (5-15m Pliocene target)"), pch="",
bty="n", cex=0.75)#, y.intersp=c(0.5,1.1,1.2))

dev.off()

#------------------------------------- Figure 4 No Noise -------------------------------------------
#round(mcmc_mean_NN[240100],2)
#round(mcmc_std_NN[240100],2)

#png(file="Scratch/Figures/Fig4_C_mat_NN.tif", family="Helvetica", units="in", width=3.4, height=3, pointsize=11, res=300)
#par(mfrow=c(1,1), mgp=c(1.5,.5,0),mar=c(4, 4, 2, 1))
#plot(date[240000], proj.mcmc.1961_1990[1,240000], type="l", col="powderblue", lwd=2,
#xlab="Year", ylab="AIS Volume loss [SLE m]",
#xlim=c(-15,106), ylim=c(-0.1,2), xaxt="n")#ylim=c(-0.1,1.35), xaxt="n")

# polygon(y_mcmc_90, plus_minusSTD, col=studies.color[2], border=NA)
#polygon(y_STD_2100, plus_minusSTD_2100_NN, col=studies.color[4], border=NA)
#lines(date[1:240100], mcmc_mean_NN[1:240100], col="darkorange2")
#abline(v=100, lty=2)
#text(85, 0.25, paste(round(mcmc_mean_NN[240100],2), "±", round(mcmc_std_NN[240100],2), sep = ""), col="darkorange2", cex=0.75)
# text(85, 0.2, "0.06 ± 0.05", col="goldenrod2", cex=0.75)

# DeConto and Pollard 2016
#highPliocene = c((1.05-0.3), 1.05, (1.05+0.3))
#lowPliocene = c((0.64-0.49), 0.64, (0.64+0.49))

#place.where = c(104, 108)
#width = 2.5  # Width of bars

#polygon(x = c((place.where[1]), (place.where[1]),  (place.where[1] - width), (place.where[1]-width)),
#y = c(highPliocene[1], highPliocene[3], highPliocene[3], highPliocene[1]),
#border=NA, col = studies.color[1])
#points((place.where[1]-1), highPliocene[2], col="red", pch="_")
#text(82.5, 1.2, "1.05 ± 0.30", col="red", cex=0.75)

#polygon(x = c((place.where[2]), (place.where[2]),  (place.where[2] - width), (place.where[2]-width)),
#y = c(lowPliocene[1], lowPliocene[3], lowPliocene[3], lowPliocene[1]),
#border=NA, col = studies.color[3])
#points((place.where[2]-1), lowPliocene[2], col="blue", pch="_")
#text(82.5, 0.6, "0.64 ± 0.49", col="blue", cex=0.75)

#ticks=c(0,20,40,60,80,100)
#axis(side=1, at=ticks, labels=expression(2000,2020,2040,2060,2080,2100))
#put.fig.letter(label="b.", location="topleft", font=2)
#legend("topleft", legend=c(expression(paste("±1", sigma)), expression(paste("±1", sigma)),
#expression(paste("±1", sigma))), pch=15, col=c(studies.color[2], studies.color[1], studies.color[3]),
#bty="n", cex=0.725, y.intersp=c(0.5,1.6,1.9))

#legend(-7, 2.12, legend=c("This study neglecting cliff
#instability", "DeConto & Pollard 2016
#(10-20m Pliocene target)", "DeConto & Pollard 2016
#(5-15m Pliocene target)"), pch="",
#bty="n", cex=0.725, y.intersp=c(0.5,1.1,1.3))
#dev.off()

#------------------------------------- Marginal pdf plot  -------------------------------------------

#pdf(file="Figures/SuppFigures/Marginalpdf_SF.pdf")
#png(file="Scratch/Figures/SuppFigures/SFig3_C_instPaleo.tif", family="Helvetica", width=text_column_width,
#    height=6, units="in",pointsize=16, res=300)
#png(file="Scratch/Figures/SuppFigures/SFig3_C_instPaleo.tif", family="Helvetica", pointsize=12, res=300)
pdf(file="Scratch/Figures/SuppFigures/Marginalpdf_SF.pdf", family="Helvetica", pointsize=16)
par(mfrow=c(4,4),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

plot(density(DAIS_chains_burnin[,1]), main="", ylab = "PDF", xlab=expression(gamma), yaxt="n", xlim = c(bound.lower[1], bound.upper[1]))
#lines(density(sub_chain[,1]), col="red")
plot(density(DAIS_chains_burnin[,2]), main="", ylab = "", xlab=expression(alpha), yaxt="n", xlim = c(bound.lower[2], bound.upper[2]))
#lines(density(sub_chain[,2]), col="red")
plot(density(DAIS_chains_burnin[,3]), main="", ylab = "", xlab=expression(mu), yaxt="n", xlim = c(bound.lower[3], bound.upper[3]))
#lines(density(sub_chain[,3]), col="red")
plot(density(DAIS_chains_burnin[,4]), main="", ylab = "", xlab=expression(nu), yaxt="n", xlim = c(bound.lower[4], bound.upper[4]))
#lines(density(sub_chain[,4]), col="red")
plot(density(DAIS_chains_burnin[,5]), main="", ylab = "PDF", xlab="P0", yaxt="n", xlim = c(bound.lower[5], bound.upper[5]))
#lines(density(sub_chain[,5]), col="red")
plot(density(DAIS_chains_burnin[,6]), main="", ylab = "", xlab=expression(kappa), yaxt="n", xlim = c(bound.lower[6], bound.upper[6]))
#lines(density(sub_chain[,6]), col="red")
plot(density(DAIS_chains_burnin[,7]), main="", ylab = "", xlab="f0", yaxt="n", xlim = c(bound.lower[7], bound.upper[7]))
#lines(density(sub_chain[,7]), col="red")
plot(density(DAIS_chains_burnin[,8]), main="", ylab = "", xlab="h0", yaxt="n", xlim = c(bound.lower[8], bound.upper[8]))
#lines(density(sub_chain[,8]), col="red")
plot(density(DAIS_chains_burnin[,9]), main="", ylab = "PDF", xlab="c", yaxt="n", xlim = c(bound.lower[9], bound.upper[9]))
#lines(density(sub_chain[,9]), col="red")
plot(density(DAIS_chains_burnin[,10]), main="", ylab = "", xlab="b0", yaxt="n", xlim = c(bound.lower[10], bound.upper[10]))
#lines(density(sub_chain[,10]), col="red")
plot(density(DAIS_chains_burnin[,11]), main="", ylab = "", xlab="slope", yaxt="n", xlim = c(bound.lower[11], bound.upper[11]))
#lines(density(sub_chain[,11]), col="red")
plot(density(DAIS_chains_burnin[,12]), main="", ylab = "", xlab=expression(sigma[P]^{2}), yaxt="n")
#lines(density(sub_chain[,12]), col="red")
plot(density(DAIS_chains_burnin[,13]), main="", ylab = "PDF", xlab=expression(sigma[I]^{2}), yaxt="n", xlim = c(bound.lower[13], bound.upper[13]))
#lines(density(sub_chain[,13]), col="red")

dev.off()

#------------------------------------- Supplementary Figure 3 -------------------------------------------
pairs.image <- function(x) {
    pairs(x, labels = c(expression(gamma),expression(alpha), expression(mu), expression(nu),
    "P0", expression(kappa), "f0", "h0", "c","b0", "slope", expression(sigma[P]^{2}), expression(sigma[I]^{2})), panel=function(x,y) {
    foo <- bin2(cbind(x,y),nbin=c(50,50))
    foo <- ash2(foo,m=c(8,8))
    image(foo,add=T,xlab="",ylab="",col=topo.colors(1000))
  }, upper.panel = NULL)
}
#colnames(d.pos_parameters, do.NULL = FALSE)
#colnames(d.pos_parameters) = c(expression(gamma),expression(alpha), expression(mu), expression(nu),
#                               "Po", expression(kappa), "fo", "ho", "co","bo", "s", expression(sigma^2))
pdf(file="Scratch/Figures/SuppFigures/SFig4_C_instPaleo.pdf", family="Helvetica", pointsize=12)
#png(file="Figures/SuppFig5_dais_mcmc.tif", family="Helvetica",pointsize=12, res=300)
pairs.image(d.pos_parameters)
par(fig = c(.5, 1, .5, 1), new=TRUE)

image.plot(zlim=c(0,5), legend.only=TRUE,col=topo.colors(1000), legend.shrink = 0.585, legend.lab=NULL,
           legend.width = 0.8, horizontal = TRUE, 
           axis.args=list(at=c(0,5),labels=c("Less likely", "More likely")))
dev.off()

#---------------------------------- Comparison to LHS analysis -------------------------------------------
#------------------------ Find the runs that pass through each constraint ---------------------------
source("Scripts/surviveTargetfunc.R")
surv.targ = matrix(c(proj.mcmc.1961_1990[1:subset_length,120000], proj.mcmc.1961_1990[1:subset_length,220000],proj.mcmc.1961_1990[1:subset_length,234000], 
                     proj.mcmc.1961_1990[1:subset_length,240002]), nrow=subset_length, ncol=4)

dais.1992_2011.NN = rep(NA, subset_length)
dais.1992_2011 = rep(NA, subset_length)
for(i in 1:subset_length){
    dais.1992_2011.NN[i] = proj.mcmc.anomaly[i,obs.years[4]] - proj.mcmc.anomaly[i,239992]
    dais.1992_2011[i] = proj.mcmc.1961_1990[i,obs.years[4]] - proj.mcmc.1961_1990[i,239992]
}

surLIG = surviveTarget(windows[1,], surv.targ[,1])
surLGM = surviveTarget(windows[2,], surv.targ[,2])
surMH = surviveTarget(windows[3,], surv.targ[,3])
sur9311trend = surviveTarget(windows[4,], dais.1992_2011)

#Find the runs that pass through the instrumental constraint
total.sur = c(surLIG,surLGM,surMH,sur9311trend)
template <- table(as.vector(total.sur))
all = names(template)[template == max(template)]
sur.all = as.numeric(all)

percent.include = c((subset_length/subset_length)*100, (length(surLIG)/subset_length)*100, 
                    (length(surLGM)/subset_length)*100, (length(surMH)/subset_length)*100, 
                    (length(sur9311trend)/subset_length)*100, (length(sur.all)/subset_length)*100)
print(percent.include)
constraints = c("No constraints","Last integlacial","Last glacial maximum","Mid-Holocene","Instrumental period","All constraints")
table.parts = matrix(c(constraints, percent.include), nrow=6, ncol=2)
write.csv(table.parts, file="Scratch/Random_out/constraint_percent_MCMC_instPaleo.csv")
#---------- NO NOISE
surv.targ = matrix(c(proj.mcmc.anomaly[1:subset_length,120000], proj.mcmc.anomaly[1:subset_length,220000],proj.mcmc.anomaly[1:subset_length,234000], 
                     proj.mcmc.anomaly[1:subset_length,240002]), nrow=subset_length, ncol=4)

surLIG = surviveTarget(windows[1,], surv.targ[,1])
surLGM = surviveTarget(windows[2,], surv.targ[,2])
surMH = surviveTarget(windows[3,], surv.targ[,3])
sur9311trend = surviveTarget(windows[4,], dais.1992_2011.NN)

#Find the runs that pass through the instrumental constraint
total.sur = c(surLIG,surLGM,surMH,sur9311trend)
template <- table(as.vector(total.sur))
all = names(template)[template == max(template)]
sur.all = as.numeric(all)

percent.include = c((subset_length/subset_length)*100, (length(surLIG)/subset_length)*100, 
                    (length(surLGM)/subset_length)*100, (length(surMH)/subset_length)*100, 
                    (length(sur9311trend)/subset_length)*100, (length(sur.all)/subset_length)*100)
print(percent.include)
table.parts.nonoise = matrix(c(constraints, percent.include), nrow=6, ncol=2)
write.csv(table.parts.nonoise, file="Scratch/Random_out/NoNoise_constraint_percent_MCMC_instPaleo.csv")
#------------------------------------- Figure 2 but with no added noise -------------------------------------------
# MCMC hindcasts & projection
png(file="Scratch/Figures/Fig2_NN_C_instPaleo.tif", family="Helvetica", width=text_column_width,
    height=single_panel_height*2, units="in",pointsize=12, res=300)
# pdf(file="Figures/Ruckertetal_dais90MCMC_Fig1.pdf", family="Helvetica", width=6.7, height=8.1, pointsize=12)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

# Last interglacial 240 kyr Bp - 2010 AD
plot(date[1], proj.mcmc.1961_1990[1,1], type="l", col="powderblue", lwd=2,
     xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(date[1], date[235000]), ylim=c(-25,10), xaxt="n")

#polygon(y_mcmc_99_NN, x_mcmc_99_NN, col=this.study.color_99, border=NA)
polygon(y_mcmc_95_NN, x_mcmc_95_NN, col=nonoise.color_95, border=NA)
polygon(y_mcmc_90_NN, x_mcmc_90_NN, col=nonoise.color_90, border=NA)
polygon(y_mcmc_90_NN, plus_minusSTD_NN, col=nonoise.color_68, border=NA)

abline(h=0, lty=3, col="black")
lines(date, mcmc_mean_NN, col="brown")

ticks=c(-200000,-150000,-100000,-50000,0)
axis(side=1, at=ticks, labels=expression(-200,-150,-100,-50,0))
put.fig.letter(label="a.", location="topleft", font=2)

axis(1, at = date[90000]:date[240010], labels = FALSE, col=boxpalette[1], lwd=2)

# Last interglacial 150 kyr Bp - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))

plot(date[90000], proj.mcmc.1961_1990[1,90000], type="l", col="powderblue", 
     lwd=2,xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]",
     xlim=c(date[90000], date[235000]), ylim=c(-5,6), xaxt="n")

#polygon(y_mcmc_99_NN, x_mcmc_99_NN, col=this.study.color_99, border=NA)
polygon(y_mcmc_95_NN, x_mcmc_95_NN, col=nonoise.color_95, border=NA)
polygon(y_mcmc_90_NN, x_mcmc_90_NN, col=nonoise.color_90, border=NA)
polygon(y_mcmc_90_NN, plus_minusSTD_NN, col=nonoise.color_68, border=NA)

lines(date, mcmc_mean_NN, col="brown")
lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="black", lwd=2)
lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="black", lwd=2)
lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="black", lwd=2)
points(last.interglacial[2],median(windows[1,]), col="black", pch=8, cex=0.585)
ticks=c(-150000,-125000,-100000, -75000, -50000, -25000, 0)
axis(side=1, at=ticks, labels=expression(-150, -125, -100,-75,-50,-25, 0))
put.fig.letter(label="b.", location="topleft", font=2)

box(col = boxpalette[1])
axis(1, at = date[kyrbp_25]:date[240010], labels = FALSE, col=boxpalette[1], lwd=2)

# Last glacial maximum 25 kyr Bp - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

plot(date[kyrbp_25], proj.mcmc.1961_1990[1,kyrbp_25], type="l", col="powderblue", 
     lwd=2,xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]",
     xlim=c(date[kyrbp_25], date[238000]), ylim=c(-21,0.5), xaxt="n")

#polygon(y_mcmc_99_NN, x_mcmc_99_NN, col=this.study.color_99, border=NA)
polygon(y_mcmc_95_NN, x_mcmc_95_NN, col=nonoise.color_95, border=NA)
polygon(y_mcmc_90_NN, x_mcmc_90_NN, col=nonoise.color_90, border=NA)
polygon(y_mcmc_90_NN, plus_minusSTD_NN, col=nonoise.color_68, border=NA)

lines(date, mcmc_mean_NN, col="brown")
lines(c(last.glacialmax[2],last.glacialmax[2]), c(windows[2,1],windows[2,2]), col="black", lwd=2)
lines(last.glacialmax, c(windows[2,1],windows[2,1],windows[2,1]), col="black", lwd=2)
lines(last.glacialmax, c(windows[2,2],windows[2,2],windows[2,2]), col="black", lwd=2)
points(last.glacialmax[2],median(windows[2,]), col="black", pch=8, cex=0.585)
ticks=c(-25000,-20000,-15000,-10000,-5000,0)
axis(side=1, at=ticks, labels=expression(-25,-20,-15,-10,-5,0))
put.fig.letter(label="c.", location="topleft", font=2)

box(col = boxpalette[1])
axis(1, at = date[kyrbp_6]:date[240010], labels = FALSE, col=boxpalette[2], lwd=2)

# Mid-Holocene6 kyr Bp - 2010 AD 
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))

plot(date[kyrbp_6], proj.mcmc.1961_1990[1,kyrbp_6], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(date[kyrbp_6], date[240010]),ylim=c(-5,3)) #ylim=c(-5,0.5))

#polygon(y_mcmc_99_NN, x_mcmc_99_NN, col=this.study.color_99, border=NA)
polygon(y_mcmc_95_NN, x_mcmc_95_NN, col=nonoise.color_95, border=NA)
polygon(y_mcmc_90_NN, x_mcmc_90_NN, col=nonoise.color_90, border=NA)
polygon(y_mcmc_90_NN, plus_minusSTD_NN, col=nonoise.color_68, border=NA)

lines(date, mcmc_mean_NN, col="brown")
lines(c(holocene[2],holocene[2]), c(windows[3,1],windows[3,2]), col="black", lwd=2)
lines(holocene, c(windows[3,1],windows[3,1],windows[3,1]), col="black", lwd=2)
lines(holocene, c(windows[3,2],windows[3,2],windows[3,2]), col="black", lwd=2)
points(holocene[2],median(windows[3,]), col="black", pch=8, cex=0.585)
put.fig.letter(label="d.", location="topleft", font=2)

box(col = boxpalette[2])
axis(1, at = date[AD_1880]:date[240010], labels = FALSE, col=boxpalette[3], lwd=2)

# Present day 1880 - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))
plot(date[AD_1880], proj.mcmc.1961_1990[1,AD_1880], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(date[AD_1880], date[240010]), xaxt="n", ylim=c(-0.04,0.03)) #ylim=c(-0.04,0.02))

#polygon(y_mcmc_99_NN, x_mcmc_99_NN, col=this.study.color_99, border=NA)
polygon(y_mcmc_95_NN, x_mcmc_95_NN, col=nonoise.color_95, border=NA)
polygon(y_mcmc_90_NN, x_mcmc_90_NN, col=nonoise.color_90, border=NA)
polygon(y_mcmc_90_NN, plus_minusSTD_NN, col=nonoise.color_68, border=NA)

lines(date, mcmc_mean_NN, col="brown")
#err_positive = windows[4,2]
#err_negative = windows[4,1]
#present = 2
#arrows(present, err_positive, present, err_negative, length=0, lwd=2, col="black")
#points(present,median(windows[4,]), col="black", pch=8, cex=0.585)
ticks=c(-120,-100,-50,0)
axis(side=1, at=ticks, labels=expression(1880,1900,1950,2000))
put.fig.letter(label="e.", location="topleft", font=2)

box(col = boxpalette[3])
axis(1, at = date[240000]:date[240020], labels = FALSE, col=boxpalette[4], lwd=2)

# Future 2010 - 2300 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
plot(date[240000], proj.mcmc.1961_1990[1,240000], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(-10,285), xaxt="n", ylim=c(-0.2,7)) #ylim=c(-0.5,4))

#polygon(y_mcmc_99_NN, x_mcmc_99_NN, col=this.study.color_99, border=NA)
polygon(y_mcmc_95_NN, x_mcmc_95_NN, col=nonoise.color_95, border=NA)
polygon(y_mcmc_90_NN, x_mcmc_90_NN, col=nonoise.color_90, border=NA)
polygon(y_mcmc_90_NN, plus_minusSTD_NN, col=nonoise.color_68, border=NA)

lines(date, mcmc_mean, col="brown")
#err_positive = windows[4,2]
#err_negative = windows[4,1]
#arrows(present, err_positive, present, err_negative, length=0, lwd=2, col="black")
ticks=c(0,50,100,150,200,250,300)
axis(side=1, at=ticks, labels=expression(2000,2050,2100,2150,2200,2250,2300))
put.fig.letter(label="f.", location="topleft", font=2)

box(col = boxpalette[4])

legend("topleft", expression(paste("±2 ", sigma, " Data constraint")), pch="I", bty="n", col="black")
legend("topleft", c("", "Mean hindcast & projection", "95% credible interval", "90% credible interval",
expression(paste("±1 ", sigma, " credible interval"))),
lty=c(NA,1,NA,NA,NA), pch=c(NA,NA,15,15,15), bty="n",
col=c("white", "brown", nonoise.color_95, nonoise.color_90, nonoise.color_68))

#legend("topleft", expression(paste("±2 ", sigma, " Data constraint")), pch="I", bty="n", col="black")
#legend("topleft", c("", "Mean hindcast & projection", "95% credible interval", "90% credible interval",
#"68% credible interval"), lty=c(NA,1,NA,NA,NA), pch=c(NA,NA,15,15,15), bty="n",
#col=c("white", "brown", this.study.color_99, this.study.color_95, this.study.color))
#legend("topleft", c("", "Mean hindcast & projection", "95% credible interval", "90% credible interval",
#expression(paste("±1 ", sigma, " credible interval"))),
#lty=c(NA,1,NA,NA,NA), pch=c(NA,NA,15,15,15), bty="n",
#col=c("white", "brown", this.study.color_99, this.study.color_95, this.study.color))

dev.off()

#------------------------------------- Supp Figure X -------------------------------------------
# Function to find SLE values in certain years 'fn.prob.proj'
# Calculate the pdf, cdf, and sf of AIS melt estimates in:
NN_mcmc.prob_proj <- fn.prob.proj(proj.mcmc.anomaly, year.pcs, subset_length, un.constr=T)
NN_LIG.mcmc <- plot.sf(NN_mcmc.prob_proj[,1], make.plot=F) # 120,000 BP (Last interglacial)
NN_LGM.mcmc <- plot.sf(NN_mcmc.prob_proj[,2], make.plot=F)
NN_MH.mcmc <- plot.sf(NN_mcmc.prob_proj[,3], make.plot=F)
#NN_2002.mcmc <- plot.sf(NN_mcmc.prob_proj[,4], make.plot=F)
NN_2002.mcmc <- plot.sf(dais.1992_2011.NN, make.plot=F)
sf.2002.mcmc <- plot.sf(dais.1992_2011, make.plot=F)

#pdf(file="Figures/Ruckertetal_daisLIG2100pdf_Fig2aa.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
png(file="Scratch/Figures/FigXa_instPaleo.tif", family="Helvetica", width=3.5,
height=8, units="in",pointsize=12, res=300)
par(mfrow=c(3,1),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1.5)) # set figure dimensions

plot(LGM.sf.mcmc$pdf, main="",lwd=1, col=this.study.color, xlab="Hindcasted AIS Volume loss", sub="during Last glacial maximum [SLE m]",
ylab="Probability Density", ylim=c(-0.03,0.19),xlim=c(min(LGM.sf.mcmc$pdf$x),8), yaxt="n")

lines(LGM.sf.mcmc$pdf, col=this.study.color, lwd=1.5)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,2], probs, width.size = 0.02, where.at = -0.02, tick.length = 0.004, line.width = 1.5, color = this.study.color)

abline(h=0.13, lty=2)
text(-2, 0.13+0.02, cex=0.95, "Last glacial maximum
Data Constraint")
text(-2, 0.13-0.02, cex=0.95, "Model Inversion
This study")

plotrange(windows[2,2], (windows[2,1] + windows[2,2])/2, windows[2,1], year=F, height=0.18, color="black")
put.fig.letter(label="a.", location="topleft", font=2)

#---
plot(MH.sf.mcmc$pdf, main="",lwd=1, col=this.study.color, xlab="Hindcasted AIS Volume loss", sub="during the Mid-Holocene [SLE m]",
ylab="Probability Density", ylim=c(-0.07,0.52),xlim=c(min(MH.sf.mcmc$pdf$x),6), yaxt="n")

lines(MH.sf.mcmc$pdf, col=this.study.color, lwd=1.5)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,3], probs, width.size = 0.07, where.at = -0.05, tick.length = 0.01, line.width = 1.5, color = this.study.color)

abline(h=0.4, lty=2)
text(3, 0.4+0.06, cex=0.95, "Mid-Holocene
Data Constraint")
text(3, 0.4-0.06, cex=0.95, "Model Inversion
This study")

plotrange(windows[3,2], (windows[3,1] + windows[3,2])/2, windows[3,1], year=F, height=0.5, color="black")
put.fig.letter(label="b.", location="topleft", font=2)

#---
plot(sf.2002.mcmc$pdf, main="",lwd=1, col=this.study.color, xlab="Hindcasted AIS Volume loss", sub="during the Instrumental period [SLE m]",
ylab="Probability Density", ylim=c(-0.2,2),xlim=c(min(sf.2002.mcmc$pdf$x),2), yaxt="n")

lines(sf.2002.mcmc$pdf, col=this.study.color, lwd=1.5)

probs = c(0.05, 0.95)
#add.hor.box(mcmc.prob_proj[,4], probs, width.size = 0.1, where.at = -0.06, tick.length = 0.01, line.width = 1.5, color = this.study.color)

add.hor.box(dais.1992_2011, probs, width.size = 0.2, where.at = -0.1, tick.length = 0.01, line.width = 1.5, color = this.study.color)

abline(h=1, lty=2)
text(1, 1+0.12, cex=0.95, "Instrumental period
Data Constraint
(relative to 1992)")#0.62
text(1, 1-0.1, cex=0.95, "Model Inversion
This study")

plotrange(windows[4,2], (windows[4,1] + windows[4,2])/2, windows[4,1], year=F, height=0.7, color="black")
put.fig.letter(label="c.", location="topleft", font=2)

dev.off()

#-------------------------------
png(file="Scratch/Figures/FigXb_instPaleo.tif", family="Helvetica", width=5.2,
height=4, units="in",pointsize=10, res=300)
par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1.5)) # set figure dimensions

plot(NN_LIG.mcmc$pdf, main="",lwd=1, col=nonoise.color_68, xlab="Hindcasted AIS Volume loss", sub="during Last interglacial [SLE m]",
ylab="Probability Density", ylim=c(-0.07,0.65),xlim=c(min(NN_LIG.mcmc$pdf$x),17) , yaxt="n")

lines(NN_LIG.mcmc$pdf, col=nonoise.color_68, lwd=1.5)

probs = c(0.05, 0.95)
add.hor.box(NN_mcmc.prob_proj[,1], probs, width.size = 0.07, where.at = -0.05, tick.length = 0.01, line.width = 1.5, color = nonoise.color_68)

abline(h=0.4, lty=2)
text(12, 0.4+0.1, cex=0.8, "Last Interglacial period
Data Constraint")
text(12, 0.4-0.15, cex=0.8, "Model Inversion
This study without
added noise")

plotrange(windows[1,2], (windows[1,1] + windows[1,2])/2, windows[1,1], year=F, height=0.55, color="black")
put.fig.letter(label="a.", location="topleft", font=2)

#---
plot(NN_LGM.mcmc$pdf, main="",lwd=1, col=nonoise.color_68, xlab="Hindcasted AIS Volume loss", sub="during Last glacial maximum [SLE m]",
ylab="Probability Density", ylim=c(-0.03,0.22),xlim=c(min(NN_LGM.mcmc$pdf$x),10), yaxt="n")

lines(NN_LGM.mcmc$pdf, col=nonoise.color_68, lwd=1.5)

probs = c(0.05, 0.95)
add.hor.box(NN_mcmc.prob_proj[,2], probs, width.size = 0.02, where.at = -0.02, tick.length = 0.004, line.width = 1.5, color = nonoise.color_68)

abline(h=0.13, lty=2)
text(-1, 0.13+0.03, cex=0.85, "Last glacial maximum
Data Constraint")
text(-1, 0.13-0.04, cex=0.85, "Model Inversion
This study without
added noise")

plotrange(windows[2,2], (windows[2,1] + windows[2,2])/2, windows[2,1], year=F, height=0.21, color="black")
put.fig.letter(label="b.", location="topleft", font=2)

#---
plot(NN_MH.mcmc$pdf, main="",lwd=1, col=nonoise.color_68, xlab="Hindcasted AIS Volume loss", sub="during the Mid-Holocene [SLE m]",
ylab="Probability Density", ylim=c(-0.07,0.57),xlim=c(min(NN_MH.mcmc$pdf$x),7), yaxt="n")

lines(NN_MH.mcmc$pdf, col=nonoise.color_68, lwd=1.5)

probs = c(0.05, 0.95)
add.hor.box(NN_mcmc.prob_proj[,3], probs, width.size = 0.07, where.at = -0.05, tick.length = 0.01, line.width = 1.5, color = nonoise.color_68)

abline(h=0.45, lty=2)
text(3, 0.45+0.08, cex=0.85, "Mid-Holocene
Data Constraint")
text(3, 0.45-0.1, cex=0.85, "Model Inversion
This study without
added noise")

plotrange(windows[3,2], (windows[3,1] + windows[3,2])/2, windows[3,1], year=F, height=0.55, color="black")
put.fig.letter(label="c.", location="topleft", font=2)

#---
plot(NN_2002.mcmc$pdf, main="",lwd=1, col=nonoise.color_68, xlab="Hindcasted AIS Volume loss", sub="during the Instrumental period [SLE m]",
ylab="Probability Density", ylim=c(-20,440),xlim=c(min(NN_2002.mcmc$pdf$x),0.02), yaxt="n")

lines(NN_2002.mcmc$pdf, col=nonoise.color_68, lwd=1.5)

probs = c(0.05, 0.95)
#add.hor.box(NN_mcmc.prob_proj[,4], probs, width.size = 15, where.at = -10, tick.length = 3, line.width = 1.5, color = nonoise.color_68)

add.hor.box(dais.1992_2011.NN, probs, width.size = 30, where.at = -20, tick.length = 5, line.width = 1.5, color = nonoise.color_68)

abline(h=320, lty=2)
text(0.015, 320+70, cex=0.85, "Instrumental period
Data Constraint
(relative to 1992)")
text(0.015, 320-70, cex=0.85, "Model Inversion
This study without
added noise")

plotrange(windows[4,2], (windows[4,1] + windows[4,2])/2, windows[4,1], year=F, height=350, color="black")
put.fig.letter(label="d.", location="topleft", font=2)

dev.off()


###################################### END ############################################

png(file="Scratch/Figures/Fig3_C_newinstPaleo2.tif", family="Helvetica", width=text_column_width,
height=single_panel_height*2, units="in",pointsize=12, res=300)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

plot(LIG.sf.mcmc$pdf, main="",lwd=1.5, col=this.study.color, xlab="Hindcasted AIS Volume loss", sub="during Last interglacial [SLE m]",
ylab="Probability Density", ylim=c(-0.07,0.42), yaxt="n",xlim=c(min(LIG.sf.mcmc$pdf$x),17))

lines(LIG.sf.mcmc$pdf, col=this.study.color, lwd=1.5)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,1], probs, width.size = 0.07, where.at = -0.05, tick.length = 0.01, line.width = 1.5, color = this.study.color)

abline(h=0.3, lty=2)
abline(v=0, lty=2, col="gray")
text(11, 0.3+0.06, cex=0.95, "Last Interglacial period
Data Constraint")
text(11, 0.3-0.06, cex=0.95, "Model Inversion
This study")

plotrange(windows[1,2], (windows[1,1] + windows[1,2])/2, windows[1,1], year=F, height=0.4, color="black")
axis(3, labels = FALSE)
put.fig.letter(label="a.", location="topleft", font=2)

#---
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
plot(LGM.sf.mcmc$pdf, main="",lwd=1, col=this.study.color, xlab="Hindcasted AIS Volume loss", sub="during Last glacial maximum [SLE m]",
ylab="Probability Density", ylim=c(-0.03,0.19),xlim=c(min(LGM.sf.mcmc$pdf$x),8), yaxt="n")

lines(LGM.sf.mcmc$pdf, col=this.study.color, lwd=1.5)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,2], probs, width.size = 0.02, where.at = -0.02, tick.length = 0.004, line.width = 1.5, color = this.study.color)

abline(h=0.13, lty=2)
abline(v=0, lty=2, col="gray")

boxed.labels(-2, 0.13+0.02, "Last glacial maximum
Data Constraint", bg="white", border=NA, cex=0.6, xpad=1)
#text(-2, 0.13+0.02, cex=0.9, "Last glacial maximum
#Data Constraint")
#text(-2, 0.13-0.02, cex=0.9, "Model Inversion
#This study")
boxed.labels(-2, 0.13-0.02, "Model Inversion
This study", bg="white", border=NA, cex=0.6, xpad=1)

plotrange(windows[2,2], (windows[2,1] + windows[2,2])/2, windows[2,1], year=F, height=0.18, color="black")
axis(3, labels = FALSE)
put.fig.letter(label="b.", location="topleft", font=2)

#---
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))
plot(MH.sf.mcmc$pdf, main="",lwd=1, col=this.study.color, xlab="Hindcasted AIS Volume loss", sub="during the Mid-Holocene [SLE m]",
ylab="Probability Density", ylim=c(-0.07,0.52),xlim=c(min(MH.sf.mcmc$pdf$x),6), yaxt="n")

lines(MH.sf.mcmc$pdf, col=this.study.color, lwd=1.5)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,3], probs, width.size = 0.07, where.at = -0.05, tick.length = 0.01, line.width = 1.5, color = this.study.color)

abline(h=0.4, lty=2)
abline(v=0, lty=2, col="gray")
boxed.labels(3, 0.4+0.06, "Mid-Holocene
Data Constraint", bg="white", border=NA, cex=0.6, xpad=1)
#text(3, 0.4+0.06, cex=0.9, "Mid-Holocene
#Data Constraint")
#text(3, 0.4-0.06, cex=0.9, "Model Inversion
#This study")
boxed.labels(3, 0.4-0.06, "Model Inversion
This study", bg="white", border=NA, cex=0.6, xpad=1)

plotrange(windows[3,2], (windows[3,1] + windows[3,2])/2, windows[3,1], year=F, height=0.5, color="black")
axis(3, labels = FALSE)
put.fig.letter(label="c.", location="topleft", font=2)

#---
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
plot(sf.2002.mcmc$pdf, main="",lwd=1, col=this.study.color, xlab="Hindcasted AIS Volume loss", sub="during the Instrumental period [SLE m]",
ylab="Probability Density", ylim=c(-4,50),xlim=c(min(sf.2002.mcmc$pdf$x),0.1), yaxt="n")

lines(sf.2002.mcmc$pdf, col=this.study.color, lwd=1.5)

probs = c(0.05, 0.95)
#add.hor.box(mcmc.prob_proj[,4], probs, width.size = 0.1, where.at = -0.06, tick.length = 0.01, line.width = 1.5, color = this.study.color)

add.hor.box(dais.1992_2011, probs, width.size = 5, where.at = -3, tick.length = 1, line.width = 1.5, color = this.study.color)

abline(h=40, lty=2)
abline(v=0, lty=2, col="gray")
text(0.05, 40+6, cex=0.9, "Instrumental period
Data Constraint
(relative to 1992)")#0.62
text(0.05, 40-4, cex=0.9, "Model Inversion
This study")

plotrange(windows[4,2], (windows[4,1] + windows[4,2])/2, windows[4,1], year=F, height=45, color="black")
axis(3, labels = FALSE)
put.fig.letter(label="d.", location="topleft", font=2)

#---
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))
plot(sf.2100.mcmc$pdf, main="", col=this.study.color,lwd=1.5,
xlab="Projected AIS Volume loss", sub="in 2100 [SLE m]",
ylab="Probability Density",
xlim=c(-0.1,1.1),
ylim=c(-0.6,11.5), yaxt="n") #ylim=c(-7,53), yaxt="n") #2.4

lines(sf.2100.mcmc$pdf, col=this.study.color, lwd=1.5)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,6], probs, width.size = 1, where.at = -.4, tick.length = 0.12, line.width = 1.5, color = this.study.color)

rect(-0.3, 10.1, 1.3, 12.0, col = "gray93", border=NA)

abline(h=7.5, lty=2)
abline(v=0, lty=2, col="gray")
text(0.65, 7.5+0.4, cex=0.9, "Expert Assessments")
text(0.65, 7.5-0.8, cex=0.9, "Model Inversion
This study")

place.where = c(8.3, 8.8, 9.3, 9.8, 10.8, 11.3, 11.8)
width = 0.4  # Width of bars

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
put.fig.letter(label="e.", location="topleft", font=2)
box()

#---
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
plot(1, 1, col="white", xlab="", ylab="", xlim=c(-10,290),ylim=c(-10,20), xaxt="n", yaxt="n")
legend("center", c("This study", "Median & 90% C.I. (Little et al. 2013)","Median & 90% C.I. (Kopp et al. 2014)",
"Median & 90% C.I. (Ritz et al. 2015)", "'Low' & 'High' estimate\n(Golledge et al. 2015)",
"Median & 90% C.I. (Church et al. 2013)", "'Low2', 'Low1', & 'High' estimate\n(Pfeffer et al. 2008)","Median & 90% C.I.\n(Bamber & Aspinall 2013)"), pch=15, bty="n", cex = 0.85, y.intersp=1.3,
col=c(this.study.color,model_assessment_colors, nonmodel_assessment_colors))#, y.intersp=c(0.4,0.5,0.5,0.5, 0.5, 0.6))

dev.off()

