##################################################################################################
#
#  -file = "LHS_plots.R"   Code written September 2014 edited March 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in the LHS workspace
#
#   -NOTE: The graphs will be saved as tif files in the current working directory.
#
###################################################################################################

#load("Scratch/Workspace/DAIS_precalibration_LHS_C.RData") # Load in the saved workspace from LHS precalibration
load("Scratch/Workspace/DAIS_precalibration_LHS_relative_2.RData")

#install.packages('ash')
library(ash)
#install.packages('fields')
library(fields)
#install.packages('RColorBrewer')
library(RColorBrewer)

mypalette <- brewer.pal(9,"YlGnBu")
boxpalette <- brewer.pal(11,"Spectral")
expert_assessment <- brewer.pal(5,"Pastel2")

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
#Set up previous study ranges for the year 2100 using RCP 8.5
#the previous ranges are in cm so divide by 100 to get meters
pfeffer = c(12.8, 14.6, 61.9)/100 #Low2, low1, and high estimate in Pfeffer et al. 2008
Bamber = c(-2, 13, 83)/100 # 5%, Median, 95% estimates in Bamber and Aspinall 2013
IPCC_AR5 = c(-6, 4, 12)/100 # 5%, Median, 95% estimates in IPCC AR5
Little = c(-8, 2.4, 13.3)/100 # 5%, Median, 95% estimates in Little et al. 2013
# Little_unweighted = c(-14.7, -1.1, 12.6)/100 # 5%, Median, 95% estimates in Little et al. 2013
Kopp = c(-11, 4, 33)/100 # 5%, Median, 95% estimates in Kopp et al. 2014
present = 2

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
#------------------------------------- Supplementary Figure 3 -------------------------------------------

png(file="Scratch/Figures/SuppFigures/S2_Fig_inst2.tif", family="Helvetica", width=text_column_width,
height=single_panel_height*2, units="in",pointsize=12, res=300)
par(mfrow=c(2,1),mgp=c(1.5,.5,0), mar=c(3.5,4,1,2)) # set figure dimensions

plot(un.sflig$pdf, main="",lwd=3, col="gray91", xlab="Projected AIS volume loss", sub="during Last interglacial [SLE m]",
ylab="Probability Density",xlim=c(-10,25), ylim=c(-0.45,0.6), yaxt="n")
#ylab="Probability Density",xlim=c(min(un.sflig$pdf$x),20), ylim=c(-0.45,0.6), yaxt="n")

lines(LIG.sflig$pdf, col=mypalette[3], lwd=2)
lines(LGM.sflig$pdf, col=mypalette[5], lwd=2)
lines(MH.sflig$pdf, col=mypalette[1], lwd=2)
lines(present.sflig$pdf, col=mypalette[7], lwd=2)
lines(all.sflig$pdf, col=mypalette[9], lwd=2)

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
text(14, 0.4+0.1, cex=0.75, "Last Interglacial period
Data Constraint")
text(14, 0.4-0.1, cex=0.75, "Model Inversion
This study")

plotrange(windows[1,2], (windows[1,1] + windows[1,2])/2, windows[1,1], year=F, height=0.55, color="black")

put.fig.letter(label="a.", location="topleft", font=2)

par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 1, 2))
plot(un.sf2100$pdf, main="",lwd=3, col="gray91", xlab="Projected AIS volume loss", sub= "in 2100 [SLE m]",
ylab="Probability Density",
xlim=c(-1.5,3.5), ylim=c(-5.5,14.5))#, yaxt="n") #ylim=c(-18,78)

lines(LIG.sf2100$pdf, col=mypalette[3], lwd=2)
lines(LGM.sf2100$pdf, col=mypalette[5], lwd=2)
lines(MH.sf2100$pdf, col=mypalette[1], lwd=2)
lines(present.sf2100$pdf, col=mypalette[7], lwd=2)
lines(all.sf2100$pdf, col=mypalette[9], lwd=2)

probs = c(0.05, 0.95)
#where = c(-0.2, -0.5, -1, -1.5, -2, -2.5)
add.hor.box(all.prob_proj[,6], probs, width.size = 1.2, where.at = -0.5, tick.length = 0.1, line.width = 2, color = mypalette[9])
#points(dais.pre.cali[sur.all,240100], -0.2, pch=15, col= mypalette[9], cex=0.75)
add.hor.box(present.prob_proj[,6], probs, width.size = 1.2, where.at = -1.5, tick.length = 0.1, line.width = 2, color = mypalette[7])
add.hor.box(MHprob_proj[,6], probs, width.size = 1.2, where.at = -2.5, tick.length = 0.1, line.width = 2, color = mypalette[1])
add.hor.box(LGMprob_proj[,6], probs, width.size = 1.2, where.at = -3.5, tick.length = 0.1, line.width = 2, color = mypalette[5])
add.hor.box(LIGprob_proj[,6], probs, width.size = 1.2, where.at = -4.5, tick.length = 0.1, line.width = 2, color = mypalette[3])
add.hor.box(unprob_proj[,6], probs, width.size = 1.2, where.at = -5.5, tick.length = 0.1, line.width = 2, color = "gray91")

abline(h=11, lty=2)
text(-1, 11+0.5, cex=0.65, "Expert Assessments")
text(-1, 11-1, cex=0.65, "Model Inversion
This study")

place.where = c(11.6,12.25,13,13.75,14.5)
#place.where = c(2.6,3.25,4,4.75,5.5)
width = 0.45  # Width of bars

polygon(x = c(Little[1], Little[3], Little[3], Little[1]),
y = c((place.where[1]), (place.where[1]),  (place.where[1] - width), (place.where[1]-width)),
border=NA, col = expert_assessment[1])
points(Little[2], (place.where[1]-0.2), col="black", pch="|", cex=0.5)

polygon(x = c(IPCC_AR5[1], IPCC_AR5[3], IPCC_AR5[3], IPCC_AR5[1]),
y = c((place.where[2]), (place.where[2]),  (place.where[2] - width), (place.where[2]-width)),
border=NA, col = expert_assessment[2])
points(IPCC_AR5[2], (place.where[2]-0.2), col="black", pch="|", cex=0.5)

polygon(x = c(Kopp[1], Kopp[3], Kopp[3], Kopp[1]),
y = c((place.where[3]), (place.where[3]),  (place.where[3] - width), (place.where[3]-width)),
border=NA, col = expert_assessment[3])
points(Kopp[2], (place.where[3]-0.2), col="black", pch="|", cex=0.5)

polygon(x = c(pfeffer[1], pfeffer[3], pfeffer[3], pfeffer[1]),
y = c((place.where[4]), (place.where[4]),  (place.where[4] - width), (place.where[4]-width)),
border=NA, col = expert_assessment[4])
points(pfeffer[2], (place.where[4]-0.2), col="black", pch="|", cex=0.5)

polygon(x = c(Bamber[1], Bamber[3], Bamber[3], Bamber[1]),
y = c((place.where[5]), (place.where[5]),  (place.where[5] - width), (place.where[5]-width)),
border=NA, col = expert_assessment[5])
points(Bamber[2], (place.where[5]-0.2), col="black", pch="|", cex=0.5)

put.fig.letter(label="b.", location="topleft", font=2)

legend.names = c(paste("Unconstrained L.H. fits (", sample_length, ")", sep=""),
paste("LIG constraint fits (", length(surLIG), ")", sep=""),
paste("LGM constraint fits (", length(surLGM), ")", sep=""),
paste("MH constraint fits (", length(surMH), ")", sep=""),
paste("Instrumental constraint fits (", length(sur9311trend), ")", sep=""),
#paste("All L.H. constrained fits (0)", sep=""),# length(sur.all), ")", sep=""),
paste("All L.H. constrained fits (", length(sur.all), ")", sep=""),
"Median & 90% C.I. (Little et al. 2013)","Median & 90% C.I. (Church et al. 2013)",
"Median & 90% C.I. (Kopp et al. 2014)", "Mean & 90% C.I. (Pfeffer et al. 2008)",
"Median & 90% C.I. (Bamber & Aspinall 2013)")

legend("topright", legend.names, cex=0.575, pch=15, #bty="n",
col=c("gray91", mypalette[3], mypalette[5], mypalette[1], mypalette[7],mypalette[9],expert_assessment))#,
#y.intersp=c(0.5,0.6,0.6,0.6, 0.6, 0.6,0.6,0.6,0.6,0.6,0.6))

dev.off()

#------------------------------------- Supplementary Figure 2 -------------------------------------------
# # LHS hindcasts & projection
png(file="Scratch/Figures/SuppFigures/S1_Fig_inst2.tif", family="Helvetica", width=text_column_width, height=single_panel_height*2, units="in",pointsize=12, res=300)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

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
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", ylim=c(-0.4,0.4), xaxt="n") #(-0.10,0.05)
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
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", xlim=c(-10,290),ylim=c(-0.4,20), xaxt="n") #ylim=c(-0.5,5)
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

