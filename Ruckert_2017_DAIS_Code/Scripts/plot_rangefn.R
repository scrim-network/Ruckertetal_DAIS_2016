#######################################################################
## -file = plot_rangefn.R
## -Input: For ploting a range [plotrange]:
##         The lower bound of the range [low.number]
##         A value between the lower and upper bound [median]
##         The upper bound of the range [high.number]
##         Plot the range vertically; the range is associate with y-axis [year; T or F]
##         Plot the range horizontally; the range is associate with x-axis  [height; T or F]
##         Color of lines and points [color]
##
## -Input: For ploting a blox plot [add.hor.box]:
##         Vector of values [data.numbers]
##         Vector of two values specifying other probabilities to be estimated and ploted [probabilities]
##         Width of the box [width.size]
##         Location where the boxplot should be drawn [where.at]
##         Width of line segments showing the specified probability [line.width]
##         Length of line segments showing the specified probability [tick.length]
##         Color of lines and points [color]
##
## -Output: Adds a range or box and whisker to an existing plot
#############################################################
## Copyright 2015 Kelsey Ruckert
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
##
## -Author: Kelsey Ruckert (klr324@psu.edu)
## -June 10 2015
#######################################################################
# Function to plot a range on an existing plot
plotrange <- function(low.number, median, high.number, year=T, height=F, color){
  
  range <- c(low.number, median, high.number)
  
  if(year){
    plot.position <- c(year, year, year)
    
    points(plot.position[1:2], c(range[1], range[3]), col=color, cex=1, pch="-")
    points(plot.position[1], range[2], col=color, cex=1, pch=8)
    segments(plot.position[1], range[1], plot.position[3], range[3], lwd=2, col=color)
  }
  
  if(height){
    plot.position <- c(height, height, height)
    
    points(c(range[1], range[3]), plot.position[1:2], col=color, cex=1, pch="|")
    points(range[2], plot.position[1], col=color, cex=1, pch=8)
    segments(range[1], plot.position[1], range[3], plot.position[3], lwd=2, col=color)
  }
}

# Function to plot a blox plot on an existing plot with adding lines showing select probabilities
add.hor.box <- function(data.numbers, probabilities, width.size, where.at, tick.length, line.width, color){
    
  quants <- quantile(data.numbers, probabilities)
  boxplot(data.numbers, add=TRUE, horizontal=TRUE, axes=FALSE, outline=FALSE, col=color, boxwex=width.size, at=where.at)
  segments(x0 = quants, y0 = rep(where.at - tick.length, 2),
           x1 = quants, y1 = rep(where.at + tick.length, 2), col = c("blue", "blue"), lwd = line.width)
}