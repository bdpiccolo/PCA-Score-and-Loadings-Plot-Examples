
# Plotting PCA Scores and Loadings in Scatterplots #

## Required packages

# ```{r PCA_packages}
library(ChemometricsWithRData)
library(plyr)
library(car)
library(maptools)
library(rgeos)

# if these pacakges are not installed in your library, then use the following code to install the packages

#if (length(setdiff(c("ChemometricsWithRData", "plyr", "car","maptools","rgeos"), rownames(installed.packages()))) > 0) {
#	install.packages(setdiff(c("ChemometricsWithRData", "plyr", "car","maptools","rgeos"), rownames(installed.packages())))
#}
# ```

## Data ##  

# We will use the wines data set from the ChemometricsWithRData package

# ```{r PCA_data}

# load data
data(wines)

# Data
vino <- wines

# Wine classes
vint <- vintages

# Look at data
head(vino)

# ```
##  Transform and scale data ##
Log transformation and autoscaling to set all variables on same scale.
# ```{r PCA_pp}
# log transform, pseudo scale - not entirely necessary
vino_log <- log(vino)

# Auto scale
vino_log_scale <- scale(vino_log)
head(vino_log_scale)
# ```
## PCA ##
# Use the singular value decomposition algorithm - prcomp() function
# ```{r PCA_PCA}

vino_PCA <- prcomp(vino_log_scale, center=FALSE)
summary(vino_PCA)
# ```

## Base Graphics (Default Settings) ##
# Score and Loading Plots with base graphics
# ```{r cars_regressions, dev='svg', warning=FALSE, fig.height=12}

PCAcolors <- c("#66c2a5","#fc8d62","#8da0cb")[as.integer(vint)]


PCAscores <- vino_PCA$x
PCAloadings <- vino_PCA$rotation


par(mfrow=c(2,1))
plot(PCAscores[,1:2],  # x and y data
     pch=21,           # point shape
     col=PCAcolors,    # point border color
     bg=PCAcolors,     # point color
     cex=1.5,          # point size
     main="Scores"     # title of plot
)
legend("topright",                                # position of legend
       legend=levels(vint),                       # legend display
       pch=21,                                    # point shape
       pt.bg=c("#66c2a5","#fc8d62","#8da0cb"),    # point colors
       pt.cex=1.5,                                # point size
       col = c("#66c2a5","#fc8d62","#8da0cb")    # point border color
)
plot(PCAloadings[,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1,               # point size
     main="Loadings"      # title of plot
)
text(PCAloadings[,1:2],             # sets position of labels
     labels=rownames(PCAloadings)   # print labels
)


# ```

# We see that the text is on top of the points, which makes the labels a little unreadable.  We can fix this using the maptools package.  Also, we can add 95% confidence ellipses to the groups in the score plot and clean up the plots using a custom function.

## ellipseplot() function ##

# ```{r cars_ellipse}

## Cut and paste into console

## Customize yaxis range to makes sure axis ticks cover data
## Axes ticks do not always cover data range in R plots - reviewer did not like!
plotat <- function(RANGE) {
	if(length(RANGE) != 2) stop("RANGE argument must have a length of 2")
	if(RANGE[1] > RANGE[2]) stop("First element in RANGE must be smaller then second element")
	prettyres <- pretty(sprintf("%.2f",RANGE[1]):sprintf("%.2f",RANGE[2]), 7)
	while((min(prettyres) < RANGE[1]) == FALSE) {
		prdiff <- prettyres[2] - prettyres[1]
		prettyres[length(prettyres) + 1] <- prettyres[1] - prdiff
		prettyres <- sort(prettyres)
	} 
	while((max(prettyres) > RANGE[2]) == FALSE) {
		prdiff <- prettyres[2] - prettyres[1]
		prettyres[length(prettyres) + 1] <- prettyres[length(prettyres)] + prdiff
		prettyres <- sort(prettyres)	
	}	
	plotticks <- as.numeric(sprintf("%.2f",prettyres))
	plotticks
}

## ellipseplot function
ellipseplot <- function(x, y, factr, 
						elev=0.95, # Ellipse probability level
						legpos=c("topright","topleft","bottomleft","bottomleft"), # Legend position
					    pcol=NULL, # manual addition of colors, must meet length of factors
						cexsize=1, # point size
						ppch=21, # Point type, must meet length of factors
						legcexsize=2, # legend font size
						legptsize=2, # legend point size
						pbgcol=TRUE,
						axissize=1, 
						linewidth=1, 
						font=1) {
	require(plyr)
	require(car)
	## Set factor levels
	if(is.factor(factr)) {
		f <- factr
	} else {
		f <- factor(factr, levels=unique(as.character(factr)))
	}
	intfactr <- as.integer(f) # Set integer vector that matches factor levels
	# Checking to make sure length of ppch equals number of factor levels
	if((length(ppch) > 1 & length(unique(intfactr)) != length(ppch))) stop("Can only increase point shape if equal to factor levels")
	
	## Get data for ellipses
	edf <- data.frame(LV1 = x, LV2=y, factr = f) # create data frame with data and factor
	ellipses <- dlply(edf, .(factr), function(x) {
		LV1 <- x[,1]
		LV2 <- x[,2]
		dataEllipse(LV1, LV2, levels=elev, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
	})
	## Get range of x and y data
	xrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,1])), min(x), max(x))))
	yrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,2])), min(y), max(y))))
	

	## Set colors for plots
	if(is.null(pcol) != TRUE) { # If colors are supplied by user
		ptcol <- pcol
		pgcol <- paste(pcol, "7e", sep="") # adds opaqueness
	} else { # Default
		pgcol <- c("#e41a1c7e","#377eb87e","#4daf4a7e","#984ea37e","#807f7d7e") # Defaults at 5 colors
		ptcol <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#807f7d") # For opaqueness
	}
	# Plotting graphic
	plot(x,y, type="n", xlab="", ylab="", main="", xlim=range(xrange), ylim=range(yrange), axes=FALSE)
	axis(1, at=xrange, labels=xrange, cex.axis=axissize,lwd=linewidth, font=font)
	axis(2, las=2, cex.axis=axissize,lwd=linewidth, font=font)
	box(lwd=linewidth, font=font)
	abline(h=0, v=0, col="gray", lty=2) # Adds lines at 0
	legpch <- c() # vector to collect legend pch data
	legcol <- c() # vector to collect legend col data
	## Not sure why I split this up, might have been an artifact of an older version.
	## Adds points, ellipse, and determines color specifications for legend 
	if(pbgcol==TRUE)  {
		for(i in 1:length(unique(intfactr))){
			points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col=ptcol[i], bg=ptcol[i],cex=cexsize)
			polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
			legpch[i] <- ppch[i]
			legcol[i] <- ptcol[i]
		}
	} else {
		for(i in 1:length(unique(intfactr))){
			points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col="black", bg=ptcol[i],cex=cexsize)
			polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
			legpch[i] <- ppch[i]
			legcol[i] <- ptcol[i]		
		}
	}
	## Legend
	legend(x=legpos, legend=levels(f), pch=legpch, 
		pt.bg=legcol, col=legcol, bty="n", border=FALSE, pt.cex=legptsize, cex=legcexsize)
}	

## Axis legends for PCA output using prcomp() function
PCAvarAxis <- function(PCA, decimal=1) {
	pcavar <- round((PCA$sdev^2)/sum((PCA$sdev^2)),3)*100   #Calculate % variance explained
	PC1var <- paste("Principal Component 1 (", pcavar[1], "%)", sep="")
	PC2var <- paste("Principal Component 2 (", pcavar[2], "%)", sep="")
	PC3var <- paste("Principal Component 3 (", pcavar[3], "%)", sep="")
	PC4var <- paste("Principal Component 4 (", pcavar[4], "%)", sep="")	
	PC5var <- paste("Principal Component 5 (", pcavar[5], "%)", sep="")		
	return(list(PC1=PC1var, PC2=PC2var, PC3=PC3var, PC4=PC4var, PC5=PC5var))
}	
# ```

## Base Graphics (Cleaner Plots) ##
# Score plot using custom ellipseplot() function
# Loadings plot using base defaults, but cleaned up...
# ```{r PCA_clean, dev='svg', warning=FALSE, fig.height=12}

## Capture % variance explained from PCA
explainPCAvar <- PCAvarAxis(vino_PCA)

par(mfrow=c(2,1))
ellipseplot(PCAscores[,1],                          # data for x-axis
            PCAscores[,2],                          # data for y-axis
            vint,                                   # factor with classes
            pcol=c("#66c2a5","#fc8d62","#8da0cb"),  # colors for plotting (must match # of factors)
            pbgcol=FALSE,                           # point borders black?
            cexsize=1.5,                            # size of points 
            ppch=c(21:23),                          # shape of points (must match # of factors)
            legpos="bottomright",                   # position of legend           
            legcexsize=1.5,                         # legend text size
            legptsize=1.5,                          # legend point size 
            axissize=1.5,                           # Set axis text size
            linewidth=1.5                           # Set axis line size
)                         
title(xlab=explainPCAvar[["PC1"]],    # % variance explained on PC1
      ylab=explainPCAvar[["PC2"]],    # % variance explained on PC2 
      main="Scores",                  # Title
      cex.lab=1.5,                    # size of label text
      cex.main=1.5                    # size of title text
)

plot(PCAloadings[,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1.5,             # point size
    # type="n",           # does not plot points
     axes=FALSE,          # does not print axes
     xlab="",             # removes x label
     ylab=""              # removes y label
)
pointLabel(PCAloadings[,1:2],             # set position of labels
           labels=rownames(PCAloadings),  # print labels
           cex=1.5                          # set size of label
) # pointLabel will try to position the text around the points
axis(1,                 # display x-axis
     cex.axis=1.5,      # set size of text
     lwd=1.5            # set size of axis line
)
axis(2,                 # display y-axis
     las=2,             # argument sets direction of text, 2 is perpendicular
     cex.axis=1.5,      # set size of text
     lwd=1.5            # set size of axis line
)
box(lwd=1.5             # line width of box surrounding plot
)
title(xlab=explainPCAvar[["PC1"]],    # % variance explained on PC1
      ylab=explainPCAvar[["PC2"]],    # % variance explained on PC2 
      main="Loading",                 # Title
      cex.lab=1.5,                    # size of label text
      cex.main=1.5                    # size of title text
)
# ```

# Author: [Brian Piccolo](https://scholar.google.com/citations?user=XnewiIUAAAAJ&hl=en&oi=ao)

# Highly recommend [Chemometrics with R Multivariate Data Analysis in the Natural Sciences and Life Sciences]() by Ron Wehrens for those interested in PCA and Multivariate Analysis. 

# Created with [knitr & RMarkdown](https://www.rstudio.com/products/rpackages/) and [knitrBootstrap](https://github.com/jimhester/knitrBootstrap).