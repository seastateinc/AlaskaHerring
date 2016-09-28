# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# R4HAM.R 
# Author: Steven Martell
# Sep 27, 2016
# Objective:
# 	This is a simple example of how I structure my R-scripts so they are
#   easy to maintain and portable.
#
#
#
# ---------------------------------------------------------------------------- #

#
# LIBRARIES
#
library(ggplot2)
library(dplyr)
library(tidyr)


#
#	OTHER R-SRCIPTS
#
source(file.path("globals.R"))
source(file.path("plotIndex.R"))
source(file.path("ham.R"))


# ---------------------------------------------------------------------------- #
# READ MODEL OUTPUTS AND STORE IN A LIST
# ---------------------------------------------------------------------------- #
# MODELDIRS will recursively look in all folders specified in the file path.
"~/Library/Mobile Documents/com~apple~CloudDocs/SSI/ADFG_HerringModel/R"
.MODELDIRS <- list.dirs(file.path("../models_2015"))

.getModelOutput <- function(fp) {
	#filter on the correlation file to ensure model converged
	if(file.exists(paste0(fp,"/ham.cor"))) {
		output <- read.admb(paste0(fp,"/ham"))
		output$Model <- fp
		return (output)
	} else {
		return (NULL)
	}
}

M <- lapply(.MODELDIRS,.getModelOutput)
M <- M[!sapply(M, is.null)]  # Remove NULL models from list (unconverged).


# ---------------------------------------------------------------------------- #
# PLOT DATA INPUTS
# ---------------------------------------------------------------------------- #

plot.all.data <- function(D) {
	# print(names(D))
	p.data <- NULL
	p.data$catch   <- plot.catch(D, nm = "data_ct_raw", y="Catch (tons)")
	p.data$sp_waa  <- plot.waa(D, nm = "data_sp_waa")
	p.data$cm_waa  <- plot.waa(D, nm = "data_cm_waa")
	p.data$sp_comp <- plot.comp(D, nm = "data_sp_comp")
	p.data$cm_comp <- plot.comp(D, nm = "data_cm_comp")

	cat("Finished plot.all.data\n")
	return( p.data )
}

# P.data <- plot.all.data(M[[1]])

data.Plots <- lapply(M,plot.all.data)


# ---------------------------------------------------------------------------- #
# PLOT MLE MODEL OUTPUT
# ---------------------------------------------------------------------------- #

plot.all.output <- function(D) {
	p.output <- NULL
	p.output$ssb <- plot.ssb(D)
	p.output$egg <- plot.eggIndex(D,fit=TRUE)  
	p.output$resd_cct <- plot.resd(D,nm="resd_catch")
	p.output$resd_egg <- plot.resd(D,nm="resd_egg_dep")
	p.output$resd_rec <- plot.resd(D,nm="resd_rec")
	p.output$resd_csp <- plot.resd(D,nm="resd_sp_comp")
	p.output$resd_ccm <- plot.resd(D,nm="resd_cm_comp")

	cat("Finished plot.all.output\n")
	return(p.output)
}

output.Plots <- lapply(M,plot.all.output)