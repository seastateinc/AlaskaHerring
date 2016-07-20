# ---------------------------------------------------------------------------- #
# Herring Age Model, or ham 
# June 29, 2016
# Author: Steven Martell
# email: steve@seastateinc.com
# phone: 808 280-2548
# ---------------------------------------------------------------------------- #
library(ggplot2)
library(dplyr)
library(tidyr)

# Read in the data from the model report, par, and cor files.
source("./globals.R")
D <- read.admb("../src/ham")

# ---------------------------------------------------------------------------- #
# DATA SECTION
# ---------------------------------------------------------------------------- #

plot.catch <- function(D=D, nm = "data_ct_raw",...) {
	df <- as.data.frame(D[[nm]])
	colnames(df) <- c("Year","Catch","log.se")
	df <- df %>% 
				mutate(std = 1.96*sqrt(log(log.se+1))) %>%
				mutate(lower=Catch-std*Catch,upper=Catch+std*Catch)

	ggplot(df,aes(Year,Catch)) +
		geom_pointrange(aes(ymin = lower, ymax = upper),size=0.5,fatten=2) + 
		labs(...)
}


plot.waa <- function(D=D, nm = "data_sp_waa",...) {
	df <- as.data.frame(D[[nm]])
	colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
	gdf <- 	gather(df,"Age","Weight.at.age",-Year) %>% 
					transform(Age=as.integer(Age)) %>% 
					mutate(Cohort=as.factor(Year-Age))

	ggplot(gdf,aes(Year,Weight.at.age,color=Cohort)) + 
		geom_line(alpha=0.90) + 
		geom_point(alpha=0.5,aes(fill=Cohort),show.legend=FALSE,size=0.5) +
		labs(x="Year",y="Weight-at-age (grams)",color="Cohort") +
		guides(col = guide_legend(ncol = 9)) +
		theme(legend.position="bottom")
}

plot.comp <- function(D=D, nm = "data_cm_comp",...) {
	df <- as.data.frame(D[[nm]])
	df[df==-9] <- NA
	colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
	gdf <- 	gather(df,"Age","Proportion",-Year) %>%
					transform(Age=as.integer(Age)) %>%
					mutate(Cohort=as.factor(Year-Age))

	ggplot(gdf,aes(Year,Age,color=Cohort)) + 
		geom_point(alpha=0.5,aes(size=abs(Proportion)),show.legend=FALSE) +
		scale_size_area(max_size=8) + 
		labs(x="Year",...)
}


plot.ssb <- function(D=D){
	qplot(D$year,D$ssb/1000,geom="line") + ylim(c(0,NA)) +
	labs(x="Year",y="Female Spawning Stock Biomass (1000 mt)")
}

plot.datafit <- function(D=D, sfx="egg_dep") {
	data <- paste0("data_",sfx)
	data <- as.data.frame(D[[data]])
	colnames(data) <- c("year","index","log.se")

	pred <- paste0("pred_",sfx)
	pred <- as.data.frame(cbind(D$year,D[[pred]]))
	colnames(pred) <- c("year","pred")

	df   <- right_join(data,pred) %>%
					mutate(std = 1.96*sqrt(log(log.se+1))) %>%
					mutate(lower=index-std*index,upper=index+std*index) %>%
					tbl_df()

	# print(head(df))

	ggplot(df,aes(year,index)) + 
	geom_pointrange(aes(ymin = lower, ymax = upper),size=0.5,fatten=2) +
	geom_line(aes(year,pred),alpha=0.8) + 
	labs(x="Year",y="Egg Deposition (trillions)")

}

plot.resd <- function(D=D, nm = "resd_cm_comp", ...) {

	# Dealing with composition data.
	if( grepl("comp",nm) ){
		df <- as.data.frame(cbind(D[['mod_syr']]:D[['mod_nyr']],D[[nm]]))
		# df[df==-9] <- NA
		colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
		gdf <- 	gather(df,"Age","Residual",-Year) %>%
						transform(Age=as.integer(Age)) %>%
						mutate(Cohort=as.factor(Year-Age))

		ggplot(gdf,aes(Year,Age,color=factor(sign(Residual)))) + 
			geom_point(alpha=0.5,aes(size=abs(Residual)),show.legend=TRUE) +
			scale_size_area(max_size=8) + 
			labs(x="Year",color="Sign",size="Residual",...)

	} else if( grepl("egg",nm) ) {
		df <- as.data.frame(cbind(D[['year']],D[[nm]]))
		colnames(df) <- c("Year","Residual")

		ggplot(df,aes(Year,Residual)) + geom_point() +
		geom_segment(aes(x = Year, xend = Year, y = 0, 
		                 yend = Residual),data=df,size=0.2)

	} else if( grepl("rec",nm) ) {
		df <- as.data.frame(cbind(D[['rec_years']],D[[nm]]))
		colnames(df) <- c("Year","Residual")

		ggplot(df,aes(Year,Residual)) + geom_point() +
		geom_segment(aes(x = Year, xend = Year, y = 0, 
		                 yend = Residual),data=df,size=0.2)
	} else if( grepl("catch",nm) ) {
		df <- as.data.frame(cbind(D[['year']],D[[nm]]))
		colnames(df) <- c("Year","Residual")

		ggplot(df,aes(Year,Residual)) + geom_point() +
		geom_segment(aes(x = Year, xend = Year, y = 0, 
		                 yend = Residual),data=df,size=0.2)
	}
}

# ---------------------------------------------------------------------------- #
# PLOTS FOR DATA SECTION
# ---------------------------------------------------------------------------- #
# prefix with d for data
# d1 » catch time series
# d2 » spawn weight-at-age
# d3 » commercial weight-at-age
# d4 » commercial age-proportions
# d5 » spawn sample age-proportions
d1 <- plot.catch(D, nm = "data_ct_raw",y="Catch (tons)")
d2 <- plot.waa(D,nm = "data_sp_waa")
d3 <- plot.waa(D,nm = "data_cm_waa")
d4 <- plot.comp(D,nm = "data_cm_comp")
d5 <- plot.comp(D,nm = "data_sp_comp")




