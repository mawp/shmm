# Spatial hidden Markov model (SHMM)
#    Copyright (C) 2015-2016  Martin Waever Pedersen, mawp@dtu.dk or wpsgodd@gmail.com
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' @name get.dat.sun
#' @title Extract sun info from psat xlsx-file.
#' @details Extract sun info from psat xlsx-file.
#' @param fn Path to input xlsx-file.
#' @param filter_thres Threshold to use in the filtering process. Smaller values give more restrictive filtering. 0.2 is a good starting point.
#' @return Dataframe (dat_s) with sun info from xlsx-file. Time of sunrise and -set along with depth at time of sunrise/set. Columns sr_rm and ss_rm contain result from filter (0=keep, 1=remove).
get.dat.sun <- function(fn, filter_thres=0.2, showPlot=TRUE){
	dat_s <- openxlsx::read.xlsx(fn, sheet="Sunrise and Sunset Times", startRow=2, cols=1:5)
	colnames(dat_s) <- c('date','sr','depth_sr', 'ss','depth_ss')
	dat_s$date <- as.POSIXct((dat_s$date*(24*60*60)), tz="GMT", origin="1899-12-30")
	dat_s$depth_sr <- as.numeric(dat_s$depth_sr)
	dat_s$depth_ss <- as.numeric(dat_s$depth_ss)
	dat_s$sr_rm <- 0
	dat_s$ss_rm <- 0

	filter_ss <- filter.sun(dat_s, 'ss',filter_thres)
	filter_sr <- filter.sun(dat_s, 'sr',filter_thres)
	
	dat_s$sr_rm[filter_sr[['rems']]] <- 1
	dat_s$ss_rm[filter_ss[['rems']]] <- 1
	
	
	dat_s$day_length[dat_s$sr_rm == 0 & dat_s$ss_rm == 0] <- dat_s$ss[dat_s$sr_rm == 0 & dat_s$ss_rm == 0] - dat_s$sr[dat_s$sr_rm == 0 & dat_s$ss_rm == 0]

	if(showPlot){
		plot(dat_s$date, dat_s$sr, ylim=c(0,1))
		points(dat_s$date, dat_s$ss)
		points(dat_s$date[dat_s$sr_rm == 1], dat_s$sr[dat_s$sr_rm == 1], col="red")
		points(dat_s$date[dat_s$ss_rm == 1], dat_s$ss[dat_s$ss_rm == 1], col="red")
		points(dat_s$date, dat_s$day_length, col="green")
	}
	
	n_rem <- sum(filter_ss[['n']],filter_sr[['n']])
	print(paste0(n_rem, " observations removed by filter"))
	
	return(dat_s)

}

#' @name filter.sun
#' @title Filter sun-info - support function for use in get.dat.sun().
#' @details Filter sun-info - support function for use in get.dat.sun().
#' @param dat_s Dataframe use internal in get.dat.sun().
#' @param direction Indicates whether sunrise (sr) or sunset (ss) should be fitlered. 
#' @param filter_thres Threshold to use in the filtering process. Smaller values give more restrictive filtering. 0.2 is a good starting point.
#' @return List of observations to be excluded incl. NAs (rems) and number of items filtered out not including NA's. 
filter.sun <- function(dat_s, direction, filter_thres){
	if(!direction %in% c('ss','sr')) {print("Direction must be either 'ss' or 'sr'"); return();}
	rems <- which(is.na(dat_s[direction]))
	
	Emax <- 1
	i <- 0
	wgt <- exp(dat_s[paste0("depth_", direction)])[,1]
	wgt[is.na(wgt)] <- exp(-5)
	while(abs(Emax) > filter_thres){
		if(length(rems) == 0){
			m <- mgcv::gam(sr ~ s(as.numeric(date)), data=dat_s, na.action=na.omit, weights=wgt, family="scat")
		} else {
			# wgt[rems,1] <- NA
			suppressWarnings(m <- mgcv::gam(sr ~ s(as.numeric(date)), data=dat_s[-rems,], na.action=na.omit, weights=wgt[-rems], family="scat"))
		}
		E <- resid(m, type="response")
		p <- predict(m)
		Emax_idx <- which.max(abs(E))
		Emax <- E[Emax_idx]
		Emax_idx_org <- seq(1:nrow(dat_s))[-rems][Emax_idx]
		rems <- c(rems, Emax_idx_org)
		i <- i+1
	}
	return(list(rems=rems, n=i))
}




#' @name get.dat.dt
#' @title Extract depth and temperature info from PSAT xlsx-file.
#' @details Extract depth and temperature info from PSAT xlsx-file.
#' @param fn Path to input xlsx-file.
#' @return Dataframe (dat_dt) with depth and temperature info from PSAT xlsx-file.
get.dat.dt <- function(fn){
	#get temperature data from psat xlsx-file
	dat_t <- openxlsx::read.xlsx(fn, sheet="Temp Data", startRow=2, cols=1:3)
	colnames(dat_t) <- c('datetime','temp_raw','temp')
	dat_t$datetime <- as.POSIXct((dat_t$datetime*(24*60*60)), tz="GMT", origin="1899-12-30")
	
	#get depth data from psat xlsx-file
	dat_d <- openxlsx::read.xlsx(fn, sheet="Press Data", startRow=2, cols=1:4)
	colnames(dat_d) <- c('datetime','press','gain','depth')
	dat_d$datetime <- as.POSIXct((dat_d$datetime*(24*60*60)), tz="GMT", origin="1899-12-30")
	
	#combine temperature and depth
	dat_dt <- merge(dat_t[c('datetime','temp')], dat_d[c('datetime','depth')], by='datetime', all=TRUE)

	return(dat_dt)
}

#' @name get.dat.sst
#' @title Estimate SST based on data from get.dat.dt().
#' @details Estimate SST based on data from get.dat.dt().
#' @param dat_dt Object from get.dat.dt().
#' @return Dataframe containing SST dat based on depth and temperature info from PSAT xlsx-file.
get.dat.sst <- function(dat_dt){
	dat_dt$day <- as.character(strptime(dat_dt$datetime, format="%Y-%m-%d", tz="GMT"))
	days <- unique(dat_dt$day)
	dat_sst <- plyr::aaply(days, 1, .fun=function(k) {
		dat_k <- dat_dt[which(dat_dt$day == k),]
		if(sum(!is.na(dat_k$depth[which(!is.na(dat_k$temp))])) > 0) {
			min_depth <- max(dat_k$depth[which(!is.na(dat_k$temp))], na.rm=TRUE)
			min_depth_idx <- which(dat_k$depth == min_depth)
			avg <- mean(dat_k$temp[min_depth_idx], na.rm=TRUE)
			med <- median(dat_k$temp[min_depth_idx], na.rm=TRUE)
			sd <- sd(dat_k$temp[min_depth_idx], na.rm=TRUE)
			avg_upper10 <- mean(dat_k$temp[dat_k$depth > -10], na.rm=TRUE)
			n <- sum(!is.na(dat_k$temp[min_depth_idx]))
			out_k <- c(min_depth, avg, med, sd, n, avg_upper10)
		} else {
			out_k <- c(NA,NA,NA, NA,NA,NA)
		}
		return(out_k)
	})
	
	dat_sst <- cbind(days,as.data.frame(dat_sst))
	colnames(dat_sst) <- c('day', 'min_depth', 'avg', 'med', 'sd', 'n', 'avg_upper10')
	dat_sst$day <- as.POSIXct(dat_sst$day, tz="UTC")
	return(dat_sst)
}
