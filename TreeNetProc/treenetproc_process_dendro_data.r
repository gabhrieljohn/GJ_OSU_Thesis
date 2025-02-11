# install if necessary
#packages <- (c("devtools","zoo","chron","dplyr","viridis", "RCurl", "DT"))
#install.packages(setdiff(packages, rownames(installed.packages())))
#devtools::install_github("treenet/treenetproc")

library(treenetproc)
library(zoo)
library(chron)
library(viridis)
library(dplyr)

# helper functions
left <-  function(string, char){substr(string, 1,char)}
right <-  function (string, char){substr(string,nchar(string)-(char-1),nchar(string))}

setwd("~/Documents/Ameriflux/Me6/dendrometers")

# readin dendrometer data (ts,)
ID = '137' #dendrometer ID number
all_data<-read.table('./data/data_92223137_2022_08_29_0.csv', header=F, sep=";", dec=',')
dendro_data = all_data[,c(2,7)]
colnames(dendro_data) <- c('datetime','value')

#fix timestamp
ts = as.POSIXct(dendro_data['datetime'][,],format="%Y.%m.%d %H:%M",tz="GMT")
dendro_data_L0<-cbind(dendro_data["value"]+0,ts,series=ID)

# level 1 processing (time align)
#?treenetproc::proc_L1
dendro_data_L1 <- proc_L1(data_L0 = dendro_data_L0,
                          reso = 60 ,
                          #input = "wide",
                          date_format ="%Y-%m-%d %H:%M:%S",
                          tz = "GMT")
head(dendro_data_L1)

# level 1 processing of temperature data
temp_data_L0 <- cbind(all_data[,c(4,6)],ts,series=ID)
colnames(temp_data_L0) <-c('value','empty','ts','series')

# align data
temp_data_L1 <- proc_L1(data_L0 = temp_data_L0,
                        reso = 60,
                        date_format ="%Y-%m-%d %H:%M:%S",
                        tz = "GMT")
head(temp_data_L1)


# process and plot
par(mfrow=c(1,1))
par(mar = c(5, 5, 5, 5))

# detect errors
dendro_data_L2 <- proc_dendro_L2(dendro_L1 = dendro_data_L1,
                                 temp_L1 = temp_data_L1,
                                 tol_jump = 15,
                                 plot = TRUE,
                                 tz="GMT")
# check the data
head(dendro_data_L2)


# plot minimum daily twd against day of year
par(mfrow=c(1,1))
par(mar = c(5, 5, 5, 5))

plot(1, 1,
     ylim=c(-1,max(dendro_data_L2$twd,na.rm=T)),
     xlim=c(120,240),
     ylab=expression("twd ("*mu*"m)"),
     xlab="Day of year",
     col="white")

col_sel<-c("cyan","darkorange","purple")

for(y in c(1:length(unique(left(dendro_data_L2$ts,4))))){
  # selected year
  sel<-dendro_data_L2[which(left(dendro_data_L2$ts,4)==unique(left(dendro_data_L2$ts,4))[y]),]
  # calc twd
  twd<-suppressWarnings(aggregate(sel$twd,list(as.Date(sel$ts)),min,na.rm=T))
  twd$doy<-as.numeric(strftime(as.Date(twd$Group.1), format = "%j"))
  
  # clean
  twd[which(twd$x=="Inf"),"x"]<-NA
  
  lines(twd$doy,twd$x,col=col_sel[y],lwd=1.5)
  twd[which(is.na(twd$x)==T),"x"]
  polygon(c(c(0,twd$doy),c(rev(twd$doy),0)),
          c(c(0,twd$x),rep(0,nrow(twd)+1)),
          col=rgb(0,0,0,0.1),
          border=rgb(0,0,0,0))
}


# plot growth against day of year
par(mfrow=c(1,1))
par(mar = c(5, 5, 5, 5))

plot(1, 1,
     ylim=c(-1,max(dendro_data_L2$gro_yr,na.rm=T)),
     xlim=c(120,240),
     ylab=expression("growth  ("*mu*"m)"),
     xlab="Day of year",
     col="white")

col_sel<-c("darkorange","purple")

for(y in c(1:length(unique(left(dendro_data_L2$ts,4))))){
  # selected year
  sel<-dendro_data_L2[which(left(dendro_data_L2$ts,4)==unique(left(dendro_data_L2$ts,4))[y]),]
  # calc twd
  gro<-suppressWarnings(aggregate(sel$gro_yr,list(as.Date(sel$ts)),max,na.rm=T))
  gro$doy<-as.numeric(strftime(as.Date(gro$Group.1), format = "%j"))
  
  # clean
  gro[which(gro$x=="Inf"),"x"]<-NA
  
  lines(gro$doy,gro$x,col=col_sel[y],lwd=1.5)
  gro[which(is.na(gro$x)==T),"x"]

}


