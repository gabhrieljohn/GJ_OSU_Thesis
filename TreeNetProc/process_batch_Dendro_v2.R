#create script for running a batch of dendrometers to run in treenetproc

# install if necessary
packages <- (c("devtools","zoo","chron","dplyr","viridis", "RCurl", "DT","lubridate"))
install.packages(setdiff(packages, rownames(installed.packages())))
devtools::install_github("treenet/treenetproc")


library(treenetproc)
library(zoo)
library(chron)
library(viridis)
library(dplyr)
library(ggplot2)
library(lubridate)
library(generics)

getwd()

# helper functions
left <-  function(string, char){substr(string, 1,char)}
right <-  function (string, char){substr(string,nchar(string)-(char-1),nchar(string))}

#d = '940'
#date <- '_2023_07_27_0'
#file = ".csv"
#path = 'Raw/Soda_dendrometry_8_2023/data_92223'

#fun_insert <- function(x,pos,insert) {
# gsub(paste0("^(.{",pos,"})(.*)$"),
#     paste0("\\1",insert,"\\2"),
#    x)
#}
#

#path = fun_insert(x= path, pos = 38, insert = d)
#new_path = fun_insert(x=path, pos = 41, insert = date)
#full_path = fun_insert(x=new_path,pos= 54, insert = file)

#full_path
#all_data_940 <-read.table(full_path, header=F, sep=";", dec=',')

#all_data_940<-read.table('Raw/Soda_dendrometry_8_2023/data_92223947_2023_07_27_0.csv', header=F, sep=";", dec=',')

#function for reading in and cleaning data
make_dendro_l0 <- function(d,date,path) {
  
  #read in data
  
  d = d
  date <- paste0("_",date,"_0",sep="")
  file = ".csv"
  path = paste0(path,'data_92223',sep = "")
  
  fun_insert <- function(x,pos,insert) {
    gsub(paste0("^(.{",pos,"})(.*)$"),
         paste0("\\1",insert,"\\2"),
         x)
  }
  
  new_path = fun_insert(x= path, pos = 38, insert = d)
  new_path2 = fun_insert(x=new_path, pos = 41, insert = date)
  full_path = fun_insert(x=new_path2,pos= 54, insert = file)
  
  
  all_data_d <-read.table(full_path, header=F, sep=";", dec=',')
  str(all_data_d)
  dendro_data_d = all_data_d[,c(2,7)]
  colnames(dendro_data_d) <- c('datetime','value')
  
  #fix timestamp
  ts = as.POSIXct(dendro_data_d['datetime'][,],format="%Y.%m.%d %H:%M",tz="GMT")
  #create dendro l0 and save
  dendro_L0 <-cbind(dendro_data_d["value"]+0,ts,series=d)
  output <- paste0("Cleaning/Soda_8_2023/l0/l0_",d,".csv", sep="")
  write.csv(dendro_L0, file = output)
  
  #create temp l0
  temp_data_L0 <- cbind(all_data_d[,c(4,6)],ts,series=d)
  colnames(temp_data_L0) <-c('value','empty','ts','series')
  str(temp_data_L0)
  temp_data_L0$value <- as.numeric(temp_data_L0$value)
  
  output2 <- paste0("Cleaning/Soda_8_2023/l0/temp_l0_",d,".csv", sep="")
  write.csv(temp_data_L0, file = output2)
  
  
}
make_dendro_l0("940","2023_07_27","Raw/Soda_dendrometry_8_2023/")

#function for time alignment
make_dendro_l1 <- function(d,input_dir,output_dir){
  #read in data
  
  d = d
  file = ".csv"
  path = paste0(path,"l0",sep = "")
  path2 = paste0(path,"temp_l0",sep = "")
  
  
  fun_insert <- function(x,pos,insert) {
    gsub(paste0("^(.{",pos,"})(.*)$"),
         paste0("\\1",insert,"\\2"),
         x)
  }
  
  new_path <- fun_insert(x=path, pos=27, insert = d)
  input <- fun_insert(x=new_path, pos=30, insert = file)
  
  new_path2 <- fun_insert(x=path2, pos=32, insert = d)
  input2 <- fun_insert(x=new_path2, pos=35, insert=file)
  
  l0 <- read.csv(file = input)
  l0 <- l0[,c(2,3,4)]
  l0$value <- as.numeric(l0$value)
  l0$ts = ymd_hms(l0$ts)
  
  
  temp_l0 <- read.csv(file = input2)
  temp_l0 <- temp_l0[,c(2,3,4,5)]
  temp_l0$value <- as.numeric(temp_l0$value)
  temp_l0$ts = ymd_hms(temp_l0$ts)
  
  
  dendro_L1 <- proc_L1(data_L0 = l0,
                       reso = 60 ,
                       #input = "wide",
                       date_format ="%Y-%m-%d %H:%M:%S",
                       tz = "GMT")
  head(dendro_L1)
  dendro_L1 <- dendro_L1 %>% group_by(ts) %>% filter(row_number() == 1)
  output <- paste0(output_dir,"l0_",d,".csv", sep="")
  write.csv(dendro_L1, file = output) 
  
  
  temp_L1 <- proc_L1(data_L0 = temp_l0,
                     reso = 60,
  )
  head(temp_L1)
  #remove duplicates
  temp_L1 <- temp_L1 %>% group_by(ts) %>% filter(row_number() == 1)
  
  output2 <- paste0("temp_l1_",d,".csv", sep="")
  write.csv(temp_L1, file = output2) 
  
  
}

#function for outlier correction

make_dendro_l2 <- function(d){
  
  #define path
  d = d
  file = ".csv"
  path = 'Cleaning/Soda_8_2023/l1/l1_'
  path2 = 'Cleaning/Soda_8_2023/l1/temp_l1_'
  
  fun_insert <- function(x,pos,insert) {
    gsub(paste0("^(.{",pos,"})(.*)$"),
         paste0("\\1",insert,"\\2"),
         x)
  }
  
  new_path <- fun_insert(x=path, pos=27, insert = d)
  input <- fun_insert(x=new_path, pos=30, insert = file)
  
  new_path2 <- fun_insert(x=path2, pos=32, insert = d)
  input2 <- fun_insert(x=new_path2, pos=35, insert=file)
  
  #read in data
  l1 <- read.csv(file = input)
  l1 <- l1[,c(2,3,4,5)]
  l1$value <- as.numeric(l1$value)
  l1$ts = ymd_hms(l1$ts)
  
  
  temp_l1 <- read.csv(file = input2)
  temp_l1 <- temp_l1[,c(2,3,4,5,6)]
  temp_l1$value <- as.numeric(temp_l1$value)
  temp_l1$ts = ymd_hms(temp_l1$ts)
  
  dendro_l2 <- proc_dendro_L2(dendro_L1 = l1,
                              temp_L1 = temp_l1,
                              tol_out = 10,
                              tol_jump = 50,
                              plot = FALSE,
                              frost_thr = 1,
                              plot_export = FALSE,
                              interpol = NULL,
                              frag_len = NULL)
  
  output <- paste0("Cleaning/Soda_8_2023/l2/l2_",d,".csv", sep="")
  write.csv(dendro_l2, file = output)
  
}  


#these are examples of using the functions to process data

make_dendro_l0("969","_2023_07_27_0")
make_dendro_l1("969")
make_dendro_l2("969")

make_dendro_l0("968","_2023_07_27_0")
make_dendro_l1("968")
make_dendro_l2("968")   

make_dendro_l0("967","_2023_07_27_0")
make_dendro_l1("967")
make_dendro_l2("967")

make_dendro_l0("966","_2023_07_27_0")
make_dendro_l1("966")
make_dendro_l2("966")   

make_dendro_l0("954","_2023_07_27_0")
make_dendro_l1("954")
make_dendro_l2("954")

make_dendro_l0("953","_2023_07_27_0")
make_dendro_l1("953")
make_dendro_l2("953")

make_dendro_l0("952","_2023_07_27_0")
make_dendro_l1("952")
make_dendro_l2("952")

make_dendro_l0("951","_2023_07_27_0")
make_dendro_l1("951")
make_dendro_l2("951")

make_dendro_l0("950","_2023_07_27_0")
make_dendro_l1("950")
make_dendro_l2("950")

make_dendro_l0("948","_2023_07_27_0")
make_dendro_l1("948")
make_dendro_l2("948")

make_dendro_l0("947","_2023_07_27_0")
make_dendro_l1("947")
make_dendro_l2("947")

make_dendro_l0("946","_2023_07_27_0")
make_dendro_l1("946")
make_dendro_l2("946")

make_dendro_l0("940","_2023_07_27_0")
make_dendro_l1("940")
make_dendro_l2("940")

make_dendro_l0("939","_2023_07_27_0")
make_dendro_l1("939")
make_dendro_l2("939")

make_dendro_l0("938","_2023_07_27_0")
make_dendro_l1("938")
make_dendro_l2("938")

make_dendro_l0("937","_2023_07_27_0")
make_dendro_l1("937")
make_dendro_l2("937")

make_dendro_l0("936","_2023_07_27_0")
make_dendro_l1("936")
make_dendro_l2("936")

make_dendro_l0("935","_2023_07_27_0")
make_dendro_l1("935")
make_dendro_l2("935")

make_dendro_l0("933","_2023_07_27_0")
make_dendro_l1("933")
make_dendro_l2("933")

make_dendro_l0("931","_2023_07_27_0")
make_dendro_l1("931")
make_dendro_l2("931")

make_dendro_l0("906","_2023_07_27_0")
make_dendro_l1("906")
make_dendro_l2("906")
