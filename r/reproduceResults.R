reproduceResults <- function()
{
  switch(Sys.info()[['sysname']],
         Windows= {pathSymbol <- '\\'},
         Linux  = {pathSymbol <- '/'},
         Darwin = {pathSymbol <- '/'})
  require(signal)
  require(forecast)
  require(pracma)
  require(bspec)
  require(fma)
  require(expsmooth)
  require(fpp2)
  require(TSA)
  require(astsa)
  require(AER)
  require(scmamp)
  source('sazed.R')
  source('testSazed.R')
  source('callSeasonLength.R')
  source('baselines.R')
  source('fastAcf.R')
  
  
  #load CRAN dataset
  cran_data <- list()
  cran_expected <- list()
  cran_dataset_names <- list()
  cran_dataset_names[["fma"]] <- c("airpass", "beer", "bricksq", "condmilk", "elec", "fancy", "hsales",
                                "hsales2", "labour", "milk", "motion", "plastics", "qsales", "ukdeaths",
                                "usdeaths", "uselec", "writing")
  cran_dataset_names[["expsmooth"]] <- c("cangas", "enplanements", "frexport", "mcopper", "ukcars", "utility", 
                                      "vehicles", "visitors")
  cran_dataset_names[["fpp2"]] <- c("a10", "ausbeer", "auscafe", "austourists", "debitcards", "elecequip", "h02",
                                 "hyndsight", "qauselec", "qcement", "qgas", "usmelec")
  cran_dataset_names[["TSA"]] <- c("airmiles", "co2", "flow", "JJ", "oilfilters", "retail", "tempdub")
  cran_dataset_names[["astsa"]] <- c("birth", "cmort", "flu", "gas", "oil", "part", "prodn", "rec",
                                  "so2", "soi", "sunspotz", "tempr", "UnempRate")
  cran_dataset_names[["AER"]] <- c("DutchSales", "UKNonDurables")
  number_cran_dataset_names <- Reduce(sum, lapply(cran_dataset_names, length))
  dataset_libraries <- c("fma", "expsmooth", "fpp2", "TSA", "astsa", "AER")
  for (a_dataset_library in dataset_libraries) {
    library(a_dataset_library, character.only = T)
    data(list = cran_dataset_names[[a_dataset_library]])
    for (a_dataset in cran_dataset_names[[a_dataset_library]]) {
      ts_data <- get(a_dataset)
      cran_data <- append(cran_data,list(as.vector(ts_data)))
      cran_expected <- append(cran_expected, frequency(ts_data))
    }
  }
  
  #load SL dataset
  sl_data <- list()
  sl_expected <- list(88,57,3,6,5,25,214,900,10,1000,20,15,190,107,52,1000,20,35,50000,3,
                   10000,10,96,120,19,180,1332,45,19,20,140,18,30,1500,740,20,1590,10,29,5,
                   c(20,60),c(3,10),c(50,90),c(60,300),c(70,700),c(40,3498),
                   c(8,32),c(40,75),c(400,1000),c(9,20),c(5,980),c(10,64),
                   c(78,100),c(6,12,18),c(3,7),c(9,24),c(7,70,700),c(400,800),
                   c(6,20),c(500,1400,2000),
                   128,128,128,128,128,156,156,156,156,156,37,37,37,37,37,8,8,8,8,8,
                   155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,
                   4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,
                   1,1,1,1,1,1,1,1,1,1,
                   12,12,12,4,12,6,6,6,12,6,4,12,12,12,12,12,3,3,12,12,
                   12,12,12,12,12,12,12,130,12,12,12,12,12,12,130,12,12,12,12,12)
  for (i in 1:20)
  {
    sl_data <- append(sl_data,list(unlist(read.table(paste('..',pathSymbol,'test1',pathSymbol,'t',i,sep='',collapse = ' '))[2],
                                          use.names = F)))}
  for (i in 1:20)
  {
    sl_data <- append(sl_data,list(unlist(read.table(paste('..',pathSymbol,'test2',pathSymbol,'t',i,sep='',collapse = ' '))[2],
                                          use.names = F)))}
  for (i in 1:20)
  {
    sl_data <- append(sl_data,list(unlist(read.table(paste('..',pathSymbol,'test3',pathSymbol,'t',i,sep='',collapse = ' '))[2],
                                          use.names = F)))}
  for (i in 1:20)
  {
    sl_data <- append(sl_data,list(unlist(read.table(paste('..',pathSymbol,'test4',pathSymbol,'t',i,sep='',collapse = ' '))[2],
                                          use.names = F)))}
  for (i in 1:20)
  {
    sl_data <- append(sl_data,list(unlist(read.table(paste('..',pathSymbol,'test5',pathSymbol,'t',i,sep='',collapse = ' '))[2],
                                          use.names = F)))}
  for (i in 1:15)
  {
    sl_data <- append(sl_data,list(unlist(read.table(paste('..',pathSymbol,'test6',pathSymbol,'t',i,sep='',collapse = ' '))[2],
                                          use.names = F)))}
  for (i in 1:10)
  {
    sl_data <- append(sl_data,list(unlist(read.table(paste('..',pathSymbol,'test7',pathSymbol,'t',i,sep='',collapse = ' '))[2],
                                          use.names = F)))}
  for (i in 1:20)
  {
    sl_data <- append(sl_data,list(unlist(read.table(paste('..',pathSymbol,'testR',pathSymbol,'t',i,sep='',collapse = ' '))[1],
                                          use.names = F)))}
  for (i in 1:20)
  {
    sl_data <- append(sl_data,list(unlist(read.table(paste('..',pathSymbol,'testC',pathSymbol,'t',i,sep='',collapse = ' '))[1],
                                          use.names = F)))}
  
  #plot sunspots
  dev.new()
  plot(sunspot.month)
  
  #plot vehicles and correlation
  dev.new()
  par(mfrow=c(1,2))
  plot(vehicles)
  
  y <- preprocessTs(vehicles)
  y <- computeAcf(y)
  n <- length(y)
  signs <- y
  signs[which(signs < 0)] <- -1
  signs[which(signs > 0)] <- 1
  
  zero_distance_raw <- which(signs[2:n] + signs[1:(n-1)] < 2 & signs[2:n] + signs[1:(n-1)] > -2)
  interpolation <- y[zero_distance_raw] / (-1*(y[zero_distance_raw+1]-y[zero_distance_raw]))
  zero_distance_exact <- zero_distance_raw+interpolation
  z_to_z <- diff(zero_distance_exact)
  
  if (length(z_to_z) < 2)
  {
    return(1)
  }
  dens <- density(z_to_z,kernel = 'epanechnikov')
  plot(dens)
  
  #plot ensemble example
  dev.new()
  par(mfrow=c(1,3))
  y <- ts(sl_data[[157]])
  plot(y)
  s <- spec.pgram(y,detrend = F,plot = F)
  plot(s$freq,s$spec,ylim=c(0,30),type='l')
  plot(acf.fft(y),type='l')
  
  #plot package usage example
  dev.new()
  results <- c()
  for (season_length in 3:300)
  {
    results <- c(results, sazed(sin((1:1800)*2*pi/season_length)))
  }
  plot(3:300, type='l', col='black', xlab='Expected season length', ylab='Computed season length')
  lines(results, col='blue')
  
  #plot dataset histograms
  dev.new()
  par(mfrow=c(1,2))
  hist(log(unlist(lapply(cran_data,length))),breaks = 20)
  hist(log(unlist(cran_expected)))
  
  dev.new()
  par(mfrow=c(1,2))
  hist(log(unlist(lapply(sl_data,length))),breaks = 50,xlim = c(2,14))
  hist(log(unlist(sl_expected)),breaks = 50,xlim = c(0,12))
  
  
  #compute results
  findfrequency_results_cran <- c()
  seasonLength_results_cran  <- c()
  sazed_down_results_cran <- c()
  sazed_diff_results_cran <- c()
  sazed_alt_results_cran <- c()
  s_results_cran <- c()
  sa_results_cran <- c()
  ze_results_cran <- c()
  aze_results_cran <- c()
  zed_results_cran <- c()
  azed_results_cran <- c()
  
  findfrequency_results_sl <- c()
  seasonLength_results_sl  <- c()
  sazed_down_results_sl <- c()
  sazed_diff_results_sl <- c()
  sazed_alt_results_sl <- c()
  s_results_sl <- c()
  sa_results_sl <- c()
  ze_results_sl <- c()
  aze_results_sl <- c()
  zed_results_sl <- c()
  azed_results_sl <- c()
  
  for (test in cran_data)
  {
    r <- tryCatch(findfrequency(test),error = function(e){return(1)})
    findfrequency_results_cran <- c(findfrequency_results_cran,r)
    r <- tryCatch(callSeasonLength(test),error = function(e){return(1)})
    seasonLength_results_cran  <- c(seasonLength_results_cran,r)
    r <- tryCatch(sazed(test,method = 'down'),error = function(e){return(1)})
    sazed_down_results_cran <- c(sazed_down_results_cran,r)
    r <- tryCatch(sazed(test,method='diff'),error = function(e){return(1)})
    sazed_diff_results_cran <- c(sazed_diff_results_cran,r)
    r <- tryCatch(sazed(test,method = 'alt'),error = function(e){return(1)})
    sazed_alt_results_cran <- c(sazed_alt_results_cran,r)
    r <- tryCatch(S(test),error = function(e){return(1)})
    s_results_cran <- c(s_results_cran,r)
    r <- tryCatch(Sa(test),error = function(e){return(1)})
    sa_results_cran <- c(sa_results_cran,r)
    r <- tryCatch(ze(test),error = function(e){return(1)})
    ze_results_cran <- c(ze_results_cran,r)
    r <- tryCatch(aze(test),error = function(e){return(1)})
    aze_results_cran <- c(aze_results_cran,r)
    r <- tryCatch(zed(test),error = function(e){return(1)})
    zed_results_cran <- c(zed_results_cran,r)
    r <- tryCatch(azed(test),error = function(e){return(1)})
    azed_results_cran <- c(azed_results_cran,r)
  }
  for (test in sl_data)
  {
    r <- tryCatch(findfrequency(test),error = function(e){return(1)})
    findfrequency_results_sl <- c(findfrequency_results_sl,r)
    r <- tryCatch(callSeasonLength(test),error = function(e){return(1)})
    seasonLength_results_sl  <- c(seasonLength_results_sl,r)
    r <- tryCatch(sazed(test,method = 'down'),error = function(e){return(1)})
    sazed_down_results_sl <- c(sazed_down_results_sl,r)
    r <- tryCatch(sazed(test,method='diff'),error = function(e){return(1)})
    sazed_diff_results_sl <- c(sazed_diff_results_sl,r)
    r <- tryCatch(sazed(test,method = 'alt'),error = function(e){return(1)})
    sazed_alt_results_sl <- c(sazed_alt_results_sl,r)
    r <- tryCatch(S(test),error = function(e){return(1)})
    s_results_sl <- c(s_results_sl,r)
    r <- tryCatch(Sa(test),error = function(e){return(1)})
    sa_results_sl <- c(sa_results_sl,r)
    r <- tryCatch(ze(test),error = function(e){return(1)})
    ze_results_sl <- c(ze_results_sl,r)
    r <- tryCatch(aze(test),error = function(e){return(1)})
    aze_results_sl <- c(aze_results_sl,r)
    r <- tryCatch(zed(test),error = function(e){return(1)})
    zed_results_sl <- c(zed_results_sl,r)
    r <- tryCatch(azed(test),error = function(e){return(1)})
    azed_results_sl <- c(azed_results_sl,r)
  }
  
  #create table
  result_table <- matrix(nrow = 11, ncol = 4)
  result_table[1,1] <- determineTolerance(findfrequency_results_cran,cran_expected,0.0)
  result_table[1,2] <- determineTolerance(findfrequency_results_cran,cran_expected,0.2)
  result_table[1,3] <- determineTolerance(findfrequency_results_sl,sl_expected,0.0)
  result_table[1,4] <- determineTolerance(findfrequency_results_sl,sl_expected,0.2)
  
  result_table[2,1] <- determineTolerance(seasonLength_results_cran,cran_expected,0.0)
  result_table[2,2] <- determineTolerance(seasonLength_results_cran,cran_expected,0.2)
  result_table[2,3] <- determineTolerance(seasonLength_results_sl,sl_expected,0.0)
  result_table[2,4] <- determineTolerance(seasonLength_results_sl,sl_expected,0.2)
  
  result_table[3,1] <- determineTolerance(sazed_down_results_cran,cran_expected,0.0)
  result_table[3,2] <- determineTolerance(sazed_down_results_cran,cran_expected,0.2)
  result_table[3,3] <- determineTolerance(sazed_down_results_sl,sl_expected,0.0)
  result_table[3,4] <- determineTolerance(sazed_down_results_sl,sl_expected,0.2)
  
  result_table[4,1] <- determineTolerance(sazed_diff_results_cran,cran_expected,0.0)
  result_table[4,2] <- determineTolerance(sazed_diff_results_cran,cran_expected,0.2)
  result_table[4,3] <- determineTolerance(sazed_diff_results_sl,sl_expected,0.0)
  result_table[4,4] <- determineTolerance(sazed_diff_results_sl,sl_expected,0.2)
  
  result_table[5,1] <- determineTolerance(sazed_alt_results_cran,cran_expected,0.0)
  result_table[5,2] <- determineTolerance(sazed_alt_results_cran,cran_expected,0.2)
  result_table[5,3] <- determineTolerance(sazed_alt_results_sl,sl_expected,0.0)
  result_table[5,4] <- determineTolerance(sazed_alt_results_sl,sl_expected,0.2)
  
  result_table[6,1] <- determineTolerance(s_results_cran,cran_expected,0.0)
  result_table[6,2] <- determineTolerance(s_results_cran,cran_expected,0.2)
  result_table[6,3] <- determineTolerance(s_results_sl,sl_expected,0.0)
  result_table[6,4] <- determineTolerance(s_results_sl,sl_expected,0.2)
  
  result_table[7,1] <- determineTolerance(sa_results_cran,cran_expected,0.0)
  result_table[7,2] <- determineTolerance(sa_results_cran,cran_expected,0.2)
  result_table[7,3] <- determineTolerance(sa_results_sl,sl_expected,0.0)
  result_table[7,4] <- determineTolerance(sa_results_sl,sl_expected,0.2)
  
  result_table[8,1] <- determineTolerance(ze_results_cran,cran_expected,0.0)
  result_table[8,2] <- determineTolerance(ze_results_cran,cran_expected,0.2)
  result_table[8,3] <- determineTolerance(ze_results_sl,sl_expected,0.0)
  result_table[8,4] <- determineTolerance(ze_results_sl,sl_expected,0.2)
  
  result_table[9,1] <- determineTolerance(aze_results_cran,cran_expected,0.0)
  result_table[9,2] <- determineTolerance(aze_results_cran,cran_expected,0.2)
  result_table[9,3] <- determineTolerance(aze_results_sl,sl_expected,0.0)
  result_table[9,4] <- determineTolerance(aze_results_sl,sl_expected,0.2)
  
  result_table[10,1] <- determineTolerance(zed_results_cran,cran_expected,0.0)
  result_table[10,2] <- determineTolerance(zed_results_cran,cran_expected,0.2)
  result_table[10,3] <- determineTolerance(zed_results_sl,sl_expected,0.0)
  result_table[10,4] <- determineTolerance(zed_results_sl,sl_expected,0.2)
  
  result_table[11,1] <- determineTolerance(azed_results_cran,cran_expected,0.0)
  result_table[11,2] <- determineTolerance(azed_results_cran,cran_expected,0.2)
  result_table[11,3] <- determineTolerance(azed_results_sl,sl_expected,0.0)
  result_table[11,4] <- determineTolerance(azed_results_sl,sl_expected,0.2)
  
  print(result_table)
  
  #CD-plot
  cran_expected <- unlist(cran_expected)
  
  dev.new()
  plotCD(data.frame(
    findFrequency=determineDistance(findfrequency_results_cran, cran_expected),
    seasonLength=determineDistance(seasonLength_results_cran, cran_expected),
    s=determineDistance(s_results_cran, cran_expected),
    sa=determineDistance(sa_results_cran, cran_expected),
    ze=determineDistance(ze_results_cran, cran_expected),
    zed=determineDistance(zed_results_cran, cran_expected),
    azed=determineDistance(azed_results_cran, cran_expected),
    aze=determineDistance(aze_results_cran, cran_expected),
    sazed_down=determineDistance(sazed_down_results_cran, cran_expected),
    sazed_diff=determineDistance(sazed_diff_results_cran, cran_expected),
    sazed_alt=determineDistance(sazed_alt_results_cran, cran_expected)
  ),
  
  decreasing=F,cex = 1)
  
  
  dev.new()
  plotCD(data.frame(
    findFrequency=determineDistance(findfrequency_results_sl, sl_expected),
    seasonLength=determineDistance(seasonLength_results_sl, sl_expected),
    s=determineDistance(s_results_sl, sl_expected),
    sa=determineDistance(sa_results_sl, sl_expected),
    ze=determineDistance(ze_results_sl, sl_expected),
    zed=determineDistance(zed_results_sl, sl_expected),
    azed=determineDistance(azed_results_sl, sl_expected),
    aze=determineDistance(aze_results_sl, sl_expected),
    sazed_down=determineDistance(sazed_down_results_sl, sl_expected),
    sazed_diff=determineDistance(sazed_diff_results_sl, sl_expected),
    sazed_alt=determineDistance(sazed_alt_results_sl, sl_expected)
  ),
  
  decreasing=F,cex = 1)
  
}

determineDistance <- function(a,b)
{
  return(mapply(function(x,y){
    if (length(y)==1)
    {
      return(abs(x-y))
    }
    return(x-y[which.min(abs(x-y))])
  },a,b))
}
determineTolerance <- function(a,b,tolerance)
{
  return(length(which(mapply(function(x,y){
    if (length(y)==1)
    {
      return(x >= y*(1-tolerance) & x <= y*(1+tolerance))
    }
    closest <- y[which.min(abs(x-y))]
    return(x >= closest*(1-tolerance) & x <= closest*(1+tolerance))
  },a,b))))
}
