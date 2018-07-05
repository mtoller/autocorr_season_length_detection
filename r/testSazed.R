#' Test Function for SAZED
#' 
#' \code{testSazed} loads a time series from a file and passes it as input to \code{sazed}.
#'
#' @param datafile system path to the file to be loaded.
#' @return the season length estimated by \code{sazed}.
testSazed <- function(datafile){
  #print(datafile)
  source('sazed.R')
  if (is.character(datafile))
  {
    y <- read.table(datafile)
    if (length(y) == 2)
    {
      y <- y[2]
    }
    else
    {
      y <- y[1]
    }
    y <- as.ts(y)
  }
  else
  {
    y <- as.ts(as.vector(datafile))
  }
  
  r = tryCatch(
    {
      r = sazed(y);
    },
    error = function(cond)
    {
      r = 1;
    }
  )
  return(r);
}
#' Testing Public Datasets
#' 
#' \code{testPublicDatasets} loads time series from various R libraries and uses them to test
#' different season length detectors.
#'
#' @param tolerance the allowed deviation between expected and actual result.
#' @param testExternal if true, an external octave function will be tested as well.
#' @examples 
#' testPublicDatasets()
#' testPublicDatasets(0.2)
testPublicDatasets <- function(tolerance = 0,testExternal=F)
{
  require(signal)
  require(forecast)
  require(pracma)
  require(bspec)
  public_datasets <- list()
  public_datasets[["fma"]] <- c("airpass", "beer", "bricksq", "condmilk", "elec", "fancy", "hsales",
                                "hsales2", "labour", "milk", "motion", "plastics", "qsales", "ukdeaths",
                                "usdeaths", "uselec", "writing")
  public_datasets[["expsmooth"]] <- c("cangas", "enplanements", "frexport", "mcopper", "ukcars", "utility", 
                                      "vehicles", "visitors")
  public_datasets[["fpp2"]] <- c("a10", "ausbeer", "auscafe", "austourists", "debitcards", "elecequip", "h02",
                                 "hyndsight", "qauselec", "qcement", "qgas", "usmelec")
  public_datasets[["TSA"]] <- c("airmiles", "co2", "flow", "JJ", "oilfilters", "retail", "tempdub")
  public_datasets[["astsa"]] <- c("birth", "cmort", "flu", "gas", "oil", "part", "prodn", "rec",
                                  "so2", "soi", "sunspotz", "tempr", "UnempRate")
  public_datasets[["AER"]] <- c("DutchSales", "UKNonDurables")
  number_public_datasets <- Reduce(sum, lapply(public_datasets, length))
  
  actual_results <- c()
  findfreq_advanced_results <- c()
  seasonLength_results  <- c()
  s_results <- c()
  sa_results <- c()
  ze_results <- c()
  zed_results <- c()
  azed_results <- c()
  aze_results <- c()
  sazed_down_results <- c()
  sazed_diff_results <- c()
  sazed_alt_results <- c()
  #install.packages(c("fma", "expsmooth", "fpp2", "TSA", "astsa", "AER"))
  source('trueSeasonLength.R')
  source('baselines.R')
  source('newSeasonLength.R')
  source('sazed.R')
  dataset_libraries <- c("fma", "expsmooth", "fpp2", "TSA", "astsa", "AER")
  for (a_dataset_library in dataset_libraries) {
    library(a_dataset_library, character.only = T)
    data(list = public_datasets[[a_dataset_library]])
    for (a_dataset in public_datasets[[a_dataset_library]]) {
      ts_data <- get(a_dataset)
      ts_data_nofreq <- ts(as.vector(ts_data), frequency = 1)
      actual_results <- c(actual_results, frequency(ts_data))
      findfreq_advanced_results <- c(findfreq_advanced_results, findfrequency(ts_data_nofreq))
      if (testExternal)
      {
        seasonLength_results = c(seasonLength_results,callSeasonLength(ts_data_nofreq))
      }
      s_results = c(s_results,S(ts_data_nofreq))
      sa_results = c(sa_results,Sa(ts_data_nofreq))
      ze_results = c(ze_results,ze(ts_data_nofreq))
      zed_results = c(zed_results,zed(ts_data_nofreq))
      azed_results = c(azed_results,azed(ts_data_nofreq))
      aze_results = c(aze_results,aze(ts_data_nofreq))
      sazed_down_results = c(sazed_down_results,sazed(ts_data_nofreq,method = "down"))
      sazed_diff_results = c(sazed_diff_results,sazed(ts_data_nofreq,method = "diff"))
      sazed_alt_results = c(sazed_alt_results,sazed(ts_data_nofreq,method = "alt"))
    }
  }
  #Accuracy
  cat(paste0('findfrequency passes ', evaluateTolerance(findfreq_advanced_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  if (testExternal)
  {
    cat(paste0('seasonLength passes ', evaluateTolerance(seasonLength_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  }
  cat(paste0('s passes ', evaluateTolerance(s_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('sa passes ', evaluateTolerance(sa_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('ze passes ', evaluateTolerance(ze_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('aze passes ', evaluateTolerance(aze_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('zed passes ', evaluateTolerance(zed_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('azed passes ', evaluateTolerance(azed_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('sazed_down passes ', evaluateTolerance(sazed_down_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('sazed_diff passes ', evaluateTolerance(sazed_diff_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('sazed_alt passes ', evaluateTolerance(sazed_alt_results,actual_results,tolerance), 
             ' out of ', number_public_datasets, '\n'))
  #Friedman Rank Test - best algorithm is the one with lowest distance to actual_results (hence decreasing=F)
  #if (!require("devtools")) {
  #  install.packages("devtools")
  #}
  #devtools::install_github("b0rxa/scmamp")
  library(scmamp)
  if (testExternal)
  {
    plotCD(data.frame(#hynd_basic=abs(findfreq_basic_results - actual_results), 
                    findFrequency=abs(findfreq_advanced_results - actual_results),
                    seasonLength=abs(seasonLength_results - actual_results),
                    s=abs(s_results - actual_results),
                    sa=abs(sa_results - actual_results),
                    ze=abs(ze_results - actual_results),
                    zed=abs(zed_results - actual_results),
                    azed=abs(azed_results - actual_results),
                    aze=abs(aze_results - actual_results),
                    sazed_down=abs(sazed_down_results - actual_results),
                    sazed_diff=abs(sazed_diff_results - actual_results),
                    sazed_alt=abs(sazed_alt_results - actual_results)
                    ),
                    
         decreasing=F,cex = 1)
  }
  else
  {
    plotCD(data.frame(#hynd_basic=abs(findfreq_basic_results - actual_results), 
      findFrequency=abs(findfreq_advanced_results - actual_results),
      s=abs(s_results - actual_results),
      sa=abs(sa_results - actual_results),
      ze=abs(ze_results - actual_results),
      zed=abs(zed_results - actual_results),
      azed=abs(azed_results - actual_results),
      aze=abs(aze_results - actual_results),
      sazed_down=abs(sazed_down_results - actual_results),
      sazed_diff=abs(sazed_diff_results - actual_results),
      sazed_alt=abs(sazed_alt_results - actual_results)
    ),
    
    decreasing=F,cex = 1)
  }
}
#' Testing the Diverse Dataset
#' 
#' \code{testDiverse} loads predefined local time series and tests its argument on them.
#'
#' @param func the function for reading and computing the season length of a time series.
#' @param deviation the allowed deviation between expected and actual result.
#' @examples 
#' testDiverse(testSazed,0)
#' testDiverse(callSeasonLength,0.2)
testDiverse <- function(func,deviation)
{
  
  sum <- 0
  expected <- c(88,57,3,6,5,25,214,900,10,1000,20,15,190,107,52,1000,20,35,50000,3)
  sum <- sum + testSuite(func,expected,"../test1/t",1,deviation)
  
  expected <- c(10000,10,96,120,19,180,1332,45,19,20,140,18,30,1500,740,20,1590,10,29,5)
  sum <- sum + testSuite(func,expected,"../test2/t",2,deviation)
  
  expected <- matrix(
    c(20,60,-1,
      3,10,-1,
      50,90,-1,
      60,300,-1,
      70,700,-1,
      40,3498,-1,
      8,32,-1,
      40,75,-1,
      400,1000,-1,
      9,20,-1,
      5,980,-1,
      10,64,-1,
      78,100,-1,
      6,12,18,
      3,7,-1,
      9,24,-1,
      7,70,700,
      400,800,-1,
      6,20,-1,
      500,1400,2000)
    ,nrow = 3,ncol = 20)
  expected <- t(expected)
  sum <- sum + testMultiple(func,expected,"../test3/t",3,deviation)
  
  expected <- c(128,128,128,128,128,156,156,156,156,156,37,37,37,37,37,8,8,8,8,8)
  sum <- sum + testSuite(func,expected,"../test4/t",4,deviation)

  expected <- c(155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155)
  sum <- sum + testSuite(func,expected,"../test5/t",5,deviation)
  
  expected <- c(4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536)
  sum <- sum + testSuite(func,expected,"../test6/t",6,deviation)

  expected <- c(1,1,1,1,1,1,1,1,1,1)
  sum <- sum + testSuite(func,expected,"../test7/t",7,deviation)
  
  expected <- c(12,12,12,4,12,6,6,6,12,6,4,12,12,12,12,12,3,3,12,12)
  sum <- sum + testSuite(func,expected,"../testR/t",8,deviation)
  
  expected <- c(12,12,12,12,12,12,12,130,12,12,12,12,12,12,130,12,12,12,12,12)
  sum <- sum + testSuite(func,expected,"../testC/t",9,deviation)
  
  print(paste('Total: ', sum, '/165, ',sum/165,sep=''))
}

#' Testing a Test Suite
#' 
#' \code{testSuite} is an internal function used inside testDiverse.
#'
#' @param func the function to be tested.
#' @param expected a vector of expected results.
#' @param name path to the test folder.
#' @param number number of the test suite (used for printing only)
#' @param deviation the allowed deviation between expected and actual result.
#' @return the number of passed tests.
#' @examples 
#' expected <- c(88,57,3,6,5,25,214,900,10,1000,20,15,190,107,52,1000,20,35,50000,3)
#' passed <- testSuite(func,expected,"../test1/t",1,deviation)
testSuite <- function(func,expected,name,number,deviation)
{
  n <- length(expected)
  passed <- 0
  results <- c(1:n)
  bin <- c(1:n)
  print(paste("Test suite ",number,":",sep=''))
  for (i in 1:n)
  {
    result <- func(paste(name,i,sep='',collapse = ' '))
    results[i] <- result
    if ((result >= (expected[i]*(1-deviation))  && (result <= (expected[i]*(1+deviation)))))
    {
      passed <- passed + 1
      bin[i] <- 1
    }
    else
    {
      bin[i] <- 0
    }
  }
  
  print(paste(passed,"/",n,sep='',collapse = ' '))
  print(t(matrix(c(bin,results,expected),nrow=n,ncol=3)))
  print("")
  return(passed)
}
#' Testing a Test Suite with Multiple Possible Results
#' 
#' \code{testMultiple} is an internal function used inside testDiverse.
#'
#' @param func the function to be tested.
#' @param expected a matrix of expected results, each column containing the possible results for a
#' test
#' @param name path to the test folder.
#' @param number number of the test suite (used for printing only)
#' @param deviation the allowed deviation between expected and actual result.
#' @return the number of passed tests
#' @examples 
#' expected <- matrix(
#' c(20,60,-1,
#'  3,10,-1,
#'  50,90,-1,
#'  60,300,-1,
#'  70,700,-1,
#'  40,3498,-1,
#'  8,32,-1,
#'  40,75,-1,
#'  400,1000,-1,
#'  9,20,-1,
#'  5,980,-1,
#'  10,64,-1,
#'  78,100,-1,
#'  6,12,18,
#'  3,7,-1,
#'  9,24,-1,
#'  7,70,700,
#'  400,800,-1,
#'  6,20,-1,
#'  500,1400,2000)
#' ,nrow = 3,ncol = 20)
#' expected <- t(expected)
#' testMultiple(func,expected,"../test3/t",3,deviation)
testMultiple <- function(func,expected,name,number,deviation)
{
  n <- nrow(expected)
  passed <- 0
  results <- c(1:n)
  bin <- c(1:n)
  solution <- c(1:n)
  print(paste("Test suite ",number,":",sep=''))
  for (i in 1:n)
  {
    result <- func(paste(name,i,sep='',collapse = ' '))
    results[i] <- result
    for (j in 1:ncol(expected))
    {
      if ((result >= (expected[i,j]*(1-deviation)))  && (result <= (expected[i,j]*(1+deviation))))
      {
        passed <- passed + 1
        bin[i] <- 1
        solution[i] <- j
        break
      }
      else if (j == ncol(expected))
      {
        bin[i] <- 0
        solution[i] <- 1
      }
    }
  }
  solution <- (cbind(c(1:n),c(solution)))
  print(paste(passed,"/",n,sep='',collapse = ' '))
  print(t(matrix(c(bin,results,expected[solution]),nrow=n,ncol=3)))
  print("")
  return(passed)
}
#' Evalutate Tolerance Deviation
#' 
#' \code{evaluateTolerance} checks if its arguments are in a tolerance range about each other. 
#'
#' @param a the vector of real values
#' @param b the vector of expected values
#' @param tolerance the accepted interval around \code{b} in which \code{a} may be
#' @return the number of elements in \code{a} which are inside the tolerance interval around
#' \code{b}
#' @examples 
#' a <- c(1,10,4,5,8,0)
#' b <- c(1,9,12,6,7,-0.2)
#' evaluateTolerance(a,b,0.5)
evaluateTolerance <- function(a,b,tolerance)
{
  return(length(which(a >= b*(1-tolerance) & a <= b*(1+tolerance))))
}
