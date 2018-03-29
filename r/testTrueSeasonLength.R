compareMethods = function(y)
{
  if (is.character(y))
  {
    y = read.table(y);
    if (length(y) == 2)
    {
      y = y[2];
    }
    else
    {
      y = y[1];
    }
    y = as.ts(y);
  }
  plot(y)
  z = detrend(as.vector(y));
  z = scale(z);
  z = acf.fft(z);
  z = z[2:length(z)];
  #z = z[1:round(length(z)*2/3)];
  #plot(z,type='l');
  source('trueSeasonLength.R');
  source('sazed.R')
  cat(paste0('Expected -> ',frequency(y),'\n'));
  #cat(paste0('Method1: periodigram -> ',spectral_analysis(y),'\n'));
  #z = findAndFilter(z);
  #cat(paste0('Method2: symbol ->',symbol(z),'\n'));
  #cat(paste0('Method3: azed->',azed(y),'\n'));
  #cat(paste0('Method4: peak->',peak(z),'\n'));
  cat(paste0('Method1: S->',S(y),'\n'));
  cat(paste0('Method2: Sa->',Sa(y),'\n'));
  cat(paste0('Method3: ze->',ze(y),'\n'));
  cat(paste0('Method4: zed->',zed(y),'\n'));
  cat(paste0('Method5: azed->',azed(y),'\n'));
  cat(paste0('Method6: sazed->',sazed(y),'\n'));
  cat(paste0('Method7: fazed->',fazed(y),'\n'));
  cat(paste0('Method8: sfs->',sfs(y),'\n'));
  cat(paste0('Method9: findfrequency->',findfrequency(y),'\n'));
  cat(paste0('Method10: find.freq->',find.freq(y),'\n'));

}
testSeasonLength = function(datafile){
  #print(datafile);
  source('seasonLength.R');
  y = read.table(datafile);
  if (length(y) == 2)
  {
    y = y[2];
  }
  else
  {
    y = y[1];
  }  
  y = as.ts(y);
  return(seasonLength(y));
}


testTrueSeasonLength = function(datafile){
  #print(datafile);
  source('trueSeasonLength.R');
  if (is.character(datafile))
  {
    y = read.table(datafile);
    if (length(y) == 2)
    {
      y = y[2];
    }
    else
    {
      y = y[1];
    }  
    y = as.ts(y);
  }
  else
  {
    y = as.ts(as.vector(datafile));
  } 
  
  return(trueSeasonLength(y));
}
testSazed = function(datafile){
  #print(datafile);
  source('sazed.R');
  y = read.table(datafile);
  if (length(y) == 2)
  {
    y = y[2];
  }
  else
  {
    y = y[1];
  }  
  y = as.ts(y);
  return(sazed(y));
}


testPublicDatasets <- function()
{
  
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
  findfreq_basic_results <- c()
  findfreq_advanced_results <- c()
  s_results = c()
  sa_results = c()
  ze_results = c()
  zed_results = c()
  azed_results = c()
  fazed_results = c()
  sazed_results = c()
  sfs_results = c()
  ensemble_results = c()
  #seasonLength_results = c()
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
      findfreq_basic_results <- c(findfreq_basic_results, find.freq(ts_data_nofreq))
      findfreq_advanced_results <- c(findfreq_advanced_results, findfrequency(ts_data_nofreq))
      s_results = c(s_results,S(ts_data_nofreq))
      sa_results = c(sa_results,Sa(ts_data_nofreq))
      ze_results = c(ze_results,ze(ts_data_nofreq))
      zed_results = c(zed_results,zed(ts_data_nofreq))
      azed_results = c(azed_results,azed(ts_data_nofreq))
      fazed_results = c(fazed_results,fazed(ts_data_nofreq))
      sazed_results = c(sazed_results,sazed(ts_data_nofreq))
      sfs_results = c(sfs_results,sfs(ts_data_nofreq))
      ensemble_results = c(ensemble_results,ensemble(ts_data_nofreq))
      #seasonLength_results = c(seasonLength_results,callSeasonLength(ts_data_nofreq))
    }
  }
  #Accuracy
  cat(paste0('find.freq passes ', length(which(findfreq_basic_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('findfrequency passes ', length(which(findfreq_advanced_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('s passes ', length(which(s_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('sa passes ', length(which(sa_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('ze passes ', length(which(ze_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('zed passes ', length(which(zed_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('azed passes ', length(which(azed_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('fazed passes ', length(which(fazed_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('sazed passes ', length(which(sazed_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('sfs passes ', length(which(sfs_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  cat(paste0('ensemble passes ', length(which(ensemble_results == actual_results)), 
             ' out of ', number_public_datasets, '\n'))
  #cat(paste0('seasonLength passes ', length(which(seasonLength_results == actual_results)), 
  #           ' out of ', number_public_datasets, '\n'))
  #Friedman Rank Test - best algorithm is the one with lowest distance to actual_results (hence decreasing=F)
  #if (!require("devtools")) {
  #  install.packages("devtools")
  #}
  #devtools::install_github("b0rxa/scmamp")
  library(scmamp)
  plotCD(data.frame(hynd_basic=abs(findfreq_basic_results - actual_results), 
                    hynd_advanc=abs(findfreq_advanced_results - actual_results),
                    s=abs(s_results - actual_results),
                    sa=abs(sa_results - actual_results),
                    ze=abs(ze_results - actual_results),
                    zed=abs(zed_results - actual_results),
                    azed=abs(azed_results - actual_results),
                    fazed=abs(fazed_results - actual_results),
                    sazed=abs(sazed_results - actual_results),
                    sfs=abs(sfs_results - actual_results),
                    ensemble=abs(ensemble_results - actual_results)
                    #seasonLength=abs(seasonLength_results - actual_results)
                    ),
                    
         decreasing=F)
}

testAll = function()
{
  sum = 0;
  expected = c(88,57,3,6,5,25,214,900,10,1000,20,15,190,107,52,1000,20,35,50000,3);
  sum = sum + testSuite(expected,"../test1/t",1);
  
  expected = c(10000,10,96,120,19,180,1332,45,19,20,140,18,30,1500,740,20,1590,10,29,5);
  sum = sum + testSuite(expected,"../test2/t",2);
  
  expected = matrix(
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
    ,nrow = 3,ncol = 20);
  expected = t(expected);
  sum = sum + testMultiple(expected,"../test3/t",3);
  
  expected = c(128,128,128,128,128,156,156,156,156,156,37,37,37,37,37,8,8,8,8,8);
  sum = sum + testSuite(expected,"../test4/t",4);

  expected = c(155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155);
  sum = sum + testSuite(expected,"../test5/t",5);
  
  expected = c(4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32786,65536);
  sum = sum + testSuite(expected,"../test6/t",6);

  expected = c(1,1,1,1,1,1,1,1,1,1)
  sum = sum + testSuite(expected,"../test7/t",7);
  
  expected = c(12,12,12,4,12,6,6,6,12,6,4,12,12,12,12,12,3,3,12,12)
  sum = sum + testSuite(expected,"../testR/t",8);
  
  expected = c(12,12,12,12,12,12,12,130,12,12,12,12,12,12,130,12,12,12,12,12)
  sum = sum + testSuite(expected,"../testC/t",9);
  
  print(paste('Total: ', sum, ' out of 165',sep=''));
  print('original seasonLength passes 122');
  print('findFrequency passes 83');
}


testSuite = function(expected,name,number)
{
  n = length(expected);
  passed = 0;
  results = c(1:n);
  bin = c(1:n);
  print(paste("Test suite ",number,":",sep=''));
  for (i in 1:n)
  {
    result = testSazed(paste(name,i,sep='',collapse = ' '));
    results[i] = result;
    if ((result >= (expected[i]*0.8))  && (result <= (expected[i]*1.2)))
    {
      passed = passed + 1;
      bin[i] = 1;
    }
    else
    {
      bin[i] = 0;
    }
  }
  
  print(paste(passed," out of ",n,sep='',collapse = ' '));
  print(t(matrix(c(bin,results,expected),nrow=n,ncol=3)));
  #print(paste(c("Passed:   ",bin),sep='',collapse = ' '));
  #print(paste(c("Results:  ",results),sep='',collapse = ' '));
  #print(paste(c("Expected: ", expected),sep = '',collapse = ' '));
  print("");
  return(passed);
}

testMultiple = function(expected,name,number)
{
  n = nrow(expected);
  passed = 0;
  results = c(1:n);
  bin = c(1:n);
  solution = c(1:n);
  print(paste("Test suite ",number,":",sep=''));
  for (i in 1:n)
  {
    result = testSazed(paste(name,i,sep='',collapse = ' '));
    results[i] = result;
    for (j in 1:ncol(expected))
    {
      if ((result >= (expected[i,j]*0.8))  && (result <= (expected[i,j]*1.2)))
      {
        passed = passed + 1;
        bin[i] = 1;
        solution[i] = j;
        break;
      }
      else if (j == ncol(expected))
      {
        bin[i] = 0;
        solution[i] = 1;
      }
    }
  }
  solution = (cbind(c(1:n),c(solution)));
  print(paste(passed," out of ",n,sep='',collapse = ' '));
  print(t(matrix(c(bin,results,expected[solution]),nrow=n,ncol=3)));
  #print(paste(c("Passed:   ",bin),sep='',collapse = ' '));
  #print(paste(c("Results:  ",results),sep='',collapse = ' '));
  #print(paste(c("Expected: ", expected),sep = '',collapse = ' '));
  print("");
  return(passed);
}

testNew = function(y)
{
  source('newSeasonLength.R')
  if (is.character(y))
  {
    y = read.table(y);
    if (length(y) == 2)
    {
      y = y[2];
    }
    else
    {
      y = y[1];
    }
    y = as.ts(y);
  }
  
  cat(paste0('Expected: ', frequency(y),'\n'));
  plot(y)
  return(newSeasonLength(y));
}
