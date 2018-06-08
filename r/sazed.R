sazed1 = function(y)
{
  library(signal);
  library(forecast);
  library(pracma);
  library(bspec);
  if (!is.ts(y))
  {
    print('y must be a time series')
    return(1);
  }
  #if (frequency(y) != 1)
  #{
  #  print('y already has a frequency');
  #  return(frequency(y));
  #}
  if (anyNA(y))
  {
    print('NA is not supported')
    return(1);
  }
  if (any(is.infinite(y)) || any(is.complex(y)))
  {
    print('complex numbers are not supported');
    return(1);
  }
  if (var(y) == 0)
  {
    print('y has no variance');
    return(1);
  }
  y = as.ts(detrend(as.vector(y)));
  y = scale(y);
  
  autocorrelation = computeAcf(y);
  ##S
  n = length(autocorrelation);
  periodigram = spec.pgram(autocorrelation,detrend=FALSE,plot = FALSE);
  
  if (round(n*2/pi) >= 4)
  {
    welch = welchPSD(as.ts(autocorrelation),round(n/(pi/2)));
  }
  else
  {
    welch = welchPSD(as.ts(autocorrelation),n);
  }
  period1 = round(1/(periodigram$freq[which.max(periodigram$spec)]));
  period2 = round(1/(welch$frequency[which.max(welch$power)]));
  
  if (min(period1,period2) > 2/pi * max(period1,period2))
  {
    return(round((period1+period2)/2));
  }
  ##spectral density estimation failed, compute and apply filter
  cutoff = (max(periodigram$spec))/n;
  cutoff = (cutoff + max(welch$power)/n)*pi*2;
  cutoff = min(c(1,cutoff));
  cutoff = max(c(0.05,cutoff));
  if (cutoff < 1)
  {
    b = butter(3,cutoff);
    y = filter(b$b,b$a,y);
    #filter dephases?!
  }
  
  ##azed
  return(azed(y,FALSE));
  
}
Sa = function(y,preprocess=TRUE)
{
  if (preprocess)
    y = preprocessTs(y);
  return(S(computeAcf(y),FALSE));
}
S = function(y,preprocess=TRUE)
{
  if (length(y) < 6)
  {
    return(1);
  }
  if (preprocess)
  {
    y = preprocessTs(y);
  }
  n = length(y);
  periodigram = spec.pgram(y,detrend=FALSE,plot = FALSE);
  if (n >= 6)
  {
    welch = welchPSD(as.ts(y),round(n*2/pi),windowfun = welchwindow);
  }
  else
  {
    welch = welchPSD(as.ts(y),n,windowfun = welchwindow);
  }
  #multitap = spec.mtm(as.ts(y),plot = FALSE)
  period1 = round(1/(periodigram$freq[which.max(periodigram$spec)]))
  period2 = round(1/(welch$frequency[which.max(welch$power)]))
  #period2 = round(1/multitap$freq[which.max(multitap$spec)])
  
  return(round((period1+period2)/2));
  
}
azed = function(y,preprocess=TRUE)
{
  if (preprocess)
    y = preprocessTs(y);
  return(zed(computeAcf(y),FALSE));
}
zed = function(y,preprocess=TRUE)
{
  if (preprocess)
  {
    y = preprocessTs(y);
  }
  n = length(y);
  signs = y;
  signs[which(signs < 0)] = -1;
  signs[which(signs > 0)] = 1;
  
  zero_distance_raw = which(signs[2:n] + signs[1:(n-1)] < 2 & signs[2:n] + signs[1:(n-1)] > -2);
  interpolation = y[zero_distance_raw] / (-1*(y[zero_distance_raw+1]-y[zero_distance_raw]));
  zero_distance_exact = zero_distance_raw+interpolation;
  z_to_z = diff(zero_distance_exact)
  #print(z_to_z)
  if (length(z_to_z) < 2)
  {
    return(1);
  }
  dens = density(z_to_z,kernel = 'epanechnikov')
  #zero_matrix = matrix(0,length(z_to_z),2);
  #zero_matrix[,1] = z_to_z;
  #plot(dens,main = "",xlab = 'Zero distances',ylab = 'Density');
  #print(z_to_z);
  if (!is.list(dens) && is.nan(dens))
  {
    return(1);
  }
  #plot(dens)

  return(round(dens$x[which.max(dens$y)]*2));
  #return(round(median(dens$x))*2)
}

aze = function(y,preprocess=TRUE)
{
  if (preprocess)
    y = preprocessTs(y);
  return(ze(computeAcf(y),FALSE));
}

ze = function(y,preprocess=TRUE)
{
  if (preprocess)
  {
    y = preprocessTs(y);
  }
  
  n = length(y);
  signs = y;
  signs[which(signs < 0)] = -1;
  signs[which(signs > 0)] = 1;
  
  zero_distance_raw = which(signs[2:n] + signs[1:(n-1)] < 2 & signs[2:n] + signs[1:(n-1)] > -2);
  interpolation = y[zero_distance_raw] / (-1*(y[zero_distance_raw+1]-y[zero_distance_raw]));
  zeros = zero_distance_raw+interpolation;
  
  return(round(mean(diff(zeros)))*2);
}
computeAcf = function(y)
{
  autocorrelation = as.ts(acf.fft(y));
  autocorrelation = autocorrelation[2:length(autocorrelation)];
  
  #this is optional and might improve the results,
  factor = 2/3;
  return(autocorrelation[1:round(length(autocorrelation)*factor)]);
}
preprocessTs = function(y)
{
  return(scale(as.ts(detrend(as.vector(y)))));
}
fazed = function(y)
{
  b = butter(3,0.8);
  y = filter(b$b,b$a,y);
  return(azed(y));
}
sfs = function(y, preprocess=TRUE, depth = 1)
{
  if (depth > 5)
  {
    return(1);
  }
    
  if (preprocess)
  {
    y = preprocessTs(y);
  }
  autocorrelation = computeAcf(y);
  ##S
  n = length(autocorrelation);
  periodigram = spec.pgram(autocorrelation,detrend=FALSE,plot = FALSE);
  
  if (round(n*2/pi) >= 4)
  {
    welch = welchPSD(as.ts(autocorrelation),round(n/(pi/2)));
  }
  else
  {
    welch = welchPSD(as.ts(autocorrelation),n);
  }
  period1 = round(1/(periodigram$freq[which.max(periodigram$spec)]));
  period2 = round(1/(welch$frequency[which.max(welch$power)]));
  
  if (min(period1,period2) > 2/pi * max(period1,period2))
  {
    return(round((period1+period2)/2));
  }
  ##spectral density estimation failed, compute and apply filter
  cutoff = (max(periodigram$spec))/n;
  cutoff = (cutoff + max(welch$power)/n)*2*pi;
  cutoff = min(c(1,cutoff));
  cutoff = max(c(0.05,cutoff));
  if (cutoff < 1)
  {
    b = butter(2,cutoff);
    y = filter(b$b,b$a,y);
    #filter dephases?!
  }
  return(sfs(y,FALSE,depth = depth+1));
}
sazed2 <- function(y,iter=0)
{
  library(signal);
  library(forecast);
  library(pracma);
  library(bspec);
  if (anyNA(y))
  {
    print('NA is not supported')
    return(1);
  }
  if (any(is.infinite(y)) || any(is.complex(y)))
  {
    print('complex numbers are not supported');
    return(1);
  }
  if (var(y) == 0)
  {
    print('y has no variance');
    return(1);
  }

  results = c();
  results = c(results,S(y));
  results = c(results,Sa(y));
  results = c(results,ze(y));
  results = c(results,aze(y));
  results = c(results,zed(y));
  results = c(results,azed(y));
  #results = c(results,callSeasonLength(y));
  #results = c(results,results[6]);
  #results = c(results,sazed(y));
  #results = c(results,fazed(y));
  #results = c(results,sfs(y));
  #results = c(results,findfrequency(y));
  #results = c(results,results[7]);
  #results = c(results,find.freq(y));
  #sanity check?
  #means of seasons must equal mean of detrended time series
  #write list of constants
  print(results)
  results = results[which(!is.infinite(results))];
  results = results[which(!is.na(results))];
  if (var(results) == 0)
  {
    return(results[1]);
  }
  unique_results = unique(results);
  unique_results = unique_results[which(unique_results != 2)];
  unique_results = unique_results[which(!is.infinite(unique_results))];
  tab = tabulate(match(results,unique_results));
  if (max(tab) == 1)
  {
    iter = iter + 1;
    #dens = density(results,kernel = 'epanechnikov')
    if (mod(iter,2) == 1)
    {
      return(2*sazed2(downsample(y,2),iter));
    }
    return(sazed2(diff(y,lag = 2),iter));
    #return(round(dens$x[which.max(dens$y)]));
    #return(sazed2(2*downsample(diff(y,lag = 2),2)));
    #return(2*sazed2(downsample(y,2)));
    #return(round(max(kmeans(results,c(min(results),max(results)))$centers)))
    
  }
  else if (length(tab[which(tab == max(tab))]) > 1)
  {
    #print(tab);
    #print(results);
    sorted = sort(unique_results,decreasing = TRUE)
    #print(sorted[which.max(tabulate(match(sort(results,decreasing = TRUE),sorted)))]);
    return(sorted[which.max(tabulate(match(sort(results,decreasing = TRUE),sorted)))]);
  }

  vote = unique_results[which.max(tabulate(match(results,unique_results)))];
  return(vote);
  
  #check1
  #for (result1 in unique_results)
  #{
  #  for (result2 in unique_results)
  #  {
  #    if (result1 == result2)
  #      next;
  #    if (mod(result1,result2) == 0 & vote==result2 & result1 < median(results))
  #    {
  #      cat(paste0("Result1: ",result1, " result2: ",result2, " mode: ",vote," median: ",median(results),"\n"));
  #      vote = result1;
  #    }
  #  }
  #}
}

downsample = function(data, window_size=2)
{
  library(zoo);
  return(ts(as.ts(rollapply(zoo(data),width=2,by=2,FUN=mean)),frequency=1))
  #n = length(data);
  #result = c();
  #i = 1;
  #while (i < n)
  #{
  #  result = c(result,mean(data[i:(i+window_size-1)]));
  #  i = i+window_size;
  #}
  #return(as.ts(result[which(!is.na(result))]));
}