sazed = function(y)
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
  if (preprocess)
  {
    y = preprocessTs(y);
  }
  n = length(y);
  periodigram = spec.pgram(y,detrend=FALSE,plot = FALSE);
  
  if (round(n*2/pi) >= 4)
  {
    welch = welchPSD(as.ts(y),round(n*2/pi));
  }
  else
  {
    welch = welchPSD(as.ts(y),n);
  }
  period1 = round(1/(periodigram$freq[which.max(periodigram$spec)]))
  period2 = round(1/(welch$frequency[which.max(welch$power)]))
  
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
  
  zero_distance_raw = which(signs[2:n] + signs[1:(n-1)] == 0);
  interpolation = y[zero_distance_raw] / (-1*(y[zero_distance_raw+1]-y[zero_distance_raw]));
  zero_distance_exact = zero_distance_raw+interpolation;
  z_to_z = sort(zero_distance_exact[2:n]-zero_distance_exact[1:(n-1)]);
  
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
  return(round(dens$x[which.max(dens$y)]*2));
  #return(round(median(dens$x))*2)
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
  
  zero_distance_raw = which(signs[2:n] + signs[1:(n-1)] == 0);
  interpolation = y[zero_distance_raw] / (-1*(y[zero_distance_raw+1]-y[zero_distance_raw]));
  zeros = zero_distance_raw+interpolation;
  return(round(mean(diff(round(zeros))))*2);
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
ensemble = function(y)
{
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
  results = c(results,zed(y));
  results = c(results,azed(y));
  results = c(results,sazed(y));
  results = c(results,fazed(y));
  results = c(results,sfs(y));
  results = c(results,findfrequency(y));
  results = c(results,find.freq(y));
  #sanity check?
  #means of seasons must equal mean of detrended time series
  #write list of constants
  unique_results = unique(results);
  return(unique_results[which.max(tabulate(match(results,unique_results)))]);
}