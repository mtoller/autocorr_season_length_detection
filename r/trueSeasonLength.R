skewness = function(x)
{
  n = length(x);
  m = mean(x);
  s = std(x);
  return(sum(((x-m)/s)^3)/n);
}

excess = function(x)
{
  n = length(x);
  m = mean(x);
  s = std(x);
  return(sum(((x-m)/s)^4)/n-3);
}

modus = function(x)
{
  unique_x = unique(x);
  return(unique_x[which.max(tabulate(match(x,unique_x)))]);
}
modusCount = function(x)
{
  return(max(tabulate(match(x,unique(x)))));
}

harmonic = function(n)
{
  return(cumsum(1/c(1:n)));
}

findAndFilter = function(y,m)
{
  b = butter(5,m);
  y = filter(b$b,b$a,y);
  return(y);
}

spectral_analysis = function(y)
{
  autocorrelation = acf.fft(y);
  autocorrelation = autocorrelation[2:length(autocorrelation)];
  autocorrelation = autocorrelation[1:round(length(autocorrelation)*2/3)];
  n = length(autocorrelation);
  periodigram = spec.pgram(autocorrelation,detrend=FALSE,plot = FALSE);
  
  if (round(n*2/pi) >= 4)
  {
    welch = welchPSD(as.ts(autocorrelation),round(n*2/pi));
  }
  else
  {
    welch = welchPSD(as.ts(autocorrelation),n);
  }
  
  
  
  cutoff = (max(periodigram$spec))/n;
  cutoff = (cutoff + max(welch$power)/n)*2*pi;
  cutoff = min(c(1,cutoff));
  cutoff = max(c(0.05,cutoff));
  #plot(p$freq,p$spec,type='l');
  #plot(w$frequency,w$power,type='l');
  #print(max(p$spec)/na);
  #print(max(w$power)/na);
  #print(round(1/(p$freq[which.max(p$spec)])));
  #print(round(1/(w$frequency[which.max(w$power)])));
  #print(m);
  period1 = round(1/(periodigram$freq[which.max(periodigram$spec)]))
  period2 = round(1/(welch$frequency[which.max(welch$power)]))
  if (min(period1,period2) > 0.8 * max(period1,period2))
  {
    return(c(cutoff,round((period1+period2)/2)))
  }
  return(c(cutoff,1));
}

azed = function(y, harmonic_analysis = FALSE, zero_discard=FALSE,skew_correction=FALSE)
{
  autocorrelation = acf.fft(y);
  autocorrelation = autocorrelation[2:length(autocorrelation)];
  autocorrelation = autocorrelation[1:round(length(autocorrelation)*2/3)];
  #plot(y);
  #plot(a,type = 'l');
  n = length(autocorrelation);
  if (harmonic_analysis == TRUE)
  {
    autocorrelation = autocorrelation * (1/harmonic(round(n*2/3)));
  }
  
  signs = autocorrelation;
  signs[which(signs < 0)] = -1;
  signs[which(signs > 0)] = 1;
  
  zero_distance_raw = which(signs[2:n] + signs[1:(n-1)] == 0);
  interpolation = autocorrelation[zero_distance_raw] / (-1*(autocorrelation[zero_distance_raw+1]-autocorrelation[zero_distance_raw]));
  zero_distance_exact = zero_distance_raw+interpolation;
  z_to_z = sort(zero_distance_exact[2:n]-zero_distance_exact[1:(n-1)]);
  
  if (length(z_to_z) < 2)
  {
    return(NaN);
  }
  
  if (zero_discard == TRUE)
  {
    z_to_z = z_to_z[which(z_to_z >= 3)];
  }
  
  dens = density(z_to_z,kernel = 'epanechnikov')
  if (skew_correction && skewness(dens$y) > 0)
  {
    dens$y = dens$y * c(1:round(length(dens$y)))/sum(1:length(dens$y))*length(dens$y);
  }
  zero_matrix = matrix(0,length(z_to_z),2);
  zero_matrix[,1] = z_to_z;
  
  
  #plot(zero_matrix);
  #plot(d);
  
  return(dens);
}

trueSeasonLength = function(y){
  library(signal);
  library(forecast);
  library(pracma);
  library(bspec);
  source('fastAcf.R')
  source('fastMa.R')
  if (!is.ts(y))
  {
    print('y must be a time series')
    return(1);
  }
  if (frequency(y) != 1)
  {
    print('y already has a frequency');
    return(frequency(y));
  }
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
  
  #if (length(y) < 22) {
  #  return(1);
  #}
  #par(mfrow=c(4,2));
  y = detrend(as.vector(y));
  y = scale(y);
  spectral_result = spectral_analysis(y);
  result1 = spectral_result[2];
  cutoff = spectral_result[1];
  if (result1 != 1)
  {
    return(result1);
  }
  if (cutoff < 1)
  {
    y = findAndFilter(y,cutoff);
  }
  dens = azed(y);
  if (!is.list(dens) && is.nan(dens))
  {
    return(1);
  }
  
  result2 = round(dens$x[which.max(dens$y)]*2);
  
  
  
  #-------------------------------------
  #ma = forecast::ma(a,s,centre = FALSE)
  #entropy and compression
  #kalman
  #forecasting/prediction
  #backproject with shapelet stauchen
  #implement old algorithm in R
  #anisotropic diffusion
  
  

  #print(mean(a));
  #print(std(a));
  #print(skewness(a));
  #print(excess(a));
  #ma = as.ts(fastMa(a,round(s/2)));
  #if (skewness(d$y) > 0)
  #{
  #  print(d$x[which.max(d$y)]-d$x[which.max(old_dy)]);
  #  print(max(abs(d$y-old_dy)));
  #  print(d$x[which.max(abs(d$y-old_dy))]);
  #}
  
  #s = 2*(round(median(z_to_z)));
  #print(d);
  #print(which.max(d$y));
  #print(sum(d$y));
  #print(d$x);

  #print(d$y*100);
  #interval_borders = unique(round((which(diff(sign(diff(d$y)))==-2)+1)/512*nrow(m)));
  #print(interval_borders);
  #print(z_to_z[interval_borders])
  
  #plot(m);
  #for(border in interval_borders)
  #{
  #  lines(c(z_to_z[border],z_to_z[border]),c(-1,1),col='red');
  #}
  ##noise measure, scales with noise AND short periodicities
  #t = a[2:(na)] - a[1:(na-1)];
  #t[which(t < 0)] = -1;
  #t[which(t > 0)] = 1;
  #inflections = which(t[2:length(t)] + t[1:(length(t)-1)]==0)+1;
  #d1 = diff(inflections);
  #d2 = c(diff(inflections,1,2),1);
  #w = c((na-1):1)/sum(1:(na-1))*(na-1);
  #z = which(d2 == 0);
  #count = sum(d1[z]*w[inflections[z+1]]);
  
  #print(inflections);
  #print(d1);
  #print(d2);
  #print(count);
  #print(na);
  #print(length(w));
  #print(z);
  #print(inflections[z]);
  #print(w[inflections[z]]);
  
  
  
  
  #frequency_altitude = length(inflections) /(na);
  #stability = count/na;
  #noise_partition = 1-stability;
  #print('high frequency partition');
  #print(frequency_altitude);
  #print('stability');
  #print(stability);
  #print("noise partition");
  #print(noise_partition);
  
  #-----------------------------
  return(result2);
}