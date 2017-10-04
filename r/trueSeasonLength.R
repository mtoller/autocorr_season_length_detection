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


trueSeasonLength = function(y){
  library(signal);
  library(forecast);
  library(pracma);
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
  n = nrow(y);
  na = n-1;
  #y = detrend(y);
  y = scale(y);
  a = acf.fft(y);
  
  a = a[2:length(a)];
  
  signs = a;
  
  signs[which(signs < 0)] = -1;
  signs[which(signs > 0)] = 1;
  
  delta = which(signs[2:na] + signs[1:(na-1)] == 0);
  eta = a[delta] / (-1*(a[delta+1]-a[delta]));
  zeros = delta+eta;
  z_to_z = sort(zeros[2:na]-zeros[1:(na-1)]);
  m = matrix(0,length(z_to_z),2);
  m[,1] = z_to_z;
  
  #ma = forecast::ma(a,s,centre = FALSE)
  #entropy and compression
  #kalman
  #forecasting/prediction
  #backproject with shapelet stauchen
  #implement old algorithm in R
  #anisotropic diffusion
  par(mfrow=c(2,2));
  plot(y);
  #plot(a,type = 'l');
  hist(a,c(-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),freq = FALSE)
  print(mean(a));
  print(std(a));
  print(skewness(a));
  print(excess(a));
  #hist(y,c(-3.0,-2.0,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0));
  #ma = as.ts(fastMa(a,round(s/2)));
  plot(m);
  #print(z_to_z);
  d = density(m[,1],kernel = 'epanechnikov')
  plot(d);
  #s = 2*(round(median(z_to_z)));
  #print(d);
  print(which.max(d$y));
  print(d$x);
  s = round(d$x[which.max(d$y)]*2);
  if (is.na(s))
  {
    s = round(n/10);
  }
  print(d$y*100);
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
  
  ##
  #plot(a);
  print('end');
  return(s);
}