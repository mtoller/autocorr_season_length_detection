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

aze = function(y, harmon = FALSE, zero_discard=FALSE,skew_correction=FALSE)
{
  n = nrow(y);
  na = n-1;
  
  a = acf.fft(y);
  a = a[2:length(a)];
  a = a[1:round(length(a)*2/3)];
  
  if (harmon == TRUE)
  {
    a = a * (1/harmonic(round(na*2/3)));
  }
  
  signs = a;
  signs[which(signs < 0)] = -1;
  signs[which(signs > 0)] = 1;
  
  delta = which(signs[2:na] + signs[1:(na-1)] == 0);
  eta = a[delta] / (-1*(a[delta+1]-a[delta]));
  zeros = delta+eta;
  z_to_z = sort(zeros[2:na]-zeros[1:(na-1)]);
  
  if (zero_discard == TRUE)
  {
    z_to_z = z_to_z[which(z_to_z >= 3)];
  }
  
  d = density(z_to_z,kernel = 'epanechnikov')
  if (skew_correction && skewness(d$y) > 0)
  {
    d$y = d$y * c(1:round(length(d$y)))/sum(1:length(d$y))*length(d$y)#*skewness(d$y)*IQR(d$x) * length(d$y)/10;
  }
  zero_matrix = matrix(0,length(z_to_z),2);
  zero_matrix[,1] = z_to_z;
  par(mfrow=c(2,2));
  plot(y);
  plot(a,type = 'l');
  
  plot(zero_matrix);
  plot(d);
  
  return(d);
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
  y = detrend(y);
  y = scale(y);
  #b = butter(2,0.1);
  #y = filter(b$b,b$a,y);
  dens = aze(y);
  s = round(dens$x[which.max(dens$y)]*2);
  
  
  
  #-------------------------------------
  #ma = forecast::ma(a,s,centre = FALSE)
  #entropy and compression
  #kalman
  #forecasting/prediction
  #backproject with shapelet stauchen
  #implement old algorithm in R
  #anisotropic diffusion
  
  #hist(a,c(-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),freq = FALSE)

  #print(mean(a));
  #print(std(a));
  #print(skewness(a));
  #print(excess(a));
  
  #hist(y,c(-3.0,-2.0,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0));
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
  return(s);
}