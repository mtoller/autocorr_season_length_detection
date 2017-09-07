trueSeasonLength = function(y){
  library(signal);
  library(forecast);
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
  y = scale(y);
  
  a = acf.fft(y);
  
  a = a[2:length(a)];
  
  signs = a;
  
  signs[which(signs < 0)] = -1;
  signs[which(signs > 0)] = 1;
  
  delta = which(signs[2:n] + signs[1:(n-1)] == 0);
  eta = a[delta] / (-1*(a[delta+1]-a[delta]));
  zeros = delta+eta;
  z_to_z = sort(zeros[2:n]-zeros[1:(n-1)]);
  m = matrix(0,length(z_to_z),2);
  m[,1] = z_to_z;
  s = 2*(round(median(z_to_z)));
  #ma = forecast::ma(a,s,centre = FALSE)
  
  print(abs(length(ma)-length(a)));
  par(mfrow=c(3,2));
  plot(y);
  plot(a,type = 'l');
  hist(a,c(-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),freq = FALSE)
  #hist(y,c(-3.0,-2.0,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0));
  ma = as.ts(fastMa(a,round(s/2)));
  plot(ma);
  plot(m);
  plot(a-ma);
  #plot(a);
  print('end');
  return(s);
}