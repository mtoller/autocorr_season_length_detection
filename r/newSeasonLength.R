newSeasonLength = function(y)
{
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
    #print('y already has a frequency');
    #return(frequency(y));
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
  #plot(y[1:100],type='l');
  autocorrelation = acf.fft(y);
  autocorrelation = autocorrelation[2:length(autocorrelation)];
  #autocorrelation = autocorrelation[1:round(length(autocorrelation)*2/3)];
  
  
  zeros = findAZeros(autocorrelation);
  p = findAPeaks(autocorrelation);
  v = findAValleys(autocorrelation);
  #print(zeros);
  #print(p);
  #print(v);
  z = round(10*zeros);
  p = 10*p;
  v = 10*v;
  
  symbols = rep.int(-2,10*length(autocorrelation));
  
  symbols[z] = 0;
  symbols[p] = 1;
  symbols[v] = -1;
  #print(symbols[which(symbols != -2)]);
  #print(which(symbols != -2));
  
  max_len = floor(length(autocorrelation)/3); # maximum possible season length
  min_len = 2; #minimum possible season length
  
  best_offset = -1;
  best_candidate_len = -1;
  min_error = Inf;
  #cat(paste0('max_len: ', max_len,'\nlength_zeros: ',length(zeros),'\n'));
  for (offset in (0:max_len))
  {
    if (TRUE)
      break;
    for (candidate_len in (min_len:max_len))
    {
      number_of_seasons = idivide(length(autocorrelation),candidate_len,rounding = 'floor')
      if (offset >= candidate_len)
      {
        number_of_seasons = number_of_seasons - 1;
      }
      candidate_sequence = rep.int(offset,number_of_seasons);
      candidate_sequence = candidate_sequence + 1:number_of_seasons * candidate_len;
      
      if (length(candidate_sequence) > length(zeros))
      {
        shorter = zeros;
        longer = candidate_sequence;
      }
      else
      {
        shorter = candidate_sequence;
        longer = zeros;
      }
      
      candidate_matrix = t(matrix(longer,length(longer),length(shorter)));
      difference_matrix = abs(candidate_matrix - shorter);
      min_of_rows = apply(difference_matrix,1,function(x) which.min(x[]));
      closest_sequence = t(candidate_matrix)[min_of_rows];
      
      error = mean((closest_sequence-shorter)^2);
      error = error ^ (1+abs(length(shorter)-length(longer)));
      if (error < min_error)
      {
        min_error = error;
        best_offset = offset;
        best_candidate_len = candidate_len;
        #print(candidate_sequence);
        #print(closest_sequence);
        
      }
      
    }
  }
  
  #cat(paste0('min_error: ', log(min_error),'\nbest_offset: ', best_offset, "\nbest_candidate_len: ",
  #           best_candidate_len),'\n');
  #print(candidate_sequence);
  #print(closest_sequence);
  r1 = round(mean(diff(round(zeros))))*2;
  r2 = round(mean(diff(round(p/10))));
  r3 = round(mean(diff(round(v/10))));
  #print(r1);
  #print(r2);
  #print(r3);
  return(r1);
}

findAZeros = function(y)
{
  n = length(y);
  signs = y;
  signs[which(signs < 0)] = -1;
  signs[which(signs > 0)] = 1;
  
  zero_distance_raw = which(signs[2:n] + signs[1:(n-1)] == 0);
  interpolation = y[zero_distance_raw] / (-1*(y[zero_distance_raw+1]-y[zero_distance_raw]));
  zero_distance_exact = zero_distance_raw+interpolation;
  #print('zeros');
  #print(zero_distance_exact);
  
  return(zero_distance_exact);
}
findAPeaks = function(y)
{
  difference = diff(y,1);
  difference[which(difference > 0)] = 0.5;
  difference[which(difference < -0)] = -0.5;
  #difference[which(difference <= epsilon)] = 0;
  #difference[which(difference >= -epsilon)] = 0;
  difference = diff(difference,1);
  return(sort(which(difference == 1)));
}
findAValleys = function(y)
{
  difference = diff(y,1);
  difference[which(difference > 0)] = 0.5;
  difference[which(difference < -0)] = -0.5;
  #difference[which(difference <= epsilon)] = 0;
  #difference[which(difference >= -epsilon)] = 0;
  difference = diff(difference,1);
  return(sort(which(difference == -1)));
}