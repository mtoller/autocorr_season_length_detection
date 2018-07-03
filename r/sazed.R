Sa <- function(y,preprocess=T)
{
  if (preprocess)
    y <- preprocessTs(y)
  return(S(computeAcf(y),F))
}
S <- function(y,preprocess=T)
{
  if (length(y) < 6)
  {
    return(1)
  }
  if (preprocess)
  {
    y <- preprocessTs(y)
  }
  n <- length(y)
  periodigram <- spec.pgram(y,detrend=F,plot=F)
  if (n >= 6)
  {
    welch <- welchPSD(as.ts(y),round(n*2/pi),windowfun = welchwindow)
  }
  else
  {
    welch <- welchPSD(as.ts(y),n,windowfun = welchwindow)
  }
  period1 <- round(1/(periodigram$freq[which.max(periodigram$spec)]))
  period2 <- round(1/(welch$frequency[which.max(welch$power)]))
  
  return(round((period1+period2)/2))
  
}
azed <- function(y,preprocess=T)
{
  if (preprocess)
    y <- preprocessTs(y)
  return(zed(computeAcf(y),F));
}
zed <- function(y,preprocess=T)
{
  if (preprocess)
  {
    y <- preprocessTs(y)
  }
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

    if (!is.list(dens) && is.nan(dens))
  {
    return(1)
  }

  return(round(dens$x[which.max(dens$y)]*2))
}

aze <- function(y,preprocess=T)
{
  if (preprocess)
    y <- preprocessTs(y)
  return(ze(computeAcf(y),F))
}

ze = function(y,preprocess=T)
{
  if (preprocess)
  {
    y <- preprocessTs(y)
  }
  
  n <- length(y)
  signs <- y
  signs[which(signs < 0)] <- -1
  signs[which(signs > 0)] <- 1
  
  zero_distance_raw <- which(signs[2:n] + signs[1:(n-1)] < 2 & signs[2:n] + signs[1:(n-1)] > -2)
  interpolation <- y[zero_distance_raw] / (-1*(y[zero_distance_raw+1]-y[zero_distance_raw]))
  zeros <- zero_distance_raw+interpolation;
  
  result = round(mean(diff(zeros)))*2
  if (is.numeric(result) & is.finite(result))
  {
    return(result)
  }
  return(0)
  
}
computeAcf <- function(y)
{
  autocorrelation <- as.ts(acf.fft(y))
  autocorrelation <- autocorrelation[2:length(autocorrelation)]
  factor <- 2/3
  return(autocorrelation[1:round(length(autocorrelation)*factor)])
}
preprocessTs <- function(y)
{
  return(scale(as.ts(detrend(as.vector(y)))))
}


sazed <- function(y,iter=0,method="alt")
{
  require(signal)
  require(forecast)
  require(pracma)
  require(bspec)
  if (anyNA(y))
  {
    print('NA is not supported')
    return(1)
  }
  if (any(is.infinite(y)) || any(is.complex(y)))
  {
    print('complex numbers are not supported')
    return(1)
  }
  if (var(y) == 0)
  {
    print('y has no variance')
    return(1)
  }

  results <- c()
  results <- c(results,S(y))
  results <- c(results,Sa(y))
  results <- c(results,ze(y))
  results <- c(results,aze(y))
  results <- c(results,zed(y))
  results <- c(results,azed(y))
  
  results <- results[which(!is.infinite(results))];
  results <- results[which(!is.na(results))];
  
  if (var(results) == 0)
  {
    return(results[1])
  }
  unique_results <- unique(results)
  unique_results <- unique_results[which(unique_results != 2)]
  unique_results <- unique_results[which(!is.infinite(unique_results))]
  tab <- tabulate(match(results,unique_results))
  if (max(tab) == 1)
  {
    if (method == "down")
    {
      return(2*sazed(downsample(y,2),method = "down"))
    }
      
    else if(method == "diff")
    {
      return(sazed(diff(y,lag = 1),method = "diff"));
    }
    
    iter <- iter + 1

    if (mod(iter,2) == 1)
    {
      return(2*sazed(downsample(y,2),iter))
    }
    return(sazed(diff(y,lag = 1),iter,method="alt"))
  }
  else if (length(tab[which(tab == max(tab))]) > 1)
  {
    sorted = sort(unique_results,decreasing = TRUE)
    return(sorted[which.max(tabulate(match(sort(results,decreasing = TRUE),sorted)))]);
  }

  vote = unique_results[which.max(tabulate(match(results,unique_results)))];
  return(vote);
}

downsample = function(data, window_size=2)
{
  require(zoo);
  return(ts(as.ts(rollapply(zoo(data),width=2,by=2,FUN=mean)),frequency=1))
}