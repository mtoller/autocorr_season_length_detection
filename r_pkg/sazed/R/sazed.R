#' sazed: A package for for estimating the season length of a seasonal
#' time series.
#'
#' The sazed package provides the main function to compute season length,
#' sazed, which is an ensemble of many season length estimation methods,
#' also included in this package.
#'
#' @docType package
#' @name sazed
#' @import stats
NULL

#' Compute the SA component of the SAZED ensemble
#'
#' \code{Sa} computes the autocorrelation of its argument, and then derives the
#' season length from its spectral density.
#'
#' @param y The input time series.
#' @param preprocess If true, y is detrended and z-normalized before
#'   computation.
#' @return The SA season length estimate of y.
#' #@examples
#' season_length <- 26
#' y <- sin(1:400*2*pi/season_length)
#' Sa(y)
#' Sa(y, preprocess = FALSE)
#' @export
Sa <- function(y,preprocess=T)
{
  if (preprocess)
    y <- preprocessTs(y)
  return(S(computeAcf(y),F))
}
#' Compute the S component of the SAZED ensemble
#'
#' \code{S} computes the spectral density of its argument, and then derives the
#' season length from it.
#'
#' @param y The input time series.
#' @param preprocess If true, y is detrended and z-normalized before
#'   computation.
#' @return The S season length estimate of y.
#' #@examples
#' season_length <- 26
#' y <- sin(1:400*2*pi/season_length)
#' S(y)
#' S(y, preprocess = FALSE)
#' @export
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
    welch <- bspec::welchPSD(as.ts(y),round(n*2/pi),windowfun = bspec::welchwindow)
  }
  else
  {
    welch <- bspec::welchPSD(as.ts(y),n,windowfun = bspec::welchwindow)
  }
  period1 <- round(1/(periodigram$freq[which.max(periodigram$spec)]))
  period2 <- round(1/(welch$frequency[which.max(welch$power)]))

  return(round((period1+period2)/2))

}
#' Compute the AZED component of the SAZED ensemble
#'
#' \code{azed} computes the autocorrelation of its argument, and then derives the
#' season length from its the autocorrelations zero density.
#'
#' @param y The input time series.
#' @param preprocess If true, y is detrended and z-normalized before
#'   computation.
#' @return The AZED season length estimate of y.
#' #@examples
#' season_length <- 26
#' y <- sin(1:400*2*pi/season_length)
#' azed(y)
#' azed(y, preprocess = FALSE)
#' @export
azed <- function(y,preprocess=T)
{
  if (preprocess)
    y <- preprocessTs(y)
  return(zed(computeAcf(y),F));
}
#' Compute the ZED component of the SAZED ensemble
#'
#' \code{zed} computes the zero density of its argument, and then derives the
#' season length from it.
#'
#' @param y The input time series.
#' @param preprocess If true, y is detrended and z-normalized before
#'   computation.
#' @return The ZED season length estimate of y.
#' #@examples
#' season_length <- 26
#' y <- sin(1:400*2*pi/season_length)
#' zed(y)
#' zed(y, preprocess = FALSE)
#' @export
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
#' Compute the AZE component of the SAZED ensemble
#'
#' \code{aze} estimates the season length of its argument from the mean autocorrelation zero
#' distance
#'
#' @param y The input time series.
#' @param preprocess If true, y is detrended and z-normalized before
#'   computation.
#' @return The AZE season length estimate of y.
#' #@examples
#' season_length <- 26
#' y <- sin(1:400*2*pi/season_length)
#' aze(y)
#' aze(y, preprocess = FALSE)
#' @export
aze <- function(y,preprocess=T)
{
  if (preprocess)
    y <- preprocessTs(y)
  return(ze(computeAcf(y),F))
}
#' Compute the ZE component of the SAZED ensemble
#'
#' \code{ze} estimates the season length of its argument from the mean zero distance
#'
#' @param y The input time series.
#' @param preprocess If true, y is detrended and z-normalized before
#'   computation.
#' @return The ZE season length estimate of y.
#' #@examples
#' season_length <- 26
#' y <- sin(1:400*2*pi/season_length)
#' ze(y)
#' ze(y, preprocess = FALSE)
#' @export
ze <- function(y,preprocess=T)
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
#' Compute and shorten autocorrelation
#'
#' \code{computeAcf} computes the autocorrelation function of its argument and discards
#' the zero lag and all lags greater than 2/3 of the argument's length
#'
#' @param y The input time series.
#'
#' @return The shortened autocorrelation
#' #@examples
#' season_length <- 26
#' y <- sin(1:400*2*pi/season_length)
#' computeAcf(y)
#' #@export
computeAcf <- function(y)
{
  autocorrelation <- as.ts(acf.fft(y))
  autocorrelation <- autocorrelation[2:length(autocorrelation)]
  factor <- 2/3
  return(autocorrelation[1:round(length(autocorrelation)*factor)])
}
#' Preprocess Time Series for SAZED ensemble
#'
#' \code{preprocessTs} detrends and z-normalizes its argument.
#'
#' @param y The input time series.
#'
#' @return The detrended and z-normalized time series.
#' #@examples
#' season_length <- 26
#' y <- sin(1:400*2*pi/season_length)
#' preprocessTs(y)
#' #@export
preprocessTs <- function(y)
{
  return(scale(as.ts(pracma::detrend(as.vector(y)))))
}
#' SAZED Ensemble
#'
#' \code{sazed} estimates a time series' season length by computing 6 different estimates
#' and taking a majority vote.
#'
#' @param y The input time series.
#' @param iter The recursion depth.
#' @param method The method used for breaking ties. One of \code{c("alt","diff","down")}.
#'
#' @return The season length of the input time series.
#' @examples
#' season_length <- 26
#' y <- sin(1:400*2*pi/season_length)
#' sazed(y)
#' @export
sazed <- function(y,iter=0,method="alt")
{
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

    if (pracma::mod(iter,2) == 1)
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
#' Downsample Time Series
#'
#' \code{downsample} samples down a time series with a rolling mean.
#'
#' @param data The input time series.
#' @param window_size The size of the rolling mean window used.
#' @return The downsampled time series.
#' @importFrom zoo zoo rollapply
downsample <- function(data, window_size=2)
{
  return(ts(as.ts(rollapply(zoo(data),width=2,by=2,FUN=mean)),frequency=1))
}
