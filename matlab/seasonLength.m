## Copyright (C) 2016 Maximilian Toller
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{retval} =} seasonLength(
##
## @var y: time series data
##
## @var x: x values corresponding to y
##
## @var degree: degree of the linear regression
##
## @var butter1: parameter b of a butterworth filter
##
## @var butter2: parameter a of a butterworth filter
##
## @var force_return: if not zero, then the function will return a value
##  even if no season length was found
##
## )
##
## Computes the season length of time series y
##
## @seealso{}
## @end deftypefn




## Author: Maximilian Toller
## Created: 2016-05-02

function [return_value] = seasonLength (y,x=0, degree=1, butter1=2,butter2=0.001, force_return = 0)
  pkg load tsa;
  pkg load signal;
  
  if x == 0
    x = [1:rows(y)]';
  endif
  
  #choosing degree of trend
  [X1, theta1] = linReg(x,y,1);
  [X2, theta2] = linReg(x,y,2);
  mse1 = meansq((y-X1*theta1));
  mse2 = meansq((y-X2*theta2));
  if log(abs(mse1-mse2)) > e^2
    degree=2;
  endif
  
  #expand data
  if (and(and(rows(x) > 2,rows(x) < 10000), abs(x(2,1)-x(1,1)) == 1))
    [x,y] = expandData(x,y);
  endif
  
  
  #create butterworth filter
  [b,a] = butter(butter1,butter2);
  
  
  #removing white noise
  ys = filter(b,a,y);
  
  no_filter = 0;
  do 
    stop = 1;
    #calculating autocorrelation
    yt = acorf(detrend(ys, degree)',rows(ys)-1)';

    
    #plot(yt,'k');
    #hold on;
    #plot(detrend(yt),'r')
    
    #expanding data by 1 row for regression
    #yt = [yt;yt(end)];
    
    #plot(yt,'g')
    #plot(detrend(yt),'c')
    
    #uiwait;
    #hold off;
    

    
    #solving linear regression
    #[X,theta] = linReg(x,yt,1);
    
    #intersecting regression with autocorrelation
    #alpha = find(abs(yt-X*theta) < 0.001);
    alpha = find(abs(detrend(yt)) < 0.001);
    
    #calculating distances between points of intersection
    delta = sort(alpha(2:end,1) - alpha(1:end-1,1));
   
      
    #removing numeric inaccuracies
    delta = delta(find(delta > 1));
    
    
    if and(rows(delta) < 6,no_filter == 0)
      stop = 0;
      no_filter = 1;
      ys = y;
      continue;
    endif
  until stop == 1  
  #finding thresholds for relevant data sequence
  threshold = findThreshold(delta,0.05, force_return); 
    
  if threshold == 0
    return_value = 0;
    return
  endif
  
  delta = delta(find(and(delta >= threshold(1,1),delta <= threshold(2,1))));
  
  #removing potetial outliers
  return_value = 2 * mean(quantile(delta,[0.1,0.25,0.5,0.75,0.9],1,METHOD=8));
  
  if and(return_value >= 0.5*rows(y), force_return == 0)
    return_value = 0;
    return;
  endif
  
  #scaling result to x-axis
  return_value *= abs(x(2,1) - x(1,1));
  
  
  #rounding to nearest integer value
  return_value = round(return_value);
  
  #there is no season length <= 1
  if return_value <= 1
    return_value = 0;
  endif
  
endfunction
