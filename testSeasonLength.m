## Copyright (C) 2016 dev
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
## @deftypefn {Function File} {@var{retval} =} testSeasonLength (
##
## @var time_series: name of the file where the time series is stored
##
## @var dim: degree of the linear regression
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
## Tests the function seasonLength with a time series from the test cases
##
## @seealso{}
## @end deftypefn



## Author: dev <dev@dev-PC>
## Created: 2016-05-02



function [retval] = testSeasonLength(time_series, dim=1, butter1=2,butter2=0.001, force_return = 0)
  pkg load tsa;
  
  #preparing data
  
  ts = load(time_series);
  x = ts(:,1);
  y = ts(:,2);

  y_max = max(y);
  y_min = min(y);
  [_ i_max] = max(detrend(y));

  #getting seasonLength

  s = seasonLength(y,x,dim,butter1,butter2,force_return);

  #evaluating results
  
  
  if and(s == 0, force_return == 0)
    retval = 0
    disp('No seasonality or too little data');
    plot(x,y);
    return;
  endif
  
  #plotting results
  s = round(s);
  plot(x,y);
  hold on;

  line_pos = i_max - floor(i_max/s)*s;
  
  if line_pos == 0
    line_pos = 1;
  endif
  
  for i = x(line_pos,1):s:x(end,1)
    plot([i,i],[y_min,y_max],'r');
  endfor

  hold off;

  retval = s;

endfunction
