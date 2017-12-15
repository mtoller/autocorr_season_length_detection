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
## @deftypefn {Function File} {@var{retval} =} expandData (
##
## @var x x-values of the time series, choose [1:n]' if you have only y
##
## @var y y-values of the time series
##
## )
##
## Interpolates a time series to increase its size 100 fold
##
##
## @seealso{}
## @end deftypefn

## Author: Maximilian Toller
## Created: 2016-05-02

function [ex,ey] = expandData (x,y)
  if rows(x) != rows (y)
    disp("dimension mismatch");
    return_value = 0;
    return;
  endif
  
  ex = [x(1,1):0.01:x(end,1)]';
  ey = interp1(x,y,ex);
endfunction
