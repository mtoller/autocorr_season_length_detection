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
## @deftypefn {Function File} {@var{retval} =} findThreshold (@var{input1}, @var{input2})
## 
##  choses the best range of delta values
##  returns begin and end of that interval
##
## @seealso{}
## @end deftypefn

## Author: Maximilian Toller
## Created: 2016-05-06

function [return_value] = findThreshold (delta, limit, force_return = 0)

  if rows(delta) == 0
    return_value = 0
    return
  endif
  quotients = delta(2:end,1) ./ delta(1:end-1,1);
  
  Gamma = find(abs(quotients(2:end,1) - quotients(1:end-1,1)) > limit);
  if and(rows(find(Gamma == 1)) == 0, rows(quotients) > 0)
    Gamma = [1;Gamma];
  endif
  
  Gamma(end+1,1) = rows(quotients);
  
  
  if and(rows(Gamma) <= 2, rows(delta) <= 2)
    return_value(1,1) = delta(1,1);
    return_value(2,1) = delta(end,1);
    return;
  endif
  
  [m,index] = max(Gamma(2:end,1)-Gamma(1:end-1,1));
  
  return_value(1,1) = delta(Gamma(index)+1,1);
  return_value(2,1) = delta(Gamma(index+1)+1,1);
endfunction
