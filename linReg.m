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
## @deftypefn {Function File} {@var{retval} =} linReg (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Maximilian Toller
## Created: 2016-05-06
function [X,theta] = linReg(x,y,degree=1)
  X = [ones(rows(y),1)];

  for i=1:degree
    X(:,end+1) = x.^i;
  endfor

  if log10(rcond(X'*X)) < -10
    X = 0;
    theta = 0;
    return;
  endif
  theta = (inv(X' * X))* X' * y;
endfunction