function [a] = callSeasonLength(filename='../r/temp')
  ts = load(filename);
  m = seasonLength(ts);
  save('-ascii','../r/temp','m');
  a = m;
endfunction