callSeasonLength = function(y)
{
  #return(1)
  if (is.character(y))
  {
    y = read.table(y);
    if (length(y) == 2)
    {
      y = y[2];
    }
    else
    {
      y = y[1];
    }
    y = as.ts(y);
  }
  write(y,file = 'temp',ncolumns = 1);
  system('cd ../matlab;octave --silent --no-gui --eval "callSeasonLength()"');
  system('cd ../r');
  result = as.numeric(read.table('temp')[1])
  if (result == 0)
  {
    result = 1;
  }
  return(result);
}