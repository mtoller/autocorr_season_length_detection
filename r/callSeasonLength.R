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
  write(y,file = paste0('..',pathSymbol,'matlab',pathSymbol,'temp'),ncolumns = 1);
  system(paste0('cd ..',pathSymbol,'matlab;octave --silent --no-gui --eval "callSeasonLength()"'));
  system(paste0('cd ..',pathSymbol,'r'));
  result = as.numeric(read.table(paste0('..',pathSymbol,'matlab',pathSymbol,'temp'))[1])
  if (result == 0)
  {
    result = 1;
  }
  return(result);
}