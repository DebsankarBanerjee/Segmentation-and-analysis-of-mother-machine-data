function [endCellDiv, endExt, cellIntrusion] = end_cell_check(nlength,Lstart,Lend,cent_y,L0,l_end0,divf,pdelta,max_dl_tot,k)

  endExt = 0;
  cellIntrusion = 0;

  if (k==1) 
  endCellDiv = 0;
  return;
  end

  er = 1e-4;

  n=numel(nlength);
  n0=numel(L0);

  k1 = 0;
  dum0 = 0;
  dum1 = 0;
  
  %-------------------------------------
  for kk=1:min(n,n0)

  if (kk+k1 > n) break; end

  if ( nlength(kk+k1) < divf*L0(kk) )
  dum0 = kk;
  dum1 = kk+k1;
  k1 = k1 + 1;
  end

  end
  %-------------------------------------

  if(dum0 == n0 || dum1+1 == n)
  endCellDiv = 1;
  else
  endCellDiv = 0;
  end

  %-------------------------------------
  if (endCellDiv > er)

  if ( dum1 == n ) 
  endExt = 1; 
  else
  endExt = 0;
  end

  end

  %-------------------------------------
 
  dl = sum(nlength) - sum(L0);		%% do not use absolute value 
 
  if ( dl > max_dl_tot && cent_y(end)>l_end0 )
  cellIntrusion = 1;
  end


  endExt  
  endCellDiv
  cellIntrusion










