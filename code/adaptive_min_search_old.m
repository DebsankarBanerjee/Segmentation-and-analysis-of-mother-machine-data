function [miny] = adaptive_min_search(CC,miny,miny0,negch_y,agf,divf,pdelta,nlength,Lcsum0,Lcsum,pk_prom,pk_dist,k,n1);

  er = 1e-4;
  n_merge = 0;
  k_merge = 0;
  cellMerged = 0;

  npk = 10;
  dpk = (0.5*pk_prom)/10.0 ;

  %--------------- first frame -------------
  if (k==1) 
  miny = miny;
  return;
  end
  %-----------------------------------------


  s  = regionprops(CC,'Centroid');
  centroids = cat(1, s.Centroid);
  cent_y=centroids(:,2);
  cent_y=sort(cent_y);


  n0 = numel(Lcsum0);
  n = numel(Lcsum);

%  for kk=1:n0

%  p1 = miny0(kk)-2*pdelta;
%  p2 = miny0(kk)+2*pdelta;

  







  Lcsum0
  Lcsum


  k1=0;
  k0=0;

  nn=min(n0,n);

  for kk=1:nn						%%---- main for

  if(kk+k1 > nn || kk+k0 > nn) break; end

%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  if( Lcsum(kk+k1) <= Lcsum0(kk+k0) - 2*pdelta )		%%---- main if

  if( abs(Lcsum(kk+k1+1) - Lcsum0(kk+k0)) < pdelta )		%%---- merging if
  fprintf('Cell division: new cell index = %d and %d \n', kk+k1, kk+k1+1)
  k1=k1+1;
  elseif( Lcsum(kk+k1+1) > Lcsum0(kk+k0) + 2*pdelta )
  fprintf('Cell div and cell merging: cell index= %d \n', kk+k1+1)

  %%%%%%%%%%%%%%%%%%%% search for minima %%%%%%%%%%%%%%%%%%%%
  
  min1=miny(kk+k1);
  dum1=miny(1:kk+k1);

  if(min1==miny(end))
  min2=n1;
  dum2=[];
  else
  min2=miny(kk+k1+1);
  dum2=miny(kk+k1+1:end);
  end

  j1 = min1 + pdelta ;
  j2 = min2 - pdelta ;

  new_chy = negch_y(j1:j2)
  
  for jj=1:npk
  pkp = pk_prom - jj*dpk 
  jj	
  [peaksize,peakpos] = findpeaks(new_chy,'MinPeakProminence',pkp)		%,'MinPeakDistance',pk_dist
  if(~isempty(peakpos))
  fprintf('New peak \n')
  new_miny = j1 + peakpos
  break
  end
  end

  miny=[];
  fprintf('New min created \n')
  miny = [dum1; new_miny; dum2]

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end								%%---- merging if end

  end								%%---- main if end

%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  if( Lcsum(kk+k1) >= Lcsum0(kk+k0) + 2*pdelta )		%%---- main if

  fprintf('merging without div at cell = %d \n', kk+k1)
  
  %%%%---------------- no of merged cells -----------------%%

  km=1;
  for jj=kk+k0:nn
  dlc = Lcsum0(jj)-Lcsum(kk+k1);
  if(abs(dlc) < 2*pdelta)
  nmerge=km;
  break;
  end
  km=km+1;
  end

  k0 = k0+(nmerge-1)
  fprintf('No of cell merged : %d \n', nmerge)


  n_newmin = nmerge -1;

  new_miny = zeros(1,n_newmin);


  %%%%%%%%%%%%%%%%%%%% search for minima %%%%%%%%%%%%%%%%%%%%
  
  min1=miny(kk+k1-1);
  dum1=miny(1:kk+k1-1);

  if(min1==miny(end))
  min2=n1;
  dum2=[];
  else
  min2=miny(kk+k1);
  dum2=miny(kk+k1:end);
  end

  j1 = min1 + pdelta 
  j2 = min2 - pdelta 

  new_chy = negch_y(j1:j2)
  
  for jj=1:npk
  pkp = pk_prom - jj*dpk;
  [peaksize,peakpos] = findpeaks(new_chy,'MinPeakProminence',pkp);		%,'MinPeakDistance',pk_dist
  nsize = length(peakpos);


  if(abs(nsize-n_newmin)<er)
  fprintf('New peak  at jj=%d, pkp=%f, pheight=%f\n', jj, pkp)
  peaksize
  new_miny = j1 + sort(peakpos)
  break
  end
  end

  miny=[];
  fprintf('New min created \n')
  miny = [dum1; new_miny; dum2]

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end								%%---- merging if end

  end								%%---- main if end
  



  end							%%---- main for end
























