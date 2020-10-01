function [miny] = adaptive_min_search(imth,CC,miny,miny0,negch_y,agf,divf,pdelta,nlength,Lcsum0,Lcsum,pk_prom,pk_dist,xr_1,xr_2,meanxav,k,n1);

  er = 1e-4;
  nmergeflag = 1;

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


  n0 = numel(Lcsum0)
  n = numel(Lcsum)


  Lcsum0
  Lcsum


  k1=0;
  k0=0;

  nn=min(n0,n);

  Lsum_add=0.0;

  for kk=1:nn			%%------------------------- main for

  [kk, nn]

  if(kk+k1+1 > nn || kk+k0 > nn) break; end

%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (Lcsum(kk+k1) - Lcsum0(kk+k0) > 0 && Lcsum(kk+k1) - Lcsum0(kk+k0) < 2*pdelta)
  Lsum_add = (Lcsum(kk+k1) - Lcsum0(kk+k0))
  elseif (Lcsum(kk+k1) - Lcsum0(kk+k0) > 2*pdelta)
  Lsum_add = Lsum_add + 1
  end
  
  
  if( Lcsum(kk+k1) <= Lcsum0(kk+k0) - 2*pdelta )		%%---- main if

  if( abs(Lcsum(kk+k1+1) - Lcsum0(kk+k0)) < pdelta + Lsum_add )		%%---- merging if
  fprintf('Cell division: new cell index = %d and %d \n', kk+k1, kk+k1+1)
  k1=k1+1;
  elseif( Lcsum(kk+k1+1) > Lcsum0(kk+k0) + 2*pdelta + Lsum_add )
  fprintf('Cell div and cell merging: cell index= %d \n', kk+k1+1)
  fprintf('Lcsum (n+1) = %d and Lcsum0(n) = %d \n', Lcsum(kk+k1+1), Lcsum0(kk+k0))

  %%%%%%%%%%%%%%%%%%%% search for minima %%%%%%%%%%%%%%%%%%%%
  % searches for minima with varying peak-prominence to compensate 
  % intensity fluctuation
  
  if(kk+k1==1)
  dum1=[];
  dum2=miny;
  else  
  dum1=miny(1:kk+k1);
  dum2=miny(kk+k1+1:end);
  end

  oldmin = miny0(kk+k0)

  j1 = oldmin - 3*pdelta 
  j2 = oldmin + 3*pdelta 

  new_chy = negch_y(j1:j2)

  k1=k1+1;


  %j1 = min1 + pdelta ;
  %j2 = min2 - pdelta ;

  %new_chy = negch_y(j1:j2)
  
  for jj=1:npk							%%---- peak find loop
  pkp = pk_prom - jj*dpk 
  jj	
  [peaksize,peakpos] = findpeaks(new_chy,'MinPeakProminence',pkp)		%,'MinPeakDistance',pk_dist
  if(~isempty(peakpos))
  fprintf('New peak \n')
  new_miny = j1 + peakpos
  break
  end
  end								%%---- peak find loop

  %% single column search - when central projection fails
  if (isempty(peakpos))
  [newmin_sc] = single_column_search(imth,oldmin,pdelta,pk_prom,npk,dpk,xr_1,xr_2,meanxav,k,n1);
  new_miny = newmin_sc;
  end
  %% ----------------------------------------------------

  miny=[];
  fprintf('New min created \n')
  miny = [dum1; new_miny; dum2]

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end								%%---- merging if end

  end								%%---- main if end




%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




  [Lcsum(kk+k1), Lcsum0(kk+k0) + 2*pdelta + Lsum_add]

  if( Lcsum(kk+k1) >= Lcsum0(kk+k0) + 2*pdelta + Lsum_add )		%%---- main if
  
  %%%%---------------- no of merged cells -----------------%%

  dlc_array = []; 
  nmerge=0;
  km=1;
  for jj=kk+k0:n0
  jj
  dlc = Lcsum0(jj)-Lcsum(kk+k1)
  dlc_array = [dlc_array; dlc];
  if(abs(dlc) < Lsum_add + 2*pdelta)
  nmerge=km
  break;
  end
  km=km+1;
  end

  if (nmerge<er)		%% when nmerge search failed
  cellind=kk+k1;
  cellind0=kk+k0;
  nmergeflag = 0;
  break;
  end

  k0 = k0+(nmerge-1)

  fprintf('merging without div at cell = %d \n', kk+k1)
  fprintf('No of cell merged : %d \n', nmerge)


  n_newmin = nmerge -1;

  new_miny = zeros(1,n_newmin);


  %%%%%%%%%%%%%%%%%%%% search for minima %%%%%%%%%%%%%%%%%%%%

  if(kk+k1==1)
  min1=1;
  dum1=[];
  else  
  min1=miny(kk+k1-1);
  dum1=miny(1:kk+k1-1);
  end

  if(min1==miny(end))
  min2=n1;
  dum2=[];
  else
  min2=miny(kk+k1);
  dum2=miny(kk+k1:end);
  end


  for ii=1:n_newmin						%%---- minima search loop

  %%---- search about lost minima -------

  oldmin = miny0(kk+k0-1 - (ii-1))

  j1 = oldmin - 3*pdelta 
  j2 = oldmin + 3*pdelta 

  new_chy = negch_y(j1:j2)
  
  for jj=1:npk							%%----- peak find loop
  pkp = pk_prom - jj*dpk;
  [peaksize,peakpos] = findpeaks(new_chy,'MinPeakProminence',pkp);					%,'MinPeakDistance',pk_dist
  nsize = length(peakpos)


  if(abs(nsize-1)<er)
  fprintf('New peak  at jj=%d, pkp=%f, pheight=%f\n', jj, pkp)
  peaksize
  new_miny(1,ii) = j1 + sort(peakpos)
  break
  end
  end								%%----- peak find loop

  %% single column search - when central projection fails
  if (nsize<er)
  [newmin_sc] = single_column_search(imth,oldmin,pdelta,pk_prom,npk,dpk,xr_1,xr_2,meanxav,k,n1);
  new_miny(1,ii) = newmin_sc;
  end
  %% ----------------------------------------------------

  end								%%---- minima search loop

  miny=[];
  fprintf('New min created \n')
  miny = [dum1; new_miny'; dum2]

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  

  end								%%---- main if end
  



  end							%%---- main for end






%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




%%%%%%%%%%%%%%%%% cell merging and division: case 2  %%%%%%%%%%%%%%%%%%%%

  if(nmergeflag<er)							%%---- main case:2 if

  fprintf('cell div and merging: case 2 at cell = %d \n', cellind)

  [mindlc, dlc_ind] = min(abs(dlc_array))

  l1 = Lcsum0(cellind0+dlc_ind-1)
  l2 = Lcsum(cellind)

  if(l1>l2)
  nmerge = dlc_ind;
  else
  nmerge = dlc_ind + 1;
  end

  k0 = k0+(nmerge-1)

  fprintf('merging with div case 2 at cell = %d \n', kk+k1)
  fprintf('No of cell merged : %d \n', nmerge)


  n_newmin = nmerge -1;

  new_miny = zeros(1,n_newmin);


  %%%%%%%%%%%%%%%%%%%% search for minima %%%%%%%%%%%%%%%%%%%%

  if(kk+k1==1)
  min1=1;
  dum1=[];
  else  
  min1=miny(kk+k1-1);
  dum1=miny(1:kk+k1-1);
  end

  if(min1==miny(end))
  min2=n1;
  dum2=[];
  else
  min2=miny(kk+k1);
  dum2=miny(kk+k1:end);
  end


  for ii=1:n_newmin						%%---- minima search loop

  %%---- search about lost minima -------

  oldmin = miny0(kk+k0-1 - (ii-1))

  j1 = oldmin - 3*pdelta 
  j2 = oldmin + 3*pdelta 

  new_chy = negch_y(j1:j2)
  
  for jj=1:npk							%%----- peak find loop
  pkp = pk_prom - jj*dpk;
  [peaksize,peakpos] = findpeaks(new_chy,'MinPeakProminence',pkp);					%,'MinPeakDistance',pk_dist
  nsize = length(peakpos)


  if(abs(nsize-1)<er)
  fprintf('New peak  at jj=%d, pkp=%f, pheight=%f\n', jj, pkp)
  peaksize
  new_miny(1,ii) = j1 + sort(peakpos)
  break
  end
  end								%%----- peak find loop

  %% single column search - when central projection fails
  if (nsize<er)
  [newmin_sc] = single_column_search(imth,oldmin,pdelta,pk_prom,npk,dpk,xr_1,xr_2,meanxav,k,n1);
  new_miny(1,ii) = newmin_sc;
  end
  %% ----------------------------------------------------

  end								%%---- minima search loop

  miny=[];
  fprintf('New min created \n')
  miny = [dum1; new_miny'; dum2]






  end								%%---- main case:2 if



%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




%%%%%%%%%%%%%%%%% last cell merging %%%%%%%%%%%%%%%%%%%%

  Last=nlength(n)
  if(n0>2)
  l01=Lcsum0(n0)-Lcsum0(n0-1)
  l02=Lcsum0(n0-1)-Lcsum0(n0-2)
  elseif(n0==2)
  l01=Lcsum0(2)-Lcsum0(1)
  l02=Lcsum0(1)
  elseif(n0<2)
  return
  end
 
  [Lcsum(n), Lcsum0(n0-1) + 2*pdelta + Lsum_add]

  if( Lcsum(n) >= Lcsum0(n0-1) + 2*pdelta + Lsum_add && abs(Last-(l01 + l02)) <= 2*pdelta)		%%------- main end cell merge if
  fprintf('end cell merging at cell = %d \n', n)
  

  %%%%%%%%%%%%%%%%%%%% search for minima %%%%%%%%%%%%%%%%%%%%


  min1=miny(end);
  dum1=miny;
  min2=n1;
  dum2=[];

  n_newmin=1;


  for ii=1:n_newmin						%%---- minima search loop

  %%---- search about lost minima -------

  oldmin = miny0(end - ii+1)

  j1 = oldmin - 3*pdelta 
  j2 = oldmin + 3*pdelta 

  new_chy = negch_y(j1:j2)
  
  for jj=1:npk							%%---- peak find loop
  pkp = pk_prom - jj*dpk;
  [peaksize,peakpos] = findpeaks(new_chy,'MinPeakProminence',pkp);		%,'MinPeakDistance',pk_dist
  nsize = length(peakpos)


  if(abs(nsize-1)<er)
  fprintf('New peak  at jj=%d, pkp=%f, pheight=%f\n', jj, pkp)
  peaksize
  new_miny(1,ii) = j1 + sort(peakpos)
  break
  end

  if(nsize>1+er)						%%---- multiple peak if-end
  dumpk=zeros(1,nsize);
  for jk=1:nsize
  dumpk(jk) = abs(peakpos(jk)-oldmin);
  end
  [dumv,dumi] = min(dumpk);
  fprintf('New peak-multiple at jj=%d, pkp=%f, pheight=%f\n', jj, pkp)
  new_miny(1,ii) = j1 + peakpos(dumi)
  break
  end								%%---- multiple peak if-end

  end								%%---- peak find loop

  %% single column search - when central projection fails
  if (nsize<er)
  [newmin_sc] = single_column_search(imth,oldmin,pdelta,pk_prom,npk,dpk,xr_1,xr_2,meanxav,k,n1);
  new_miny(1,ii) = newmin_sc;
  end
  %% ----------------------------------------------------

  end								%%---- minima search loop

  miny=[];
  fprintf('New min created \n')
  miny = [dum1; new_miny'; dum2]


  end								%%------- main end cell merge if-end
  
  
end












