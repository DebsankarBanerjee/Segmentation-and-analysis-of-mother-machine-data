function [miny] = division_loss(CC,miny,l,l0,agf,divf,min_lb,m1,m2)

%% ------------ prevent loss of cell div due to fluctuation in min value ----------

 miny0 = miny;

 er = 1e-4;
 ncell = CC.NumObjects ;
 endcell=numel(l)



  s = regionprops(CC,'Centroid');
  cent = cat(1, s.Centroid);
  cent_y = cent(:,2);
  cent_y=sort(cent_y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Find division loss position (if any)  %%%%%%%%%%

 if(isempty(l0))		%%% IF(1)
  
  miny=miny;
  
 else				%%% ELSE(1)

  mm = min(numel(l),numel(l0));  

  kl0 = 0;
  kl = 0;
  kdivloss = 1;
  divloss = zeros(endcell,1);
    
  for kk=1:mm			%%% FOR-IF(2)

   if(kk+kl>numel(l) || kk+kl0>numel(l0))
   break
   end

   if ( l(kk+kl) > agf*l0(kk+kl0) )
   divloss(kdivloss)=kk+kl;
   kdivloss = kdivloss + 1;
   kl0 = kl0 + 1;
   elseif (l(kk+kl) < divf*l0(kk+kl0))
   kl = kl + 1;
   end
  fprintf('divloss, kk, kl, kl0, kdivloss = %d %d %d %d\n',kk, kl, kl0, kdivloss)
  end				%%% FOR-IF(2)

  
  %---------- correct miny -----------------

  divloss(divloss<er)=[];

  if(~isempty(divloss))		%%% IF(3)

  divloss

  

  for kk=1:numel(divloss)	%%% FOR(4)

   fprintf('minima position corrected, for = %d, cell index = %d\n',kk,divloss(kk))

   %%% IF(5)--------------------------------
   if(divloss(kk)==1)		
   dum_miny = miny;
   miny = [];
   dum2 = dum_miny(1:end);
   new_miny = round( cent_y(divloss(kk)) ) ;				%- m1;
   miny = [new_miny; dum2]
   kk
   end				
   %----------------------------------------
   %%% IF(6)--------------------------------
   if(divloss(kk)>1 && divloss(kk)<endcell-er)
   dum_miny0 = miny0;
   dum_miny = miny;
   miny = [];
   dum1 = dum_miny(1:divloss(kk)-1);
   dum2 = dum_miny(divloss(kk):end);
   new_miny = round( cent_y(divloss(kk)) ) ;				%- m1;
   miny = [dum1; new_miny; dum2]
   kk
   end
   %----------------------------------------
   %%% IF(7)--------------------------------
   if(divloss(kk)==endcell)
   dum_miny = miny;
   miny = [];
   dum1 = dum_miny(1:end);				% 1:divloss(kk)-1
   new_miny = round( cent_y(divloss(kk)) ) ;				%- m1;
   miny = [dum1; new_miny]
   kk
   end
   %----------------------------------------

  end				%%% FOR(4)

  end				%%% IF(3)

  %-----------------------------------------


 end				%%% IF(1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% check for over segmentation %%%%%%%%%%%%%%%%%%%%%%%



  
































