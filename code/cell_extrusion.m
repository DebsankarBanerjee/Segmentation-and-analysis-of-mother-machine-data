function [cell_id, Next] = cell_extrusion(cmy, cmy0, cell_id, ncell, Next, Ndiv, pdelta) 

%%%%%%%%%%%%%%% quantify cell extrusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   cmy
   cmy0

   eps = 1E-4;		% small number : numerical tolerence

   ncell0 = numel(cmy0);
   [mxc, dum] = size(cell_id);
   Next = 0;
   id_extrusion = 0;
   cindex = cell_id(:,1);

   del_n = ncell - ncell0 - Ndiv;

%%%%%%%%%% if no cell division has occurred %%%%%%%%%%%

   


  
   

    if(id_extrusion < eps && del_n < 0)
    
    for i=1:abs(del_n)			%----- del_n loop

    for kk=1:mxc
    if( cell_id(kk,1)==max(cindex) - (i-1) )
    id_extrusion = kk;    
    Next = Next + 1;
    fprintf('------------test  OK ,id_extrusion =%d, current index=%d \n', kk, max(cindex) - (i-1))
    break;
    end
    end

    cell_id(id_extrusion,1) = 0;
    end		%----- del_n loop

    end

%    cell_id











