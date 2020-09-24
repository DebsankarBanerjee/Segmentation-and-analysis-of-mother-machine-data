function [cell_id, ncid, Ndiv] = cell_division(L, L0, ma_L, ma_L0, cell_id, time, cent_y0, ncell, ncid, divf,fil_length,fil_div, pdelta) 


%%%%%%%%%%%%%%% quantify cell division %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds the cell id of dividing cells, assign new cell ids, records 
% time of division, division size, birth size etc in the matrix cell_id


   addi = 0;		% shift in current index due cell-div
   Ndiv = 0;
   ncell0 = numel(L0);
   [mxc, dum] = size(cell_id);
	

   for i=1:numel(L0)		%------------------- start of i-loop : loop on old cells

   if(i+addi>numel(L))
      break
    end

%---------- locate cell division -----------

    if( L(i+addi) <= divf*L0(i) || ( L0(i)>fil_length && L0(i)-L(i+addi)>fil_div ) )	%----------- cell-div if loop 

    Ndiv = Ndiv + 1;

%---------------------------- find parent id and kill the parent cell 
    for kk=1:mxc
    if(cell_id(kk,1)==i+addi)
    id_parent = kk;
    fprintf('------------test  OK ,parent  id=%d, index=%d \n', kk,i)
    break;
    end
    end

    cell_id(id_parent,1) = 0;
    cell_id(id_parent,5) = ma_L0(i);
    cell_id(id_parent,3) = time;
    cell_id(id_parent,7) = i;
    cell_id(id_parent,8) = cent_y0(i);

%------------------------------- SHIFT CURRENT INDICES

    for kk = 1:mxc
    if(cell_id(kk,1)>i+addi)		%>i
    cell_id(kk,1) = cell_id(kk,1) + 1;
    end
    end

%%%%%%%%%%  NEW CELLS  %%%%%%%%

%---------------------------- FIRST DAUGHTER
   
    cell_id(ncid+1,1) = i+addi;
    cell_id(ncid+1,2) = time;
    cell_id(ncid+1,4) = ma_L(i+addi);
    cell_id(ncid+1,6) = id_parent;
    

%---------------------------- SECOND DAUGHTER
    
    cell_id(ncid+2,1) = i+addi+1;
    cell_id(ncid+2,2) = time;
    cell_id(ncid+2,4) = ma_L(i+addi+1);
    cell_id(ncid+2,6) = id_parent;


%    cell_id

%---------------------------- 
%---------------------------- new identity of cell
    ncid = ncid + 2;    
    addi = addi + 1
    end			     %----------- cell-div if loop

   end   %------------------- end of i-loop
   


