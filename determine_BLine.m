function [imth,BLine,calBLineFlag,l_sum] = determine_BLine(imth,Lend,Lstart,BLine,pdelta,bottom_line,calBLineFlag,l_sum,max_dbline,endCellDiv,n1,k)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calBLineFlag = 0, when no cells are at/on bottom_line, it becomes 1 when 
  % a cell is near bottom_line and remains 1 until that cell has been extracted
  % by flow. This dynamic calculation of bottom line enables correct tracking of
  % the last cell.

  

  l = Lend - Lstart;

  er = 1E-4 ;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ( calBLineFlag < er ) 		%----------- main if

  %------- clear from new BLine ----
  kstart = n1 - BLine;

  for kk = kstart:n1
  imth(kk,:) = 0;
  end

  fprintf('Recalibrate Bottom Line at = %d , calBLineFlag = %d \n',kstart, calBLineFlag)

  calBLineFlag = 1;
  
  if(endCellDiv < er)
  l_sum = sum( l(1:end-1) );
  else
  l_sum = sum( l(1:end-2) );
  end
  
  else					%----------- main if-else

  %------- Find new BLine -----------

  newBLfound = 0;

%%%%%%%%%%%%%%%%%%%% method 1 %%%%%%%%%%%%%%%%%%%%%%

  %for kk=1:numel(Lstart)	%%%--- for-if

  %dum =  (n1 - BLine) - Lstart(kk)
  
  %-------- nearest l-start --------

  %if ( abs(dum) < 2*pdelta )
  %BLine = n1 - Lstart(kk)
  %newBLfound = 1;
  %fprintf('Nearest l-start selected as new BLine \n')
  %break;
  %end

  %end				%%%--- for-if

%%%%%%%%%%%%%%%%%%%% method 2 %%%%%%%%%%%%%%%%%%%%%%

  if ( newBLfound < er)		%%%--- if

  l_sum_k = [];
  l_sum

  for kk=1:numel(Lstart)	%%%--- for
  l_sum_k(kk) = sum(Lend(1:kk) - Lstart(1:kk)) 
  end				%%%--- for-end

  %-------- find by l_sum ----------

  dl_sum = abs(l_sum - l_sum_k);
  dl_sum
  [mindl_sum,idx] = min(dl_sum) 

  if ( idx < numel(Lstart) )
  dl_BLine = abs( BLine - (n1 - Lstart(idx+1)) ) 
  end

  if (idx < numel(Lstart) && dl_BLine < max_dbline)
  BLine = n1 - Lstart(idx+1)
  new_l_sum = l_sum_k(idx);
  l_sum = new_l_sum;
  newBLfound = 1;
  fprintf('Method 2 : New BLine found by l_sum \n')
  end

  end				%%%--- if-end

%%%%%%%%%%%%%%%%%%%% method 3 %%%%%%%%%%%%%%%%%%%%%%

  %-------- new BLine not found, assume the last l-start ----
  if ( newBLfound < er && endCellDiv < er)
  BLine = n1 - Lstart(end)
  fprintf('Method 3 : New BLine assumed -> last l-start \n')
  elseif ( newBLfound < er && endCellDiv > er )
  BLine = n1 - Lend(end)
  fprintf('Method 3 : New BLine assumed -> last l-end (last cell divided) \n')
  end



  %------- clear from new BLine ----
  kstart = n1 - BLine;

  for kk = kstart:n1
  imth(kk,:) = 0;
  end

  fprintf('Recalibrate Bottom Line at = %d , BLine = %d, calBLineFlag = %d \n',kstart, BLine, calBLineFlag)


  end					%----------- main if-else-end














  



