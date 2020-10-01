function [cell_id] = initiate_cell(cell_id, cn)

  for kk=1:cn
  cell_id(kk,1) = kk;		%----- current location
  %cell_id(kk,7) = 1;		%----- alive
  end

