  clc;
  clear all;
  
  % cell_one data : current_time / length-projection / length-major-axis
  
  er=1E-4;
  dx = 50/310;
  dt = 1;

  pk_dist = 20;
  pk_prom = 1;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%% get file names %%%%%%%%%%%%%%%%%%
  
  fid = fopen('dataname');
  txt = textscan(fid,'%s','delimiter','\n');
  c = txt{1} ;
  
  
  
  ndata=numel(c)
  
  cell_growth = [];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%% LOOP ON DATA (TIF FILES) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for in=1:ndata

  fid10=fopen(['../../result/allgrowth_' char(c(in)) '.txt'],'w');

  datafile=['../../data/', char(c(in)),'.txt'];				% ,'.txt'

  lt=load(datafile);

  l=lt(:,3);
  t=lt(:,1);

  [nt, dum] = size(lt);
  for kk=1:nt-1 
  dlt(kk) = ( lt(kk+1,3) - lt(kk,3) ) / ( lt(kk+1,1) - lt(kk,1) ) ;
  end

  negdlt = -1*dlt;

  [dps,dpy] = findpeaks(negdlt,'MinPeakDistance',pk_dist,'MinPeakProminence',pk_prom);	

  %--------------- segments --------------
  nseg = numel(dpy)-1;

  

  %---------------------------------------
  for kk=1:nseg			%% for-segment

  t0 = t(dpy(kk)+1);
  x0 = l(dpy(kk)+1);

  m1 = dpy(kk)+1;
  m2 = dpy(kk+1)

  tseg = t(m1:m2)-t(m1);
  lseg = log(l(m1:m2));

  %plot (t(m1:m2)-t(m1), l(m1:m2), 'r-')
  plot (tseg, lseg, 'r-')
  hold on

  %---------------------------------------

  P = polyfit(tseg,lseg,1);
  lfit = P(1)*tseg + P(2);
  gr = P(1);

  %---------------------------------------
  % export Lb, Ld, Tau, alpha(gr)
  seg_data = [x0, l(m2), (m2-m1)*dt, gr ]
  
  cell_growth = [cell_growth ; seg_data];  

  end				%% for-segment


  end				%% for-data



  [ng, dum]=size(cell_growth);
  
  for kk=1:ng

  fprintf(fid10,'%f %f %f %f\n', cell_growth(kk,1), cell_growth(kk,2), cell_growth(kk,3), cell_growth(kk,4) );

  end



























