function [miny] = min_search(imth,xr_1,xr_2,pk_dist,pk_prom)



  x_av=mean(imth(:,xr_1:xr_2),2);
  meanxav = mean(x_av);
  ch_y = x_av/mean(x_av);	
  
  %%%%%%%%%%%%%%%%%%% Local min of average column  %%%%%%%%%%%

  fprintf('New miny after segmentation cleanup: \n')

  negch_y = -ch_y + max(ch_y) + 0.1;
  [dps,dpy] = findpeaks(negch_y,'MinPeakDistance',pk_dist,'MinPeakProminence',pk_prom);	

  miny = dpy

