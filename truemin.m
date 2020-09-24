function [imbw] = truemin(im_truemin,true_miny)

  %%%%%%%%%%%%%%%%% Segmentation + binary image %%%%%%%%%%%
  %%% segment the cells at the min-positions

  [dum , ns] = size(im_truemin); 

  tm = [];

  for kk=1:ns
  ch_y = im_truemin(:,kk);
  negch_y = -ch_y + max(ch_y) + 0.1;
  [dps,dpy] = findpeaks(negch_y,'MinPeakDistance',pk_dist,'MinPeakProminence',pk_prom);	
  for jj=1:numel(dpy)
  tm = [tm; dpy(jj), xr_1 + kk];
  end
  end




