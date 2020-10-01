function [imbw] = channel_segment(miny,imth,imbw,mincellsize,er,m1)

  %%%%%%%%%%%%%%%%% Segmentation + binary image %%%%%%%%%%%
  %%% segment the cells at the min-positions

%%%%%%%%%% disconnect cells from miny positions %%%%%%%%%%%%%%

  for kk=1:numel(miny)
    imth(miny(kk),:)=0;				%imth(m1+miny(kk),:)=0;		% m1+miny(kk) = kk'th minimum position
  end

%%%%%%%%%% disconnect cells via cell boundary (from edge detection)  %%%%%%%%%%

%  imedge = edge(imth,'Sobel');
%  imth(imedge>er)=0;


  imbw(imth<=er)=0;
  imbw(imth>er)=1;



  %---------------------- remove very small components

  imbw = bwareaopen(imbw,mincellsize);


