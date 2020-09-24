% This function plots the raw image, first segmentation and final segmented image with
% cells fitted as ellipses. 

function [] = plot_channel(IM,imbw,imbw0,centroids2,end_line,BLine,n1,n2,k,NumberImages,i,c)

  minI = 1000;
  maxI = 6000;

  
  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 16])

%----------- BG substracted raw image -----------

  subplot(1,3,1)
  imagesc(IM)	
  colormap copper
  hold on;
  line([1,n2], [end_line, end_line], 'Color', 'red','LineWidth',2);
  hold off
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

%----------- First segmented image -----------  

  subplot(1,3,2)
  imagesc(imbw0)			%imagesc(flipud(im))
  colormap copper
  hold on
  plot(centroids2(:,1), centroids2(:,2),'or','MarkerFaceColor','red','MarkerSize',12)
  hold on
  line([1,n2], [n1-BLine, n1-BLine], 'Color', 'red','LineWidth',2);
  hold off
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

%----------- Final segmented image & Ellipse fit to cells -----------

  t = linspace(0,2*pi,50);

  s = regionprops(imbw,{'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation'})

  subplot(1,3,3)
  imagesc(IM)			%imagesc(flipud(im))
  colormap copper
  caxis([minI maxI])
  hold on
  plot(centroids2(:,1), centroids2(:,2),'or','MarkerFaceColor','k','MarkerSize',12)
  hold on
  
  for kk = 1:length(s)
    a = s(kk).MajorAxisLength/2;
    b = s(kk).MinorAxisLength/2;
    Xc = s(kk).Centroid(1);
    Yc = s(kk).Centroid(2);
    phi = deg2rad(-s(kk).Orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    plot(x,y,'-','Color',[.8 .4 .8],'Linewidth',4)
  end
  hold off  

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r50',['../plot/imbw_',char(c(i)),'_',num2str(k,'%03d'),'_of_',num2str(NumberImages,'%03d'),'.png'])


























