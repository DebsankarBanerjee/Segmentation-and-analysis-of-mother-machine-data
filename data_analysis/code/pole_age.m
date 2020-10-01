clear all
clc

er = 1e-4;

  fid10=fopen(['../result/opm.txt'],'w');
  fid11=fopen(['../result/npd.txt'],'w');
  fid12=fopen(['../result/position_growth.txt'],'w');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pole age effect on bacterial growth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % data : tdiv(ndiv), tau(ndiv), dx*delL(ndiv), dx*Lbir(ndiv), dx*Ldiv(ndiv), ind_div(ndiv), C_div(ndiv)*dx
  datafile=['../result/alltau.txt'];

  c=load(datafile);

  tau = c(:,2);
  dl = c(:,3);
  lb = c(:,4);
  ld = c(:,5);

  indx = c(:,6);
  cy = c(:,7);

  % growth rate calculation----------------------------------------------

  n = numel(tau);
  gr = zeros(n,1);

  for kk=1:n
  gr(kk) = (1.0/tau(kk))*log(ld(kk)/lb(kk));
  end 

  % find 'old pole mother'(opm) and 'new pole daughter'(npd) ----------------------

  nold=0;
  nnew=0;
  data_opm = [];
  data_npd = [];

  for kk=1:n

  if (indx(kk)-1 < er)
  nold = nold + 1;
  data_opm = [data_opm; lb(kk) dl(kk) gr(kk) tau(kk)];
  end

  if (indx(kk)-2 < er)
  nnew = nnew + 1;
  data_npd = [data_npd; lb(kk) dl(kk) gr(kk) tau(kk)];
  end

  end

  % no of generation - same as the index in opm, npd
  X1 = [1:nold];
  X2 = [1:nnew];

  % plot opm vs npd : pole age effect (PAE)------------------------------

  % delta_L plot --------------------------------------------------------

  Y1 = data_opm(:,2);
  Y2 = data_npd(:,2);

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
  
  
  plot(X1,Y1,'red-', X2,Y2,'blue-','LineWidth',3) 

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/PAE_dl.png'])



  % growth rate plot --------------------------------------------------------

  Y1 = data_opm(:,3);
  Y2 = data_npd(:,3);

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
  
  
  plot(X1,Y1,'red-', X2,Y2,'blue-','LineWidth',3) 

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/PAE_gr.png'])



  % division time plot --------------------------------------------------------

  Y1 = data_opm(:,4);
  Y2 = data_npd(:,4);

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
  
  
  plot(X1,Y1,'red-', X2,Y2,'blue-','LineWidth',3) 

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/PAE_tau.png'])


  % birth length plot --------------------------------------------------------

  Y1 = data_opm(:,1);
  Y2 = data_npd(:,1);

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
  
  
  plot(X1,Y1,'red-', X2,Y2,'blue-','LineWidth',3) 

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/PAE_bl.png'])

  % distribution plots -------------------------------------------------------

  nbins = 20;

  % ----------- histogram : lb ---------------
  
  X=data_opm(:,1);
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  P1 = N;
  n1 = (edges(1:end-1) + edges(2:end))/2; 


  X=data_npd(:,1);
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  P2 = N;
  n2 = (edges(1:end-1) + edges(2:end))/2; 


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 8])

  plot(n1,P1,'ko-',n2,P2,'ro-','MarkerFaceColor','black','MarkerSize',12)
  %set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/PAE_dist_lb.png'])

  %for kk=1:numel(nt)
  %fprintf(fid10,'%f %f\n', nt(kk), Pt(kk) );
  %end


  % ----------- histogram : dl ---------------
  
  X=data_opm(:,2);
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  P1 = N;
  n1 = (edges(1:end-1) + edges(2:end))/2; 


  X=data_npd(:,2);
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  P2 = N;
  n2 = (edges(1:end-1) + edges(2:end))/2; 


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 8])

  plot(n1,P1,'ko-',n2,P2,'ro-','MarkerFaceColor','black','MarkerSize',12)
  %set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/PAE_dist_dl.png'])



  % write data ----------------------------------------------------------------

  for kk=1:nold
  fprintf(fid10,'%f %f %f %f\n', data_opm(:,1), data_opm(:,2), data_opm(:,3), data_opm(:,4) );
  end

  for kk=1:nnew
  fprintf(fid11,'%f %f %f %f\n', data_npd(:,1), data_npd(:,2), data_npd(:,3), data_npd(:,4) );
  end



  %%%%%%%%%%%%%%%%%%%% effects of division position %%%%%%%%%%%%%%%%%

  X = cy;

  % delta_L plot --------------------------------------------------------

  Y = dl;

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/position_dl.png'])


  % growth rate plot --------------------------------------------------------

  Y = gr;

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/position_gr.png'])


  % division time plot --------------------------------------------------------

  Y = tau;

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/position_tau.png'])


  % birth length plot --------------------------------------------------------

  Y = lb;

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/position_lb.png'])



  % write data ---------------------------------------------------------------


  for kk=1:n
  fprintf(fid12,'%f %f %f %f %f\n', cy(kk), lb(kk), dl(kk), gr(kk), tau(kk) );
  end













