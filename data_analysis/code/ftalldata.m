  clc
  clear all

  er=1E-4;
  maxgr = 0.5;

  fid10=fopen(['../../result/Ptau.txt'],'w');
  fid11=fopen(['../../result/PdL.txt'],'w');
  fid12=fopen(['../../result/PLb.txt'],'w');
  fid13=fopen(['../../result/PLd.txt'],'w');
  %fid11=fopen(['../../result/MoD_corr.txt'],'w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     all tau data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % data : tdiv(ndiv), tau(ndiv), dx*delL(ndiv), dx*Lbir(ndiv), dx*Ldiv(ndiv), ind_div(ndiv), C_div(ndiv)*dx
  datafile=['../../result/alltau.txt'];

  c=load(datafile);

  tau = c(:,2);
  dl = c(:,3);
  lb = c(:,4);
  ld = c(:,5);

  ndiv = numel(tau);
  


%%%%%%%%%%%%%%%%%% distribution %%%%%%%%%%%%%%%%%%%%%


% ----------- histogram : tau ---------------
  
  X=tau;

  nbins = 20;
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  Pt = N;
  nt = (edges(1:end-1) + edges(2:end))/2; 

  for kk=1:numel(nt)
  fprintf(fid10,'%f %f\n', nt(kk), Pt(kk) );
  end

% ----------- histogram : del_L ---------------
  
  X=ld - lb;

  dl0=X;

  nbins = 20;
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  PdL = N;
  ndL = (edges(1:end-1) + edges(2:end))/2; 

  for kk=1:numel(ndL)
  fprintf(fid11,'%f %f\n', ndL(kk), PdL(kk) );
  end

% ----------- histogram : lb + ld ---------------
  
  X=lb;

  nbins = 20;
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  Plb = N;
  nlb = (edges(1:end-1) + edges(2:end))/2;

  for kk=1:numel(nlb)
  fprintf(fid12,'%f %f\n', nlb(kk), Plb(kk) );
  end

  X=ld;

  nbins = 20;
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  Pld = N;
  nld = (edges(1:end-1) + edges(2:end))/2;

  for kk=1:numel(nld)
  fprintf(fid13,'%f %f\n', nld(kk), Pld(kk) );
  end

% -------------- plot -----------------------

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])

  subplot(1,2,1)
  plot(nt,Pt,'ko-','MarkerFaceColor','black','MarkerSize',12)
  set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  xlabel('Cell cycle duration','FontSize',20)
  ylabel('Frequency','FontSize',20)

  subplot(1,2,2)
  bar(nt,Pt,'r')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  xlabel('Cell cycle duration','FontSize',20)
  ylabel('Frequency','FontSize',20)

  title(['N = ' num2str(ndiv)])

  print(h,'-dpng','-r250',['../plot/Tau_distribution.png'])

% ---------------------------------------------------------------
  
  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])

  plot(nlb,Plb,'ko-','MarkerFaceColor','black','MarkerSize',12)
  hold on
  plot(nld,Pld,'ro-','MarkerFaceColor','red','MarkerSize',12)
  
  set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  xlabel('L_b (black), L_d (red)','FontSize',20)
  ylabel('Frequency','FontSize',20)

  title(['N = ' num2str(ndiv)])

  print(h,'-dpng','-r250',['../plot/Lb_Ld_distribution.png'])





%%%%%%%%%%%%%%%%%% growth dynamics %%%%%%%%%%%%%%%%%%%%%

  %---------- Timer model ----------------

  X = lb / mean(lb);
  Y = tau / mean(tau);

  %---------------------------------------

  P = polyfit(X,Y,1);
  lfit = P(1)*X + P(2);
  slope = P(1);

  %---------------------------------------

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])

  plot(X,Y,'ko','MarkerFaceColor','black','MarkerSize',8)
  hold on
  plot(X,lfit,'red-','LineWidth',2)
  %set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  xlabel('L_b/<L_b>','FontSize',20)
  ylabel('\tau/<\tau>','FontSize',20)

  title(['N = ' num2str(ndiv) ' slope = ' num2str(slope)])

  print(h,'-dpng','-r250',['../plot/Timer_model_tau_Lb.png'])






  %---------- Sizer model ----------------

  X = lb / mean(lb);
  Y = ld / mean(ld);

  %---------------------------------------

  P = polyfit(X,Y,1);
  lfit = P(1)*X + P(2);
  slope = P(1);

  %---------------------------------------

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])

  plot(X,Y,'ko','MarkerFaceColor','black','MarkerSize',8)
  hold on
  plot(X,lfit,'red-','LineWidth',2)
  %set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  xlabel('L_b/<L_b>','FontSize',20)
  ylabel('L_d/<L_d>','FontSize',20)

  title(['N = ' num2str(ndiv) ' slope = ' num2str(slope)])

  print(h,'-dpng','-r250',['../plot/Sizer_model_Ld_Lb.png'])




  %---------- Adder model ----------------

  X = lb ;
  Y = dl ;

  %---------------------------------------

  P = polyfit(X,Y,1);
  lfit = P(1)*X + P(2);
  slope = P(1);

  %---------------------------------------

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])

  plot(X,Y,'ko','MarkerFaceColor','black','MarkerSize',8)
  hold on
  plot(X,lfit,'red-','LineWidth',2)
  %set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  xlabel('L_b/<L_b>','FontSize',20)
  ylabel('\Delta_L/<\Delta_L>','FontSize',20)

  title(['N = ' num2str(ndiv) ' slope = ' num2str(slope)])

  print(h,'-dpng','-r250',['../plot/Adder_model_del_L_Lb.png'])








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     growth rate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % data : Lbir, Ldiv, tau, gr
  datafile2=['../../result/allgrowth' '.txt'];
  
  c2=load(datafile2);

  lb = c2(:,1);
  gr = c2(:,4);


% ----------- histogram : growth rate ---------------
  
  X=gr;

  nbins = 20;
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  Pgr = N;
  ngr = (edges(1:end-1) + edges(2:end))/2; 


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])

  bar(ngr,Pgr,'r')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  axis([0.0 0.04 0 0.5])
  xlabel('growth rate','FontSize',20)
  ylabel('Frequency','FontSize',20)

  title(['N = ' num2str(ndiv)])

  print(h,'-dpng','-r250',['../plot/GR_distribution.png'])



  %---------- growth rate vs Lbirth ----------------

  X = lb ;
  Y = gr ;

  %---------------------------------------

  P = polyfit(X,Y,1);
  lfit = P(1)*X + P(2);
  slope = P(1);

  %---------------------------------------

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])

  plot(X,Y,'ko','MarkerFaceColor','black','MarkerSize',8)
  hold on
  plot(X,lfit,'red-','LineWidth',2)
  %set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  axis([0 14 0.0 0.15])
  xlabel('L_b','FontSize',20)
  ylabel('growth rate','FontSize',20)

  title(['N = ' num2str(ndiv) ' slope = ' num2str(slope)])

  print(h,'-dpng','-r250',['../plot/lb_GR.png'])



























