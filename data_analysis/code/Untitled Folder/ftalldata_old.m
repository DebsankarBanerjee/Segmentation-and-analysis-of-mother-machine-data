
  er=1E-4;


  fid10=fopen(['../result/Ptau.txt'],'w');
  fid11=fopen(['../result/PdL.txt'],'w');
  fid12=fopen(['../result/PLb.txt'],'w');
  fid13=fopen(['../result/PLd.txt'],'w');
  %fid11=fopen(['../result/MoD_corr.txt'],'w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     all tau data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % data : tdiv(ndiv), tau(ndiv), dx*delL(ndiv), dx*Lbir(ndiv), dx*Ldiv(ndiv)
  datafile=['../result/alltau.txt'];

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

  subplot(1,2,2)
  bar(nt,Pt,'r')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

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

  title(['N = ' num2str(ndiv) ' slope = ' num2str(slope)])

  print(h,'-dpng','-r250',['../plot/Adder_model_del_L_Lb.png'])








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     growth rate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % data : Lbir, Ldiv, tau, gr
  datafile2=['../result/allgrowth' '.txt'];
  
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

  title(['N = ' num2str(ndiv) ' slope = ' num2str(slope)])

  print(h,'-dpng','-r250',['../plot/lb_GR.png'])








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     correlation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MDD data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tau %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % data : tau_m, tau_d1, tau_d2

  datafile=['../result/allmdd_tau.txt'];
  c1=load(datafile);

  tm = c1(:,1);
  td1 = c1(:,2);
  td2 = c1(:,3);

  %---------------------------------

  X = [tm; tm];
  Y = [td1; td2];

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc_mdd = cxy / (xst*yst) ;

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc_mdd) ' N=' num2str(nd)])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/mdd_tau_corr.png'])


 

  %----------------------------------

  X = td1;
  Y = td2;

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc_dd = cxy / (xst*yst) ;

  nd = numel(X);
 
  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc_dd) ' N=' num2str(nd)])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/dd_tau_corr.png'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% del_L and r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % data : lb_m, lb_d1, lb_d2

  datafile=['../result/allmdd_lb.txt'];
  c1=load(datafile);

  lbm = c1(:,1);
  lbd1 = c1(:,2);
  lbd2 = c1(:,3);


  datafile=['../result/allmdd_ld.txt'];
  c2=load(datafile);

  ldm = c2(:,1);
  ldd1 = c2(:,2);
  ldd2 = c2(:,3);

  %------------- del_L --------------------

  dlm = ldm - lbm;
  dld1 = ldd1 - lbd1;
  dlm = ldd2 - lbd2;

  %------------- growth rate --------------

  for kk = 1:numel(dlm)
  rm(kk) = (1/tm(kk))*log(ldm(kk)/lbm(kk));
  rd1(kk) = (1/td1(kk))*log(ldd1(kk)/lbd1(kk));
  rd2(kk) = (1/td2(kk))*log(ldd2(kk)/lbd2(kk));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %------------ del L correlations ---------------------

  X = [dlm; dlm];
  Y = [dld1; dld2];

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc_mdd = cxy / (xst*yst) ;

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
 

  title(['PCC = ' num2str(pcc_mdd) '  N=' num2str(nd)])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/mdd_del_L_corr.png'])


 

  %----------------------------------

  X = dld1;
  Y = dld2;

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc_dd = cxy / (xst*yst) ;

  nd = numel(X);
 
  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  

  title(['PCC = ' num2str(pcc_dd) '  N=' num2str(nd)])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/dd_del_L_corr.png'])




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %------------ growth rate correlations ---------------------

  X = [rm; rm];
  Y = [rd1; rd2];

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc_mdd = cxy / (xst*yst) ;

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
 

  title(['PCC = ' num2str(pcc_mdd) '  N=' num2str(nd)])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/mdd_gr_corr.png'])


 

  %----------------------------------

  X = rd1;
  Y = rd2;

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc_dd = cxy / (xst*yst) ;

  nd = numel(X);
 
  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  

  title(['PCC = ' num2str(pcc_dd) '  N=' num2str(nd)])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/dd_gr_corr.png'])



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MDC data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tau %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % data : % data : tau_m, tau_d1, tau_d2, tau_c11, tau_c12, tau_c21, tau_c22

  datafile=['../result/allmdc_tau.txt'];
  c2=load(datafile);

  c11 = c2(:,4);
  c12 = c2(:,5);
  c21 = c2(:,6);
  c22 = c2(:,7);

  all_tau_c = [c11 ; c12 ; c21 ; c22];

  nc = 1;

  for kk=1:numel(c11)		%-- for

  if(c11(kk) > er)
  if(c21(kk) > er)
  tcc(nc,1)=c11(kk);
  tcc(nc,2)=c21(kk);
  nc = nc + 1;
  end

  if(c22(kk) > er)
  tcc(nc,1)=c11(kk);
  tcc(nc,2)=c22(kk);
  nc = nc + 1;
  end
  end 

  if(c12(kk) > er)
  if(c21(kk) > er)
  tcc(nc,1)=c12(kk);
  tcc(nc,2)=c21(kk);
  nc = nc + 1;
  end
  
  if(c22(kk) > er)
  tcc(nc,1)=c12(kk);
  tcc(nc,2)=c22(kk);
  nc = nc + 1;
  end
  end 

  end				%-- for

  
  %------------------------------------

  X = tcc(:,1);
  Y = tcc(:,2);

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc_cc = cxy / (xst*yst) ;

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc_cc) ' N=' num2str(nd)])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/cc_corr.png'])





%------------------------------------

  X = tau;
 
  nbins = 20;
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  Pt = N;
  nt = (edges(1:end-1) + edges(2:end))/2;

  X=all_tau_c;

  nbins = 20;
  [N,edges] = histcounts(X,nbins,'Normalization', 'probability');
  Ptc = N;
  ntc = (edges(1:end-1) + edges(2:end))/2;


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])

  plot(nt,Pt,'ko-','MarkerFaceColor','black','MarkerSize',12)
  hold on
  plot(ntc,Ptc,'ro-','MarkerFaceColor','red','MarkerSize',12)
  
  set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  title(['N1 = ' num2str(ndiv) 'N2 = ' num2str(nd)])

  print(h,'-dpng','-r250',['../plot/Tau_cousin_corr_bias.png'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% del_L and growth rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  datafile=['../result/allmdc_lb.txt'];
  c1=load(datafile);

  lbc11 = c1(:,4);
  lbc12 = c1(:,5);
  lbc21 = c1(:,6);
  lbc21 = c1(:,7);


  datafile=['../result/allmdc_ld.txt'];
  c2=load(datafile);

  ldc11 = c2(:,4);
  ldc12 = c2(:,5);
  ldc21 = c2(:,6);
  ldc22 = c2(:,7);

  %------------- del_L --------------------

  dlc11 = ldc11 - lbc11;
  dlc12 = ldc12 - lbc12;
  dlc21 = ldc21 - lbc21;
  dlc22 = ldc22 - lbc22;

  %------------- growth rate --------------

  for kk = 1:numel(dlm)
  rm(kk) = (1/tm(kk))*log(ldm(kk)/lbm(kk));
  rd1(kk) = (1/td1(kk))*log(ldd1(kk)/lbd1(kk));
  rd2(kk) = (1/td2(kk))*log(ldd2(kk)/lbd2(kk));
  end



















