%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     correlation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


er = 1e-4;
maxgr=0.5;


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


 

  %---------------------- calculate PCC ----------------------

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
  dld2 = ldd2 - lbd2;

  %-------------- curate del_L data --------

  dlm0=dlm;

  nn=numel(dlm);

  for kk = 1:nn
  if (kk>numel(dlm)) break; end
  if (dlm(kk)<=er || dld1(kk)<=er || dld2(kk)<=er)
  dlm(kk)=[];
  dld1(kk)=[];
  dld2(kk)=[];
  end
  end



  %------------- growth rate --------------

  for kk = 1:numel(dlm)
  rm(kk) = (1/tm(kk))*log(ldm(kk)/lbm(kk));
  rd1(kk) = (1/td1(kk))*log(ldd1(kk)/lbd1(kk));
  rd2(kk) = (1/td2(kk))*log(ldd2(kk)/lbd2(kk));
  end


  %-------------- curate gr data --------

  idum=1;
  icur = [];

  for kk = 1:numel(rm)
  if (rm(kk)>maxgr || rd1(kk)>maxgr || rd2(kk)>maxgr)
  icur(idum)=kk;
  idum=idum+1;
  end
  end

  for kk=1:numel(icur)
  rm(icur(kk))=[];
  rd1(icur(kk))=[];
  rd2(icur(kk))=[];
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %------------ del L correlations ---------------------

  %---------------------- calculate PCC ----------------------

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

  X = [rm'; rm'];
  Y = [rd1'; rd2'];

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

  
  %---------------------- calculate PCC ----------------------

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

  print(h,'-dpng','-r250',['../plot/cc_tau_corr.png'])





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
  lbc22 = c1(:,7);


  datafile=['../result/allmdc_ld.txt'];
  c2=load(datafile);

  ldc11 = c2(:,4);
  ldc12 = c2(:,5);
  ldc21 = c2(:,6);
  ldc22 = c2(:,7);

  %------------- check cousin existence -------

  nc = 1;

  for kk=1:numel(c11)		%----------------------- for

  if(lbc11(kk) > er && ldc11(kk) > er)		%%------ c11 if

  if(lbc21(kk) > er && ldc21(kk) > er)
  dlcc(nc,1)=ldc11(kk)-lbc11(kk);
  dlcc(nc,2)=ldc21(kk)-lbc21(kk);
  grcc(nc,1)=(1/c11(kk))*log( ldc11(kk)/lbc11(kk) );
  grcc(nc,2)=(1/c21(kk))*log( ldc21(kk)/lbc21(kk) );
  nc = nc + 1;
  end

  if(lbc22(kk) > er && ldc22(kk) > er)
  dlcc(nc,1)=ldc11(kk)-lbc11(kk);
  dlcc(nc,2)=ldc22(kk)-lbc22(kk);
  grcc(nc,1)=(1/c11(kk))*log( ldc11(kk)/lbc11(kk) );
  grcc(nc,2)=(1/c22(kk))*log( ldc22(kk)/lbc22(kk) );
  nc = nc + 1;
  end

  end 						%%------ c11 end

  if(lbc12(kk) > er && ldc12(kk) > er)		%%------ c12 if

  if(lbc21(kk) > er && ldc21(kk) > er)
  dlcc(nc,1)=ldc12(kk)-lbc12(kk);
  dlcc(nc,2)=ldc21(kk)-lbc21(kk);
  grcc(nc,1)=(1/c12(kk))*log( ldc12(kk)/lbc12(kk) );
  grcc(nc,2)=(1/c21(kk))*log( ldc21(kk)/lbc21(kk) );
  nc = nc + 1;
  end
  
  if(lbc22(kk) > er && ldc22(kk) > er)
  dlcc(nc,1)=ldc12(kk)-lbc12(kk);
  dlcc(nc,2)=ldc22(kk)-lbc22(kk);
  grcc(nc,1)=(1/c12(kk))*log( ldc12(kk)/lbc12(kk) );
  grcc(nc,2)=(1/c22(kk))*log( ldc22(kk)/lbc22(kk) );
  nc = nc + 1;
  end

  end						%%------ c12 end 

  end				%----------------------- for

  %---------------------- curate del_L and gr data ----------

  [nn, dum]=size(dlcc);

  for kk = 1:nn
  if (dlcc(kk,1)<=er || dlcc(kk,2)<=er)
  dlcc(kk,:)=[];
  end
  end


  [nn, dum]=size(grcc);

  for kk = 1:nn
  if (grcc(kk,1)>maxgr || grcc(kk,2)>maxgr)
  grcc(kk,:)=[];
  end
  end



  %---------------------- calculate PCC ----------------------

  X = dlcc(:,1);
  Y = dlcc(:,2);

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

  print(h,'-dpng','-r250',['../plot/cc_del_L_corr.png'])


  %%----------------------------------------------------------


  X = grcc(:,1);
  Y = grcc(:,2);

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

  print(h,'-dpng','-r250',['../plot/cc_gr_corr.png'])
