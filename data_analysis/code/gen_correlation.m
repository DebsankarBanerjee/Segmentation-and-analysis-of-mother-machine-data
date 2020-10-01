%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     correlation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

er = 1e-4;
maxgr=0.15;

nRep=1000;
conf_int=95;			%% change plot title if this value is changed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Self correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % data : tdiv(ndiv), tau(ndiv), dx*delL(ndiv), dx*Lbir(ndiv), dx*Ldiv(ndiv), ind_div(ndiv), C_div(ndiv)*dx
  datafile=['../../result/alltau.txt'];

  c=load(datafile);

  tau = c(:,2);
  dl = c(:,3);
  lb = c(:,4);
  ld = c(:,5);

  % growth rate calculation----------------------------------------------

  n = numel(tau);
  gr = zeros(n,1);

  for kk=1:n
  gr(kk) = (1.0/tau(kk))*log(ld(kk)/lb(kk));
  end 



%------------------------------- lb vs del_l ----------------------------

  X = lb;
  Y = ld - lb;

  [pcc]= corr_pearson(X,Y);

  %----------------------------------------
  
  [sdr, low, high] = confidence_interval(X,Y,nRep,conf_int,pcc);
  fprintf('95p confidence interval - [%f,%f] \n', low,high)
  

  nd = numel(X);
  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc) ' 95% ci=[' num2str(low) ',' num2str(high) ']'])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  xlabel('L_b','FontSize',20)
  ylabel('\Delta_L','FontSize',20)

  print(h,'-dpng','-r250',['../plot/self_lb_dl.png'])

%------------------------------------------------------------------------


%------------------------------- gr vs del_l ----------------------------

  X = gr;
  Y = ld - lb;

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc = cxy / (xst*yst) ;

  %----------------------------------------
  
  [sdr, low, high] = confidence_interval(X,Y,nRep,conf_int,pcc);
  fprintf('95p confidence interval - [%f,%f] \n', low,high)
  

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc) ' 95% ci=[' num2str(low) ',' num2str(high) ']'])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  xlabel('growth rate','FontSize',20)
  ylabel('\Delta_L','FontSize',20)

  print(h,'-dpng','-r250',['../plot/self_r_dl.png'])

%------------------------------------------------------------------------



%------------------------------- lb vs gr -------------------------------

  X = lb;
  Y = gr;

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc = cxy / (xst*yst) ;

  %----------------------------------------
  
  [sdr, low, high] = confidence_interval(X,Y,nRep,conf_int,pcc);
  fprintf('95p confidence interval - [%f,%f] \n', low,high)
  

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc) ' 95% ci=[' num2str(low) ',' num2str(high) ']'])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)
  xlabel('L_b','FontSize',20)
  ylabel('growth rate','FontSize',20)

  print(h,'-dpng','-r250',['../plot/self_lb_r.png'])

%------------------------------------------------------------------------





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Daugher-Daughter correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % data : tau_m, tau_d1, tau_d2

  datafile=['../../result/allmdd_tau.txt'];
  c1=load(datafile);

  tm = c1(:,1);
  td1 = c1(:,2);
  td2 = c1(:,3);


  % data : lb_m, lb_d1, lb_d2

  datafile=['../../result/allmdd_lb.txt'];
  c1=load(datafile);

  lbm = c1(:,1);
  lbd1 = c1(:,2);
  lbd2 = c1(:,3);


  datafile=['../../result/allmdd_ld.txt'];
  c2=load(datafile);

  ldm = c2(:,1);
  ldd1 = c2(:,2);
  ldd2 = c2(:,3);

  %------------- del_L --------------------

  dlm = ldm - lbm;
  dld1 = ldd1 - lbd1;
  dld2 = ldd2 - lbd2;

  %------------- growth rate --------------

  for kk = 1:numel(dlm)
  rm(kk) = (1/tm(kk))*log(ldm(kk)/lbm(kk));
  rd1(kk) = (1/td1(kk))*log(ldd1(kk)/lbd1(kk));
  rd2(kk) = (1/td2(kk))*log(ldd2(kk)/lbd2(kk));
  end

%-------------------------------- D1-D2 grid -------------------------------

  n = numel(tm);
  d1data = zeros(n,3);
  d2data = zeros(n,3);

  d1data = [lbd1, dld1, rd1'];
  d2data = [lbd2, dld2, rd2'];



%------------------------------- 3x3 correlation loop ----------------------------

  for k1=1:3
  for k2=1:3

  X = d1data(:,k1);
  Y = d2data(:,k2);

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc = cxy / (xst*yst) ;

  %----------------------------------------
  
  [sdr, low, high] = confidence_interval(X,Y,nRep,conf_int,pcc);
  fprintf('95p confidence interval - [%f,%f] \n', low,high)
  

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc) ' 95% ci=[' num2str(low) ',' num2str(high) ']'])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/corr_D1_' num2str(k1) '_D2_' num2str(k2) '.png'])


  end
  end

%------------------------------------------------------------------------






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mother-Daughter correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------------- M-D grid -------------------------------

  n = numel(tm);
  mdata = zeros(n,3);
  ddata = zeros(n,3);

  mdata = [lbm, dlm, rm'];

  for kk=1:n
  x=rand;
  if (x>=0.5)
  ddata(kk,:) = d1data(kk,:);
  else
  ddata(kk,:) = d2data(kk,:);
  end
  end


 % ddata = [[lbd1; lbd2], [dld1; dld2], [rd1'; rd2']];


%------------------------------- 3x3 correlation loop ----------------------------

  for k1=1:3
  for k2=1:3

  X = ddata(:,k1);
  Y = mdata(:,k2);

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc = cxy / (xst*yst) ;

  %----------------------------------------
  
  [sdr, low, high] = confidence_interval(X,Y,nRep,conf_int,pcc);
  fprintf('95p confidence interval - [%f,%f] \n', low,high)
  

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc) ' 95% ci=[' num2str(low) ',' num2str(high) ']'])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/corr_D_' num2str(k1) '_M_' num2str(k2) '.png'])






  end
  end

%------------------------------------------------------------------------





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mother-Cousin correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---------------------------------- non-unique pair data -----------------------------------

  datafile=['../../result/allmdc_lb.txt'];
  c1=load(datafile);

  lbm = c1(:,1);

  lbc11 = c1(:,4);
  lbc12 = c1(:,5);
  lbc21 = c1(:,6);
  lbc22 = c1(:,7);

% data : tau_m, tau_d1, tau_d2, tau_c11, tau_c12, tau_c21, tau_c22

  datafile=['../../result/allmdc_tau.txt'];
  c2=load(datafile);

  cm = c2(:,1);

  c11 = c2(:,4);
  c12 = c2(:,5);
  c21 = c2(:,6);
  c22 = c2(:,7);

  %all_tau_c = [c11 ; c12 ; c21 ; c22];

  datafile=['../../result/allmdc_ld.txt'];
  c3=load(datafile);

  ldm = c3(:,1);

  ldc11 = c3(:,4);
  ldc12 = c3(:,5);
  ldc21 = c3(:,6);
  ldc22 = c3(:,7);


%-----------------------------------

  n = numel(c11);
  mdata = zeros(n,3);
  cdata = zeros(n,3);
  cdata1 = zeros(n,3);
  cdata2 = zeros(n,3);
  rm = zeros(n,1);


  for kk=1:n

%------------- M ------------------
  rm(kk) = (1/cm(kk))*log(ldm(kk)/lbm(kk)); 
  dlm(kk) = ldm(kk) - lbm(kk);
  mdata(kk,:) = [lbm(kk), dlm(kk), rm(kk)];

%------------- C1 -----------------
  x=rand;
  if(c11(kk)<er || c12(kk)<er)
  if(c11(kk)>er)
  rc1 = (1.0/c11(kk))*log(ldc11(kk)/lbc11(kk));  
  cdata1(kk,:)=[lbc11(kk), ldc11(kk)-lbc11(kk), rc1 ];
  elseif(c12(kk)>er)
  rc1 = (1.0/c12(kk))*log(ldc12(kk)/lbc12(kk));  
  cdata1(kk,:)=[lbc12(kk), ldc12(kk)-lbc12(kk), rc1 ];
  end
  elseif(c11(kk)>er && c12(kk)>er && x>=0.5)
  rc1 = (1.0/c11(kk))*log(ldc11(kk)/lbc11(kk));  
  cdata1(kk,:)=[lbc11(kk), ldc11(kk)-lbc11(kk), rc1 ];
  elseif(c11(kk)>er && c12(kk)>er && x<0.5)
  rc1 = (1.0/c12(kk))*log(ldc12(kk)/lbc12(kk));  
  cdata1(kk,:)=[lbc12(kk), ldc12(kk)-lbc12(kk), rc1 ];
  end

%------------- C2 -----------------
  x=rand;
  if(c21(kk)<er || c22(kk)<er)
  if(c21(kk)>er)
  rc2 = (1.0/c21(kk))*log(ldc21(kk)/lbc21(kk));  
  cdata2(kk,:)=[lbc21(kk), ldc21(kk)-lbc21(kk), rc2 ];
  elseif(c22(kk)>er)
  rc2 = (1.0/c22(kk))*log(ldc22(kk)/lbc22(kk));  
  cdata2(kk,:)=[lbc22(kk), ldc22(kk)-lbc22(kk), rc2 ];
  end
  elseif(c21(kk)>er && c22(kk)>er && x>=0.5)
  rc2 = (1.0/c21(kk))*log(ldc21(kk)/lbc21(kk));  
  cdata2(kk,:)=[lbc21(kk), ldc21(kk)-lbc21(kk), rc2 ];
  elseif(c21(kk)>er && c22(kk)>er && x<0.5)
  rc2 = (1.0/c22(kk))*log(ldc22(kk)/lbc22(kk));  
  cdata2(kk,:)=[lbc22(kk), ldc22(kk)-lbc22(kk), rc2 ];
  end

%--------------- C ---------------  
  x=rand;
  if(x>=0.5)
  cdata(kk,:) = cdata1(kk,:);
  else
  cdata(kk,:) = cdata2(kk,:);
  end


  end %%--- for-end ---%%



%------------------------------- 3x3 correlation loop ----------------------------

  for k1=1:3
  for k2=1:3

  X = cdata1(:,k1);
  Y = cdata2(:,k2);

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc = cxy / (xst*yst) ;

  %----------------------------------------
  
  [sdr, low, high] = confidence_interval(X,Y,nRep,conf_int,pcc);
  fprintf('95p confidence interval - [%f,%f] \n', low,high)
  

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc) ' 95% ci=[' num2str(low) ',' num2str(high) ']'])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/corr_C1_' num2str(k1) '_C2_' num2str(k2) '.png'])


  end
  end

%------------------------------------------------------------------------





%------------------------------- 3x3 correlation loop ----------------------------

  for k1=1:3
  for k2=1:3

  X = cdata(:,k1);
  Y = mdata(:,k2);

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc = cxy / (xst*yst) ;

  %----------------------------------------
  
  [sdr, low, high] = confidence_interval(X,Y,nRep,conf_int,pcc);
  fprintf('95p confidence interval - [%f,%f] \n', low,high)
  

  nd = numel(X);


  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
  
  
  scatter(X,Y,25, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.5)
  %alpha 0.35

  title(['PCC = ' num2str(pcc) ' 95% ci=[' num2str(low) ',' num2str(high) ']'])

  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  print(h,'-dpng','-r250',['../plot/corr_C_' num2str(k1) '_M_' num2str(k2) '.png'])


  end
  end

%------------------------------------------------------------------------






































