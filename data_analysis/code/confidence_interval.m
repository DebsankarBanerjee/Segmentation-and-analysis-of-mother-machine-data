function [pcc, low, high] = confidence_interval(X,Y,nRep,conf_int,rho)

%-------------------------------------------------------------

  n=numel(X);
  pcc=zeros(nRep,1);

  xydata=[X, Y];

  for k=1:nRep

  %-------------- resample data -------------
  resample = datasample(xydata,n);
  
  XX = resample(:,1);			%datasample(X,n);
  YY = resample(:,2);			%datasample(Y,n);

  %-------- PCC from resampled data ---------

  xms = XX-mean(XX);
  yms = YY-mean(YY);

  cxy = mean(xms.*yms);

  xst = std(XX);
  yst = std(YY);

  pcc(k) = cxy / (xst*yst);

  end

  drho = pcc - rho;

  sdr = sort(drho);

  idlow = ceil( ((100-conf_int)/2.0)*(nRep/100) );
  idhgh = ceil( (conf_int + (100-conf_int)/2.0)*(nRep/100) );

  prsn=rho
  
  high = rho - sdr(idlow);
  low = rho - sdr(idhgh);
