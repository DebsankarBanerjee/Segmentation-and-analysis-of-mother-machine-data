function [pcc]= corr_pearson(X,Y)

%-----------------------------------------------------------

  xms = X-mean(X);
  yms = Y-mean(Y);

  cxy = mean(xms.*yms);

  xst = std(X);
  yst = std(Y);

  pcc = cxy / (xst*yst) ;

