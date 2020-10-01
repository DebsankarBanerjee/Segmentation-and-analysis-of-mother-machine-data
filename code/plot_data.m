function []=plot_data(k)

er=1e-6;

fid = fopen('dataname');
txt = textscan(fid,'%s','delimiter','\n');
c = txt{1} ;

datafile=['../data/cell_id_', char(c), '.txt'];

cdata = load(datafile);

tb=cdata(:,2);
td=cdata(:,3);

lb=cdata(:,4);
ld=cdata(:,5);

[nn, dum]=size(cdata);

kdum=1;

for kk=1:nn
if(td(kk)>er && tb(kk)>er)
tau(kdum) = td(kk) - tb(kk);
kdum = kdum+1;
end
end


kdum=1;

for kk=1:nn
if(ld(kk)>er && lb(kk)>er)
dl(kdum) = ld(kk) - lb(kk);
kdum = kdum+1;
end
end

xi = [1:kdum-1];

size(xi)
size(tau)

  h = figure('vis', 'off');
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])

  subplot(1,2,1)
  plot(xi,tau,'ko','MarkerFaceColor','black','MarkerSize',8)
  hold on
  plot(xi,0*xi,'red-','LineWidth',2)
  %set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  title(['division time at t=' num2str(k)])


  subplot(1,2,2)
  plot(xi,dl,'ko','MarkerFaceColor','black','MarkerSize',8)
  hold on
  plot(xi,0*xi,'red-','LineWidth',2)
  %set(gca, 'YScale', 'log')
  xt = get(gca, 'XTick');
  set(gca, 'FontSize', 16)
  yt = get(gca, 'YTick');
  set(gca, 'FontSize', 16)

  title(['added length at t=' num2str(k)])

  print(h,'-dpng','-r100',['../plot/plot_',char(c),'_tau_dL','.png'])

