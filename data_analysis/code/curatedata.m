




  %-------------- curate del_L data --------

  dlmdd=[dlm, dld1, dld2];

  nn=numel(dlm);
  idum = 0;

  for kk = 1:nn
  if (dlmdd(kk-idum,1)<er || dlmdd(kk-idum,2)<er || dlmdd(kk-idum,3)<er)
  dlmdd(kk-idum,:)=[];
  idum = idum + 1;
  end
  end

  dlm=[]; dld1=[]; dld2=[]; 

  dlm=dlmdd(:,1); dld1=dlmdd(:,2); dld2=dlmdd(:,3);

  %%-----------------------------------------




  %-------------- curate gr data --------

  rmdd=[rm', rd1', rd2'];

  nn=numel(rm);
  idum = 0;

  for kk = 1:nn
  if (rmdd(kk-idum,1)>maxgr || rmdd(kk-idum,2)>maxgr || rmdd(kk-idum,3)>maxgr)
  rmdd(kk-idum,:)=[];
  idum = idum + 1;
  end
  end

  rm=[]; rd1=[]; rd2=[]; 

  rm=rmdd(:,1); rd1=rmdd(:,2); rd2=rmdd(:,3);






  %---------------------- curate del_L and gr data ----------

  [nn, dum]=size(dlcc);
  idum = 0;

  for kk = 1:nn
  if (dlcc(kk-idum,1)<=er || dlcc(kk-idum,2)<=er)
  dlcc(kk-idum,:)=[];
  idum = idum + 1;
  end
  end


  [nn, dum]=size(grcc);
  idum = 0;

  for kk = 1:nn
  if (grcc(kk-idum,1)>maxgr || grcc(kk-idum,2)>maxgr)
  grcc(kk-idum,:)=[];
  idum = idum + 1;
  end
  end




  %%-----------------------------------------
