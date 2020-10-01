  clc;
  clear all;
  
  % cell_id_array : current_index / time@Birth / time@Division / length@Birth / length@Division / parent ID / index at division / position at division
  
  er=1E-4;
  ndiv = 0;
  dx = 50/310;
  dt = 1;

  error_msg = 0;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%% get file names %%%%%%%%%%%%%%%%%%
  
  fid = fopen('dataname');
  txt = textscan(fid,'%s','delimiter','\n');
  c = txt{1} ;
  
  
  ndata=numel(c)
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%% LOOP ON DATA (TIF FILES) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for in=1:ndata

  datafile=['../../data/', char(c(in)),'.txt'];				% ,'.txt'
  
  fid10=fopen(['../../result/alltau_',char(c(in)),'.txt'],'w');

  fid21=fopen(['../../result/mdd_tau_',char(c(in)),'.txt'],'w');
  fid22=fopen(['../../result/mdd_lb_',char(c(in)),'.txt'],'w');
  fid23=fopen(['../../result/mdd_ld_',char(c(in)),'.txt'],'w');

  fid31=fopen(['../../result/mdc_tau_',char(c(in)),'.txt'],'w');
  fid32=fopen(['../../result/mdc_lb_',char(c(in)),'.txt'],'w');
  fid33=fopen(['../../result/mdc_ld_',char(c(in)),'.txt'],'w');  

  fid99=fopen(['../../result/x_error_',char(c(in)),'.txt'],'w'); 

  cid=load(datafile);
  dumcid=cid;
  
  [nr, nc] = size(cid);
  
  mdd = zeros(nr,3);

  mdc = zeros(nr,7);
  
  %---------------- current position ---------------
  cp = cid(:,1);
  
  %---------------- division & growth rate -----------------
  
  tb = cid(:,2); 
  td = cid(:,3);
  
  lb = cid(:,4); 
  ld = cid(:,5);
  
  for i=1:nr
  
   dum = tb(i) * td(i);
   if(dum > er )
   ndiv = ndiv + 1;
   tau(ndiv) = td(i) - tb(i);
   tdiv(ndiv) = td(i);
   delL(ndiv) = ld(i) - lb(i);
   Lbir(ndiv) = lb(i);
   Ldiv(ndiv) = ld(i);
   ind_div(ndiv) = cid(i,7);
   C_div(ndiv) = cid(i,8);


   fprintf(fid10,'%f %f %f %f %f %f %f\n', tdiv(ndiv), tau(ndiv), dx*delL(ndiv), dx*Lbir(ndiv), dx*Ldiv(ndiv), ind_div(ndiv), C_div(ndiv)*dx );
   end
  
  end
  
  %----------------------Initial cell no---------------------
  %------------------------- Branch no ----------------------
  Ninit = 0;
  Nbranch = 0;
  M0=[];
  for i=1:nr 
  if(tb(i)<er)
  Ninit = Ninit + 1;
  
  if(td(i)>er)
  Nbranch = Nbranch + 1;
  M0 = [M0 ; i];
  end
  
  end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % creating the mother-daughter database from cell data
  % cell ids of mother: imo, daughters: id1 and id2
  %----------------- Mother daughter : MDD unit -------------------
  
  gflag = 0;
  npair = 0;
  npoint = nr;
  
  while npoint >= 1		% ------- while loop over all cell
  
  dum = dumcid(npoint,2);
  
  %------- check for birth -------
  
  if(dum>er)	%--------- main if loop
   
   id1 = npoint;
  
  %-------------------- Mo, D1, D2 id -----------------------
   % assign mother cell id
   imo = dumcid(id1,6);
  
   % assign second daughter
   if(round(imo-dumcid(id1-1,6)) < er)
   id2 = id1-1;
   npoint = npoint-2;
   npair = npair + 1;
   else
   fprintf('cell div error')
   dum
   npoint = npoint -1;
   end
  
   mdd(npair,:) = [imo id1 id2];
  
  else
  
   npoint = npoint -1;
  
  end   		%--------- main if loop
  
  end		% ------- while loop over all cell
  
  %-------------- cure mother daughter pair data by getting rid of the
  % pairs where both daughter division was not captured
  %--------------- get rid of empty rows --------------------
  
   dum1 = mdd(:,1);
   dum1(dum1==0)=[];
   dum2 = mdd(:,2);
   dum2(dum2==0)=[];
   dum3 = mdd(:,3);
   dum3(dum3==0)=[];
   mdd = [];
   mdd=[dum1 dum2 dum3];
  
   mdd = sortrows(mdd, 1);

   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  [npair, dum]=size(mdd);

  for kk=1:npair 	% npair loop

  if(cid(mdd(kk,1),2) > er && cid(mdd(kk,2),3) > er && cid(mdd(kk,3),3) > er)
  %---- div time ----
  tm = cid(mdd(kk,1),3) - cid(mdd(kk,1),2);
  td1 = cid(mdd(kk,2),3) - cid(mdd(kk,2),2);
  td2 = cid(mdd(kk,3),3) - cid(mdd(kk,3),2);

  fprintf(fid21,'%f %f %f\n', tm*dt, td1*dt, td2*dt );
  
  %---- l@birth ----
  lmb = cid(mdd(kk,1),4);
  ld1b = cid(mdd(kk,2),4);
  ld2b = cid(mdd(kk,3),4);

  fprintf(fid22,'%f %f %f\n', lmb*dx, ld1b*dx, ld2b*dx );

  %---- l@div ----
  lmd = cid(mdd(kk,1),5);
  ld1d = cid(mdd(kk,2),5);
  ld2d = cid(mdd(kk,3),5);

  fprintf(fid23,'%f %f %f\n', lmd*dx, ld1d*dx, ld2d*dx );

   %%---------- print error -----------
   if (lmd-lmb < er || ld1d-ld1b < er || ld2d-ld2b < er )
   error_msg = 1;
   fprintf(fid99,'%d \n', kk );
   end


  end		% if-end

  end			% npair loop



   %%---------- print error -----------
   if (error_msg > er)
   for kk=1:npair
   fprintf(fid99,'%d %f %f %f\n', kk, mdd(kk,1), mdd(kk,2), mdd(kk,3) );
   end
   end




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %----------------- Mother daughter cousin : MDC unit -------------------
  nmdc = 0;

  for kk=1:npair		%---------- for every mdd pair find cousins

  c11=0;
  c12=0;
  c21=0;
  c22=0;

  %----------- D1, C11, C12 ------------
  d = mdd(kk,2);
  cc = find(cid(:,6)==d);
  if(~isempty(cc))
  c11 = cc(1);
  c12 = cc(2);
  end

  [d c11 c12]

  %----------- D2, C21, C22 ------------
  d = mdd(kk,3);
  cc = find(cid(:,6)==d);
  if(~isempty(cc))
  c21 = cc(1);
  c22 = cc(2);
  end

  [d c21 c22]

  %--------------------------------------

  %if(isempty(c11)) c11=0; end
  %if(isempty(c12)) c12=0; end
  %if(isempty(c21)) c21=0; end
  %if(isempty(c22)) c22=0; end

  %---------- make mdc -----------------
  dum1 = c11+c12;
  dum2 = c21+c22;
  if ((dum1 > er) && (dum2 > er))
  nmdc = nmdc + 1;

  mdc(nmdc,1) = mdd(kk,1);
  
  mdc(nmdc,2) = mdd(kk,2);
  mdc(nmdc,3) = mdd(kk,3);
  
  mdc(nmdc,4) = c11;
  mdc(nmdc,5) = c12;
  mdc(nmdc,6) = c21;
  mdc(nmdc,7) = c22;
  end
  

  end 			%---------- npair for loop

%-------------------------------------------------

  dum1 = mdc(:,1);
  dum1(dum1==0)=[];

  dum2 = mdc(:,2);
  dum2(dum2==0)=[];

  dum3 = mdc(:,3);
  dum3(dum3==0)=[];

  dum4 = mdc(:,4);
  dum4(dum4==0)=[];

  dum5 = mdc(:,5);
  dum5(dum5==0)=[];

  dum6 = mdc(:,6);
  dum6(dum6==0)=[];

  dum7 = mdc(:,7);
  dum7(dum7==0)=[];  
  
  
  
  mdc = [];
  mdc = [dum1 dum2 dum3 dum4 dum5 dum6 dum7];
  mdc = sortrows(mdc, 1);

  [nmdc dum] = size(mdc);

%%%%%%%%%-----------------------------------------

  for kk=1:nmdc		%%--- nmdc loop
  
  %---- div time ----

  tm = cid(mdc(kk,1),3) - cid(mdc(kk,1),2);

  td1 = cid(mdc(kk,2),3) - cid(mdc(kk,2),2);
  td2 = cid(mdc(kk,3),3) - cid(mdc(kk,3),2);

  tc11 = cid(mdc(kk,4),3) - cid(mdc(kk,4),2);
  tc12 = cid(mdc(kk,5),3) - cid(mdc(kk,5),2);
  tc21 = cid(mdc(kk,6),3) - cid(mdc(kk,6),2);
  tc22 = cid(mdc(kk,7),3) - cid(mdc(kk,7),2);

  m = [tc11, tc12, tc21, tc22];

  for jj=1:4
  if (m(jj)<er)
  m(jj) = 0;
  end 
  end %jj

  %------- length at birth and division --------

  %---- l@birth ----
  lmb = cid(mdc(kk,1),4);
  ld1b = cid(mdc(kk,2),4);
  ld2b = cid(mdc(kk,3),4);

  lc11b = cid(mdc(kk,4),4);
  lc12b = cid(mdc(kk,5),4);
  lc21b = cid(mdc(kk,6),4);
  lc22b = cid(mdc(kk,7),4);

  %---- l@div ----
  lmd = cid(mdc(kk,1),5);
  ld1d = cid(mdc(kk,2),5);
  ld2d = cid(mdc(kk,3),5);

  lc11d = cid(mdc(kk,4),5);
  lc12d = cid(mdc(kk,5),5);
  lc21d = cid(mdc(kk,6),5);
  lc22d = cid(mdc(kk,7),5);



  if(cid(mdc(kk,1),2) > er && cid(mdc(kk,2),3) > er && cid(mdc(kk,3),3) > er && (m(1)+m(2)>er) && (m(3)+m(4)>er) )

  fprintf(fid31,'%f %f %f %f %f %f %f\n', tm*dt, td1*dt, td2*dt, m(1)*dt, m(2)*dt, m(3)*dt, m(4)*dt);

  fprintf(fid32,'%f %f %f %f %f %f %f\n', lmb*dx, ld1b*dx, ld2b*dx, lc11b*dx, lc12b*dx, lc21b*dx, lc22b*dx);

  fprintf(fid33,'%f %f %f %f %f %f %f\n', lmd*dx, ld1d*dx, ld2d*dx, lc11d*dx, lc12d*dx, lc21d*dx, lc22d*dx);

  end		% if-end

  end			%%--- nmdc loop
  

  

  
  
 
  end			%%%%%%%%%%%% end data loop
  
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
