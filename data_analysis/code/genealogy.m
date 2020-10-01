



  %------------------------  Genealogy  --------------------------
  
  gmat = zeros(Nbranch*gsize,2^gen);
  gtau = zeros(Nbranch*gsize,2^gen);
  
  for i=1:Nbranch		%------- branch loop
  gen2flag=0;
  % gen - 0
  gmat(1+(i-1)*gsize,1) = M0(i);
  if(dumcid(M0(i),2)>er) gtau(1+(i-1)*gsize,1) = dumcid(M0(i),3)-dumcid(M0(i),2); end
  
  % gen - 1
  k = find(mdd(:,1)==M0(i));
  D01 = mdd(k,2);
  D02 = mdd(k,3); 
  gmat(1+(i-1)*gsize + 1,1) = D01;
  gmat(1+(i-1)*gsize + 1,2) = D02;
  
  gtau(1+(i-1)*gsize + 1,1) = dumcid(D01,3)-dumcid(D01,2);
  gtau(1+(i-1)*gsize + 1,2) = dumcid(D02,3)-dumcid(D02,2);
  
  % gen - 2
  k = find(mdd(:,1)==D01);
  if ~isempty(k)
  D11 = mdd(k,2);
  D12 = mdd(k,3); 
  gmat(1+(i-1)*gsize + 2,1) = D11;
  gmat(1+(i-1)*gsize + 2,2) = D12;
  
  gtau(1+(i-1)*gsize + 2,1) = dumcid(D11,3)-dumcid(D11,2);
  gtau(1+(i-1)*gsize + 2,2) = dumcid(D12,3)-dumcid(D12,2);
  
  else
  gen2flag = gen2flag + 1;
  end
  
  k = find(mdd(:,1)==D02);
  if ~isempty(k)
  D13 = mdd(k,2);
  D14 = mdd(k,3); 
  gmat(1+(i-1)*gsize + 2,3) = D13;
  gmat(1+(i-1)*gsize + 2,4) = D14;
  
  gtau(1+(i-1)*gsize + 2,3) = dumcid(D13,3)-dumcid(D13,2);
  gtau(1+(i-1)*gsize + 2,4) = dumcid(D14,3)-dumcid(D14,2);
  
  else
  gen2flag = gen2flag + 1;
  end
  
  if(gen2flag > 1+er)		%---- gen-3 loop 
  
  continue
  
  else
  
  % gen - 3
  k = find(mdd(:,1)==D11);
  if ~isempty(k)
  D21 = mdd(k,2);
  D22 = mdd(k,3); 
  gmat(1+(i-1)*gsize + 3,1) = D21;
  gmat(1+(i-1)*gsize + 3,2) = D22;
  
  gtau(1+(i-1)*gsize + 3,1) = dumcid(D21,3)-dumcid(D21,2);
  gtau(1+(i-1)*gsize + 3,2) = dumcid(D22,3)-dumcid(D22,2);
  
  end
  
  k = find(mdd(:,1)==D12);
  if ~isempty(k)
  D23 = mdd(k,2);
  D24 = mdd(k,3); 
  gmat(1+(i-1)*gsize + 3,3) = D23;
  gmat(1+(i-1)*gsize + 3,4) = D24;
  
  gtau(1+(i-1)*gsize + 3,3) = dumcid(D23,3)-dumcid(D23,2);
  gtau(1+(i-1)*gsize + 3,4) = dumcid(D24,3)-dumcid(D24,2);
  
  end
  
  k = find(mdd(:,1)==D13);
  if ~isempty(k)
  D25 = mdd(k,2);
  D26 = mdd(k,3); 
  gmat(1+(i-1)*gsize + 3,5) = D25;
  gmat(1+(i-1)*gsize + 3,6) = D26;
  
  gtau(1+(i-1)*gsize + 3,5) = dumcid(D25,3)-dumcid(D25,2);
  gtau(1+(i-1)*gsize + 3,6) = dumcid(D26,3)-dumcid(D26,2);
  
  end
  
  k = find(mdd(:,1)==D14);
  if ~isempty(k)
  D27 = mdd(k,2);
  D28 = mdd(k,3); 
  gmat(1+(i-1)*gsize + 3,7) = D27;
  gmat(1+(i-1)*gsize + 3,8) = D28;
  
  gtau(1+(i-1)*gsize + 3,7) = dumcid(D27,3)-dumcid(D27,2);
  gtau(1+(i-1)*gsize + 3,8) = dumcid(D28,3)-dumcid(D28,2);
  
  end
  
  end			%---- gen-3 loop 
  
  end			%------- branch loop
