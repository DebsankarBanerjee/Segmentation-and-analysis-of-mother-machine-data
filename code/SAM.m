%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	MAIN CODE FOR SEGMENTATION AND ANALYSIS OF MOTHER MACHINE DATA (SAM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; 
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% global parameters etc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = 50/310 ;		% pixel length in micro meter
dt = 1 ;		% frame rate in min

eps = 1E-6;		% small number : numerical tolerence

IwantPlot = 1;		% if frame-by-frame plot is not needed then put IwantPlot=1 ortherwise 0
plot_start_time = 0 ;
plot_end_time = 100;


sigma = 1;		% stn div of the gauss filter
gs = 3*sigma;		% size of gauss filter window = 3 x sigma

med_size = 3;		% size of median filter window = med_size x med_size pix

Rsharp = 4;		% Radius value used in imsharp
Asharp = 6;		% Amount of sharpness enhansement

R_imopen = 20;		% radius for imopen in background (BG) preparation


thr = 1500;		% threshold value to binarize image ~ 1500 
coffthr = 0.45;		% threshold for cutoff position along-y

%grad_thr = 90000;
%avgthr = 1.0;			%1.0;
%min_depth = 0.4;

pk_dist = 8;		% min peak distance tolerance in FindPeaks
pk_prom = 0.5;		% min peak prominence tolerance in FindPeaks ~ 0.45

divf = 0.65;		% division length reduction parameter, L_new < divf*L_old, used in prevention of cell div loss

fil_length = 30;	% cell size >= fil_length is considered for possible filamentation 
fil_div = 10;		% when cell_L > 40 pix, asymmetric cell div : L0 - L > fil_div, then div counted.
agf = 1.50;		% abnormal growth parameter, L_new > agf*L_old, used in prevention of cell div loss
min_lb = 10;		% minimal length at birth, to supress over segmentation
max_dbline = 15;	% max allowed change in new BLine(bottom line), while calculating from l_sum
max_dl_tot = 8;		% max change in total length of cells : used to check cell intrusion into MF-channel

bwn = 4; 		%connectivity for bw_conn_comp
mincellsize = 20;	%no of pixel, any smaller object is to be removed as noise/abbaration

min_std = 5;		% min value of std of intensity profile, used to determine if an
			% object in CC contains cells or not.

max_cell_no = 150;	% no of rows of cell_id matrix, to be fixed at some opt size

pdelta = 3;		% pix no fluctuation in level set, due to determination of level set from relative
			% quantities like intensity ratio, mean intensity at t, etc. << new cell size

max_Xshift = 8;		% max distance of object center_X from the image center (imcenter) to be allowed
			% used to remove random cell/object error from dynamic thresholding

Theta_cutoff = 30;	% objects/cells oriented with the channel are selected, any object with
			% orientation less than Theta_cutoff is removed

%end_centroid = 12;	% if cell centroid_y is within 'end centriod' range from open end then cell is extracted

%bottom_line = 20;

end_ch_high = 8000;	% high values representing clogging cells at channel end. Used to calculate bottom line. 

min_area_ratio = 0.05	% min area/ convex_area ratio allowed for cells 



%%%%%%%%%%%%%%% open files %%%%%%%%%%%%%%%
fid1=fopen('../data/divtime.txt','w');


%%%%%%%%%%%%%%% initialize %%%%%%%%%%%%%%%
cell_size = [];
divcount = 1;


%======================== DATA HANDLING MODULE =================================

%%%%%%%%%%%%%% get file names %%%%%%%%%%%%%%%%%
fid = fopen('dataname');
txt = textscan(fid,'%s','delimiter','\n');
c = txt{1} ;

fprintf('number of channels to be analysed \n')
ndata=numel(c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% LOOP ON DATA (TIF FILES) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i=1:ndata

%%%%%%%%%%%%%%% define parameters %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% open channel specific files here %%%%%%%%%%%%%%%

fid10=fopen(['../data/cell_id_',char(c(i)),'.txt'],'w');
fid20=fopen(['../data/cellone_time_',char(c(i)),'.txt'],'w');


%-----------------------------------------------------------------------------
%%%%%%%%%%%%%% extracting images %%%%%%%%%%%%%
  FileTif=['../im/', char(c(i)),'.tif'];
  InfoImage=imfinfo(FileTif);
  mImage=InfoImage(1).Width;
  nImage=InfoImage(1).Height;
  NumberImages=length(InfoImage);

  FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
  TifLink = Tiff(FileTif, 'r');
%-----------------------------------------------------------------------------
fprintf('Number of frames in this tif file : %d', NumberImages)
%-----------------------------------------------------------------------------
%%%%%%%%%%%%% image reading loop %%%%%%%%%%%%%%

for k=1:NumberImages
  TifLink.setDirectory(k);
  FinalImage(:,:,k)=TifLink.read();
   
end			
%------------------------------------------------------------- Image loop ends
  TifLink.close();
%-----------------------------------------------------------------------------

%===================== DATA HANDLING MODULE : END ============================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  n1 = nImage ;		% row = y-axis
  n2 = mImage ;		% column = x-axis

  IM = zeros(n1,n2);
  im = zeros(n1,n2);
  imth=zeros(n1,n2);
  imbw=zeros(n1,n2);


fprintf('tif images extracted : OK \n')
fprintf('**************************************************\n')
fprintf('data = %s\n',char(c(i)))
fprintf('**************************************************\n')


%%%%%%%%%%%%%%% initialization for every channel %%%%%%%%%%%%%%%

cell_id = zeros(max_cell_no,8);
ncid = 0;
L0 = [];
LRaw0 = [];
ma_L0 = [];
cent_y0 = [];
l_end0 = 0;
cy_end0 = 0;
BLine = 0;
calBLineFlag = 0;
endCellDiv = 0;
l_sum = 0;
Lcsum0 = 0;
miny0 = [];


%%%%%%%%%%%%% Calibrate bottom line %%%%%%%%%%%%%%

dum = [];
ndum = 1;

for k=1:NumberImages

  IM = FinalImage(:,:,k);
  x_av = mean(IM,2);
  high_x = x_av > end_ch_high ;

  for kk=1:numel(high_x)
  if ( high_x(kk) > eps )
  dum(ndum)=kk;
  ndum = ndum + 1;
  break;
  end
  end
   
end

if (~isempty(dum))
bottom_line = min(dum)
else
bottom_line = n1 - 3*pdelta
end





%>>>>>>>>>>>>>>>>>>>>>>>>>>  TIME LOOP  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% ANALYSIS ON STACK  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%% initialize %%%%%%%%%%%%%%%

flagLastCell=0;
flagLastCellDiv=0;

for k=1:NumberImages              % --------------------------- time loop

  %--------------------------------------------------- IM = current frame
   IM = FinalImage(:,:,k);
   time=(k-1)*dt;
   

fprintf(' >>>>>>>>>>> time = %d\n',k)
fprintf('**************************************************\n')

%%%%%%%%%%%%%%% open frame specific files here if needed %%%%%%%%%%%%%%%


%fid100=fopen(['../data/whatever_channel_no_',num2str(i,'%03d'),'_frame_no_',num2str(k,'%03d'),'.txt'],'w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% ANALYSIS ON EACH FRAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





 
  %%%%%%%%%%%%%%%%%% IMAGE ENHANCEMENT %%%%%%%%%%%%%%%
   
  %--------------------------------------------------- BG substraction 
  %[IM1] = end_channel_clear(IM,bottom_line,n1,k);
  background = imopen(IM,strel('disk',R_imopen));
  IM = IM - background;
  
  %im = medfilt2(IM, [med_size med_size]);
  
  %------------------ Gaussian filter ----------------   >>> IMPORTANT <<<
  im = imfilter(double(IM),fspecial('gaussian',[gs gs], sigma));
  
  
  %------------------ Sharpening ---------------------   >>> IMPORTANT <<<
  im = imsharpen(im,'Radius',Rsharp,'Amount',Asharp);
  
  
  %%%%%%%%%%%%%%%%%%% Thresholding  %%%%%%%%%%%%%%%%%%
  %--------------------------------------------------    >>> IMPORTANT <<<
  imth=im;
  imth(im < thr) = 0;					% using fixed threshold


  %%%%%%%%%%%%%%%%%%%%% Clear end channel %%%%%%%%%%%%%%%%%%%%
  % everything below bottomline is removed

  [imth] = end_channel_clear(imth,bottom_line,n1,k);	


  %%%%%%%% Projection along channel length %%%%%%%%%%%
  % create an 1D average intensity profile of the channel 

  x_av=mean(imth,2);
  y_av=mean(imth,1);
  x_av = x_av/mean(x_av);
  y_av = y_av/mean(y_av);

  
  [m1,m2]=cutoff_pos(x_av,coffthr);
  [m3,m4]=cutoff_pos(y_av,coffthr);

  x_mid = find(y_av==max(y_av));				

  % correct for global up/down shift in detection 
  %---- extracted part contain image ~ from 1st pole to last pole ------- 

  if (m3 - pdelta > 1)
  xr_1 = m3 - pdelta;
  else
  xr_1 = m3;
  end

  if (m4 + pdelta < n2)
  xr_2 = m4 + pdelta;
  else
  xr_2 = m4;
  end 

  x_av=mean(imth(:,xr_1:xr_2),2);
  meanxav = mean(x_av);
  ch_y = x_av/mean(x_av);	
  
  %----------------------------------------------------------------------
   




  %%%%%%%%%%%%%%%%%%% Local max and min of average column  %%%%%%%%%%%
  
  % average intensity profile is inverted and up-shifted to find intensity minima as a peaks

  %%--------------------------------------------------------

  negch_y = -ch_y + max(ch_y) + 0.1;
  [dps,dpy] = findpeaks(negch_y,'MinPeakDistance',pk_dist,'MinPeakProminence',pk_prom);	

  miny = dpy				% position of the minima

  %%--------------------------------------------------------




  %%%%%%%%%%%%%%%%% 1st Segmentation to create binary image %%%%%%%%%%%
  % after segmentation each connected component in CC represents a cell

  [imbw] = channel_segment(miny,imth,imbw,mincellsize,eps,m1);

  imbw0=imbw;

  [nlength,Lstart,Lend,Lcsum] = length_calculation(imbw,x_mid);


  CC = bwconncomp(imbw,bwn);


  %%%%%%%%%%%%%%%%%%% Adavptive minima search  %%%%%%%%%%%%%%%%

  

  [imth] = segmentation_cleanup(CC,imth,min_std,Theta_cutoff,max_Xshift,x_mid,min_area_ratio,n1);

  %[miny] = division_loss(CC,miny,nlength,LRaw0,agf,divf,min_lb,m1,m2);

  [miny] = min_search(imth,xr_1,xr_2,pk_dist,pk_prom);

  [imbw] = channel_segment(miny,imth,imbw,mincellsize,eps,m1);

  [nlength,Lstart,Lend,Lcsum] = length_calculation(imbw,x_mid);

  LRaw = nlength;


  [miny] = adaptive_min_search(imth,CC,miny,miny0,negch_y,agf,divf,pdelta,nlength,Lcsum0,Lcsum,pk_prom,pk_dist,xr_1,xr_2,meanxav,k,n1);
           
  
  [imbw] = channel_segment(miny,imth,imbw,mincellsize,eps,m1);
  [nlength,Lstart,Lend,Lcsum] = length_calculation(imbw,x_mid);


  CC = bwconncomp(imbw,bwn);
  
  s  = regionprops(CC,'Centroid');
  centroids = cat(1, s.Centroid);
  cent_y=centroids(:,2);
  cent_y=sort(cent_y);

  %%%%%%%%%%%%%%%% check end cell division and cell intrusion %%%%%%%%%%%%%%%

  [endCellDiv, endExt, cellIntrusion] = end_cell_check(nlength,Lstart,Lend,cent_y,L0,l_end0,divf,pdelta,max_dl_tot,k);  


  %%%%%%%%%%%%%%%% Remove cell intrusion %%%%%%%%%%%%%%%

  if (cellIntrusion > eps)
  [imth] = end_channel_clear(imth,Lend(end-1),n1,k);
  end

  %%%%%%%%%%%%%%%% Re-calibrate bottom line %%%%%%%%%%%%%%%

  if ( cent_y(end) <= n1 - BLine )
  calBLineFlag = 0;
  fprintf('BLine -> bottom_line. cY= %d, line = %d \n', cent_y(end), BLine)
  end

  if ( calBLineFlag > eps )
  L_check = n1 - BLine;
  else
  L_check = Lend(end);
  end

  fprintf('check for boundary cell : Lcheck = %d, calBLineFlag = %d \n', L_check, calBLineFlag)

  if( (L_check + pdelta >= bottom_line)  && calBLineFlag < eps  || ( calBLineFlag > eps && cent_y(end) > n1 - BLine ) )

  fprintf('New BLine calculation \n')

  if ( calBLineFlag < eps  &&  endCellDiv < eps )

  BLine = n1 - Lstart(end);
  
  elseif ( calBLineFlag < eps  &&  endCellDiv > eps )

  if (endExt < eps)
  BLine = n1 - Lstart(end-1);
  else
  BLine = n1 - Lstart(end);
  end

  end

  [imth,BLine,calBLineFlag,l_sum] = determine_BLine(imth,Lend,Lstart,BLine,pdelta,bottom_line,calBLineFlag,l_sum,max_dbline,endCellDiv,n1,k);

  end



  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%% Final binarization and length calculation %%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  [imbw] = channel_segment(miny,imth,imbw,mincellsize,eps,m1);
  
  [nlength,Lstart,Lend,Lcsum] = length_calculation(imbw,x_mid)

  %------- y-position of last cell centroid --------

  CC2 = bwconncomp(imbw,bwn)
  ncell = CC2.NumObjects 
  s  = regionprops(CC2,'Centroid');
  centroids2 = cat(1, s.Centroid);
  cent_y=centroids2(:,2);
  cent_y=sort(cent_y)

  cLast = max(cent_y);	

  %-------- major axis length calculation ----------

  [ma_L] = major_axis_length(CC2)	
  
  %------------------ lower pole position of the last cell -----------------------
  
  l_end = Lend(end);
  cy_end = cent_y(end);								 

  fprintf(fid20,'%f %f %f\n',time, dx*nlength(1), dx*ma_L(1) );

  %-------------------------------------------------------------------------------
  %%%%%%%%%%%%%%% plot segmented images %%%%%%%%%%%%%%%
 
  if(IwantPlot > eps)
  if(k >= plot_start_time && k <= plot_end_time)
  plot_channel(IM,imbw,imbw0,centroids2,bottom_line,BLine,n1,n2,k,NumberImages,i,c);
  end
  end


%-------------------------------------------------------------------------------------------------



  %----------------- for first frame : store L0 then stop the loop and continue here ------------- 

  if(k==1)
  L0 = nlength;
  cLast0 = cLast;
  cell_id = initiate_cell(cell_id, ncell);
  ncid = ncell;
  l_end0 = l_end;
  cent_y0 = cent_y;
  cy_end0 = cy_end;
  ncell0 = ncell;
  LRaw0 = LRaw;
  miny0 = miny;
  Lcsum0 = Lcsum;
  ma_L0=ma_L;
  continue
  end



  %===============================================================================================
  flagCellDiv=0;		% flag to indicate cell division
  flagCellExt=0;		% flag to indicate cell extrusion
  
  Ndiv=0;
  Next=0;


  nlength
  L0


  %-----------------------------------------------------------------------------------------------
  %%%%%%%%%%%%%%%%% Cases of cell division %%%%%%%%%%%%%%
  
  %--------------- check for cell division --------------
  Lmin = min(nlength);
  Lmin0 = min(L0);

  if(numel(nlength) > numel(L0))		% to determine normal cell div: number increases
    flagCellDiv=1;
  else
    for kk=1:numel(nlength)			% div + extrusion: number may remain same		         
    if(nlength(kk) < divf*L0(kk)) 
    flagCellDiv=1;
    break;
    end
    end
  end


  %fprintf('new and old min length = %d %d\n', Lmin, Lmin0)
  %fprintf('new and old cell no = %d %d\n', numel(nlength), numel(L0))
  % fprintf('deivcheck2 = %d\n',flagCellDiv)
  % fprintf('deivcheck = %d\n',flagCellDiv)
  %------------------------------------------------------
  
  %--------------- cell div has happened ----------------
  if(flagCellDiv > eps)
  
  fprintf('cell div at time = %d\n',k)
  
  [cell_id, ncid, Ndiv] = cell_division(nlength, L0, ma_L, ma_L0, cell_id, time, cent_y0, ncell, ncid, divf,fil_length,fil_div, pdelta);

  end



  %-----------------------------------------------------------------------------------------------
  %%%%%%%%%%%%%%%%% Cases of cell extrusion %%%%%%%%%%%%%%

  %--------------- check for cell extrusion --------------
  if(cLast0 > cLast + pdelta && Ndiv > eps) flagCellExt=1; end

  %if(l_end < l_end0 - pdelta) flagCellExt=1; end

  if(ncell-Ndiv < ncell0) flagCellExt=1; end

  %-------------------------------------------------------

  %--------------- cell extrusion has happened ----------------
  if(flagCellExt > eps)
 
  fprintf('cell extrusion at time = %d\n',k)
  
  [cell_id, Next] = cell_extrusion(cent_y, cent_y0, cell_id, ncell, Next, Ndiv,pdelta);

  end





%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%
  


%%%%%%%%%%%%%%%%%%% Stop in case of error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% balance cell number from division and extrusion:

dum_cn_1 = ncell;
dum_cn_2 = ncell0 + Ndiv - Next;

% balance from cell_id structure:

max_cell_index = max(cell_id(:,1));

if (abs(dum_cn_1-dum_cn_2) > eps || abs(ncell - max_cell_index) > eps)

fid90=fopen(['../data/cell_error_',char(c(i)),'.txt'],'w');
for kk=1:max_cell_no
fprintf(fid90,'%f %f %f %f %f %f %f %f\n',cell_id(kk,1),cell_id(kk,2),cell_id(kk,3),cell_id(kk,4),cell_id(kk,5),cell_id(kk,6),cell_id(kk,7),cell_id(kk,8) );
end

fprintf('cell count error at time = %d\n',k)

break;

end  


%==============================>>> STORE THE CURRENT VARIABLES FOR NEXT STEP <<<=====================================

  L0=[];
  LRaw0=[];
  ma_L0=[];
  Lcsum0=[];

  L0 = nlength
  cLast0 = cLast;
  l_end0 = l_end;
  cent_y0 = cent_y;
  cy_end0 = cy_end;
  ncell0 = ncell;
  LRaw0 = LRaw;
  Lcsum0 = Lcsum;
  miny0 = miny;
  ma_L0=ma_L;



  
end     %---------------------------- time loop on frames END

%==============================>>> WRITE THE CELL EVOLUTION DATA FOR ONE CHANNEL IN FILE <<<=====================================

[cidl,dum]=size(cell_id)
for kk=1:cidl
if (cell_id(kk,2)*cell_id(kk,3) > eps)
tau(divcount) = cell_id(kk,3) - cell_id(kk,2);
divcount = divcount + 1;
end
fprintf(fid10,'%f %f %f %f %f %f %f %f\n',cell_id(kk,1),cell_id(kk,2),cell_id(kk,3),cell_id(kk,4),cell_id(kk,5),cell_id(kk,6),cell_id(kk,7),cell_id(kk,8) );
end


end	%---------------------------- loop on data ENDs



for kk=1:numel(tau)
fprintf(fid1,'%f\n', tau(kk) );
end

mean_divtime = mean(tau);

plot_channel(IM,imbw,imbw0,centroids2,bottom_line,BLine,n1,n2,k,NumberImages,i,c);

plot_data(k);				%% only use with script file run

fprintf('number of detected cell div, mean div time = %d %f\n',divcount,mean_divtime)
fprintf('frames analyzed = %d, total frames = %d \n',k,NumberImages)
fprintf('---------------- END OF CODE --------------\n')






 


