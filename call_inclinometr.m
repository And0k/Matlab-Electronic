%use after inclinZeroCalib
str= ...
  'd:\WorkData\Cruises\_BalticSea\140301\_source\inclinometr\inclinometr.h5';
SLoad= struct('SaveDirName',[fileparts(strrep(str,'source','subproduct')) '\'],...
   'strDir',str, 'bSaveAsInclinometr',true, 'bNotSaveTimeString',true);
    %130510 c:\TEMP\130609.h5  %SaveDirName 'f:\WorkData\inclinometr,wavegage\';
    %SLoad.SaveDirName= [SaveDir SLoad.prefix sprintf('%02d', Nprobe) '\'];
%str= {'i3', 'i4', 'i6', 'i8', 'i9', 'i11', 'i12', 'i13'};
%str= {'i4', 'i7', 'i14', 'i15'};
%str= {'i2','i3','i4','i5','i6','i7'};
str= {'i3','i4','i7'}; %,'i8', 'i06','i11','i13',
FilterDat= struct('Time',struct('NoDataSeparator',1.1*dayHz));
%FilterDat.Time.min= datenum('19.05.2013 17:00:00', 'dd.mm.yyyy HH:MM:SS'); %'19.05.2013 17:00:00'
FilterDat.Time.max= datenum('04.03.2014 17:00:00', 'dd.mm.yyyy HH:MM:SS');
%FilterDat.smooth= [3 33];  FilterDat.decimate= 5;
FilterDat.Gx= struct('max',int32(32768+32639),	'min',int32(129),'MaxSpike', int32(1200));
FilterDat.Gy= FilterDat.Gx; FilterDat.Gz= FilterDat.Gx;
FilterDat.Hx= struct('max',int32(32768+4096-1), 'min',int32(32768-4096+1),'MaxSpike', int32(100));
FilterDat.Hy= FilterDat.Hx; FilterDat.Hz= FilterDat.Hx;
FilterDat.Gsum= struct('max',1+0.11,	'min',1-0.11 ,'MaxSpike', 0.05, ...
  'MaxSpikeWidth',9, 'rep2mean', false); %[G]
FilterDat.P= struct('min',int32(1));
if false %special filtering
  FilterDat.Gx= struct('max',32768+32639,	'min',129 ,'MaxSpike',80, ...
  'MaxSpikeWidth',20, 'smooth',struct('N',1), 'rep2mean', false);
  FilterDat.Gy= FilterDat.Gx; FilterDat.Gz= FilterDat.Gx;
end
if false %special filtering version 2
  FilterDat.Gx= struct('max',32768+32639,	'min',129 ,'MaxSpike',300, ...
  'MaxSpikeWidth',9, 'smooth',struct('N',3), 'rep2mean', false);
  FilterDat.Gy= FilterDat.Gx; FilterDat.Gz= FilterDat.Gx;
end
if false %true 
  str= {'WaveGage1', 'WaveGage2', 'WaveGage3'};
  FilterDat= struct('Time',FilterDat.Time, ...
    'P',struct('max',5e5, 'min',2500, 'MaxSpike', 100, 'MaxSpikeWidth',5)); %, 'MaxSpikeWidth',3[counts]
end


%% Calc cycle
if ~exist([SLoad.SaveDirName str{numel(str)} '.mat'],'file')
  SLoad.noTXT= true;
  inclinometers= 1:numel(str);
  for k= inclinometers
    fprintf('\nCalc %s\n', str{k});
    [DATA a1D ax SLog]= inclinometr2([SLoad.strDir '/' str{k}], FilterDat, SLoad);
    SLog.fs= dayHz/linTime(DATA.Time,-3);
    %load([SLoad.SaveDirName str{inclinometers(1)} '.mat']);
    boundsT= findSeparatorConstr(DATA, SLog.fs); boundsT(end)= [];
    DATA= insertNaNs(DATA, boundsT+1); %figure(); plot(DATA.Time, DATA.Vabs);
    if isfield(DATA, 'P')
      %Instrument positions are relative to an arbitrary horizontal origin so
      %are all at the same (0,0) but must be fixed absolutely vertically from the bed.
      %% Filter P:
      bGood= ~isnan(DATA.P);
      bGoodD= bGood; bGoodD(find(bGood,1,'last'))= false;
      bInd= true(size(DATA.P)); bInd(bGoodD)= diff(DATA.P(bGood))~=0;
      bInd= ~smooth(bInd, 11)<0.01; %too smooth for water
      %figure(); plot(DATA.P, 'g'); hold on; plot(find(smooth(bInd, 11)<0.01), DATA.P(smooth(bInd,11)<0.01),'r');
      %Depth from existing pressure data to approximate the mean sensor depth:
      bGoodD= true(size(DATA.P));
      bGoodD(bGood)= ~bSpike(DATA.P(bGood,1), 1, 5);
      figure('Name','Correct P'); plot(DATA.P, 'g');
      hold on; plot(find(~bGoodD), DATA.P(~bGoodD),'.r'); %
      %stopHere('Del bad depth?'); %if no then set bGoodD(:)= true;
      DATA.P(~bGoodD)= NaN;
      plot(find(~bInd), DATA.P(~bInd), '.k');
      SLog.Depth= double(nanmean(DATA.P(bInd&bGoodD,1))) + 0.2; %depth, m of the instrument
      % Correct zero depth: dP= nanmean(DATA.P(bInd&(DATA.P<SLog.Depth/2)));
      % DATA.P= DATA.P - dP; save([SLoad.SaveDirName str{k} '.mat'], 'DATA', '-mat','-v7.3','-append');
      %bGoodD(bGood)= ~bSpike(DATA.P(bGood,1), 1, 5);
      %DATA.P(~bInd)= NaN;
    else
      SLog.Depth= 12; %depth, m of the instrument
    end
    save([SLoad.SaveDirName str{k} '.mat'], 'DATA','SLog', '-mat','-v7.3');
    %save([SLoad.SaveDirName str{k} '.mat'], 'DATA','SLog','s3D','SMout','EPout', '-mat','-v7.3');
  end
end
d2r= pi/180;
%% Calc PSD
inclinometers= 1:numel(str); 
for k= inclinometers
  clear s3D
  load([SLoad.SaveDirName str{k} '.mat']); fprintf('%s loaded\n', str{k});
  if ~exist('s3D', 'var')
    ID.fs= SLog.fs; %sampling frequency, Hz
    %Estimation parameters
    EP= struct(); %all default
    % EP.dres= 180     ;%the number of directions in the calculation)
    % EP.iter= 100     ;%the number of algorithm iterations)
    % EP.smooth= 'ON'  ;%enables smoothing of the final output)
    % EP.method= 'EMEP';%the best choice for high quality estimation with a
    %3 quantity point measurement)
    EP.nfft= 1024; %256 the number of DFTs in the spectral estimation)
    %spectral matrix structure
    SM.freqs= ((ID.fs/EP.nfft):(ID.fs/EP.nfft):(ID.fs/3))'; %0.01:0.01:0.5; %the frequency bins for the output)
    SM.dirs= -178:2:180;  %the directional bins for the output)
    SM.xaxisdir= 90; %[90] direction of the x-axis

    %   figure('Name', 'Waves'); ha1= subplot(2,1,1); ha2= subplot(2,1,2);
    %   %       displayArea(bBad, ha1, FilterDat.VControl.VabsMean.max*1.5, ...
    %   %         a1D.Time+ mean(diff(a1D.Time))/2, 'FaceColor','k', 'EdgeColor','None', ...
    %   %         'BaseValue',-FilterDat.VControl.VabsMean.max*1.5); legend(ha1,'bad');
    %   %       ADV.Time= a1D.Time;
    %   ploD(DATA, {'Vabs'}, [.5 0 0; 0 .5 .5], ha1);
    %   ploD(DATA, {'Direction'}, [.5 0 0; 0 .5 .5], ha2);
    %   linkaxes([ha1,ha2],'x')
    
    %Del. all nans (may be insrted by insertNaNs):
    %   last= numel(DATA.Time);
    %   strFields= fieldnames(DATA); strFields(strcmp(strFields,'Time'))= [];
    %   bInd= true(last,1);
    %   for cur= 1:numel(strFields); bInd= bInd&isnan(DATA.(strFields{cur})); end
    %   boundsT= [find(bInd); last+1];
    %bInd= ~bInd; DATA= structfun(@(x) x(bInd), DATA, 'UniformOutput',false); %del all NaN
    boundsT= findSeparatorConstr(DATA, SLog.fs);
    %ok if all(isnan(DATA.Vabs(boundsT))) because this separation is manually added
    last= 1; ic= 1;
    s3D= zeros(numel(SM.freqs),numel(SM.dirs),numel(boundsT)); %spectrum
    SLog.H= zeros(size(boundsT)); SLog.Tp= SLog.H; SLog.DTp=SLog.H; SLog.Dp= SLog.H;       %statistics
    %% cycle intervals
    for cur= boundsT'-1 %cur= boundsT(1);
      ind= last:cur;
      %instrument data structure
      if isfield(DATA, 'P'); bInd= ~isnan(DATA.P(ind)); end
      if isfield(DATA, 'P')&&sum(bInd)>(EP.nfft*0.8)
        ID.datatypes= {'pres' 'velx' 'vely'};   %datatype codes)
        ID.layout = [ 0     0     0;      ...   %x positions)
          0     0     0;      ...   %y positions)
          0.2   0.2   0.2];         %z positions)
        ID.data= zeros(numel(ind), numel(ID.datatypes)); %instrument data. see ID.datatypes
        ID.data(:,1)= DATA.P(ind);
        if any(~bInd); ID.data(:,1)= rep2mean(ID.data(:,1)); end
        [ID.data(:,2), ID.data(:,3)]= pol2cart(pi/2-DATA.Direction(ind)*d2r, ...
          DATA.Vabs(ind)); %Ve, Vn
      else
        ID.datatypes= {'velx' 'vely'};
        ID.layout = [ 0     0;         ...   %x positions)
          0     0;         ...   %y positions)
          0.2   0.2   ];         %z positions)
        ID.data= zeros(numel(ind), numel(ID.datatypes)); %instrument data. see ID.datatypes
        [ID.data(:,1), ID.data(:,2)]= pol2cart(pi/2-DATA.Direction(ind)*d2r, ...
          DATA.Vabs(ind)); %Ve, Vn
      end
      ID.depth= SLog.Depth;
      
      [SMout,EPout]= dirspec(ID,SM,EP,{'MESSAGE',0, 'PLOTTYPE',0}); %4 - polar plot(compass angles)
      s3D(:,:,ic)= SMout.S;
      %Hsig		Signficant wave height
      %Tp			Peak period
      %DTp		Direction of spectral peak
      %Dp			Dominant direction
      [SLog.H(ic),SLog.Tp(ic),SLog.DTp(ic),SLog.Dp(ic)]= infospec(SMout);
      
      last=cur+2; ic =ic+1;
    end
    SMout.S= [];
    save([SLoad.SaveDirName str{k} '.mat'], 'DATA','s3D','SMout','EPout','SLog','-mat','-v7.3'); %,'-append'
  end
end %s3D calc
if false %true
%Figures cycle
for k= inclinometers
  load([SLoad.SaveDirName str{k} '.mat']); fprintf('%s figures\n', str{k});
  dirPSDimg= [SLoad.SaveDirName str{k} '\PSD\'];
  if(~isdir(dirPSDimg)); mkdir(dirPSDimg); end %fprintf('dir "%s" created', dirPSDimg);
   
  bGood= squeeze(any(any(s3D))); cur= find(bGood,1); %good spectrums, and which use first
  SMout.S= real(log10(real(s3D(:,:,cur))));
  strTitle=['Directional spectrum estimated using ' EPout.method ' method, log_{10}(m^2s/deg)'];
  
  hf= figure(); h= axes('Color','none', 'NextPlot','replacechildren');
  %hT= polar(h, [0 2*pi], [SM.freqs(1) 0.8*SM.freqs(end)]); delete(hT);
  hC= plotspec(SMout, 4, h); set(h,'nextplot', 'replacechildren'); %example plot
  zmax= floor(2*log10(max(real(s3D(:)))))/2; zmin= zmax-5;
  caxis(h, [zmin zmax]); set(hC, 'LevelList', zmin:0.5:zmax); %set(hC, 'LevelListMode','auto'); get(hC, 'LevelList')
  title(h, strTitle);  
  %%
  last= 1; ic= 1; boundsT= findSeparatorConstr(DATA, SLog.fs);
  for cur= boundsT'-1 %cur= boundsT(1);
    ind= last:cur;
    if bGood(ic)
      %ind= last:cur;
      SMout.S= log10(s3D(:,:,ic));
      hC= plotspec(SMout, 4, h, get(hC, 'LevelList'), 'LineStyle', get(hC, 'LineStyle')); %4 - polar plot(compass angles)
      %title(h, strTitle);,
      if ishandle(hf)
        strFileName= datestr(DATA.Time(last),'yymmdd_HHMM');  %[datestr(DATA.Time(end)  ,'dd_HHMM')];
        set(hf, 'Name',strFileName);
        print(hf, '-djpeg',[dirPSDimg strFileName, '.jpg'], '-noui');
      else
        stopHere();
      end
    end
    last= cur+2; ic= ic +1;
  end
  %%
  close(hf);
end
end
%% Mean Velocity values
%inclinometers= 3
for k= inclinometers
  load([SLoad.SaveDirName str{k} '.mat']); fprintf('calc %s mean\n', str{k});
  ic= 1; last= 1; boundsT= findSeparatorConstr(DATA, SLog.fs);
  ind= zeros(size(boundsT)); icc= numel(ind);
  Mean= struct('Ve',ind, 'Vn',ind, 'Vabs',ind, 'Vdir',ind, 'Time', ...
    DATA.Time([1; boundsT(1:end-1)-1]), ...
    'Ve5',zeros(icc,5), 'Vn5',zeros(icc,5), 'Vabs5',zeros(icc,5), 'Vdir5',zeros(icc,5));
  
  for cur= boundsT'-1
    ind= last:cur;
    [ID.data(:,1), ID.data(:,2)]= pol2cart(pi/2-DATA.Direction(ind)*d2r, ...
      DATA.Vabs(ind)); %Ve, Vn
    Mean.Ve(ic)= nanmean(ID.data(:,1));
    Mean.Vn(ic)= nanmean(ID.data(:,2));
    [Mean.Vdir(ic), Mean.Vabs(ic)]= cart2pol(Mean.Vn(ic), Mean.Ve(ic));
    Mean.Vdir(ic)= Mean.Vdir(ic)/d2r;
    icc= 1; cLast= 1; ccur= last-cur;
    for ccur= round(linspace(1,ccur,5))
      cInd= cLast:ccur;
      Mean.Ve5(ic, icc)= nanmean(ID.data(cInd,1));
      Mean.Vn5(ic, icc)= nanmean(ID.data(cInd,2));
      [Mean.Vdir5(ic, icc), Mean.Vabs5(ic, icc)]= cart2pol( ...
        Mean.Vn5(ic, icc), Mean.Ve5(ic, icc));
      Mean.Vdir5(ic, icc)= Mean.Vdir5(ic, icc)/d2r;
      cLast= ccur+1; icc= icc+1;
    end
    last= cur+2; ic= ic +1;
  end
  %save([SLoad.SaveDirName str{k} '.mat'], 'DATA','SLog', '-mat','-v7.3');
  save([SLoad.SaveDirName str{k} '.mat'], 'DATA', 'Mean','SLog','s3D', ...
      'SMout','EPout', '-mat','-v7.3');
    %add atribute for 2D grids: vsz_range [min_x, min_y, max_x, max_y];  
end
  %'HdLPF0p08f1FIR_ls50'
%       %HdLPF10f400FIR_ls100 'HdLPF7f500FIRwbl100' 'HdLPF0p02f1FIR_wk50'
%   strFilter= 'HdLPF0p08f1FIR_ls50';
%   a1D.Ve= FilterApplyCorrectly(ID.data(:,2), strFilter, {'saveMeanTrend', 'Up&Down'}); %filter up 
%   a1D.Vn= FilterApplyCorrectly(ID.data(:,3), strFilter, {'saveMeanTrend', 'Up&Down'});
%   ax2= subplot(2,1,2);
%   a1D.Time= DATA.Time(ind);
%   figure('Name', 'SmoothVeVn'); ha1= subplot(2,1,1); ha2= subplot(2,1,2);
% %       displayArea(bBad, ha1, FilterDat.VControl.VabsMean.max*1.5, ...
% %         a1D.Time+ mean(diff(a1D.Time))/2, 'FaceColor','k', 'EdgeColor','None', ...
% %         'BaseValue',-FilterDat.VControl.VabsMean.max*1.5); legend(ha1,'bad');
% %       ADV.Time= a1D.Time;
%   ploD(a1D, {'Ve'}, [.5 0 0; 0 .5 .5], ha1);
%   ploD(a1D, {'Vn'}, [.5 0 0; 0 .5 .5], ha2);
%   linkaxes([ha1,ha2],'x')

%%
% saveTxt(rmfield(DATA, fldSecondary), [SLog.SaveDirName ...
%       SLog.strFileAdd 'IHPR']);


