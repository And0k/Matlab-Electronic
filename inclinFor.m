%use after inclinZeroCalib
str= ...
  'd:\WorkData\Cruises\_BalticSea\140301\_source\inclinometr\inclinometr.h5';
SLoad= struct('SaveDirName','c:\TEMP\', 'strDir',str, ...
  'bSaveAsInclinometr',true, 'bNotSaveTimeString',true);
    %130510 c:\TEMP\130609.h5  %SaveDirName 'f:\WorkData\inclinometr,wavegage\';
    %SLoad.SaveDirName= [SaveDir SLoad.prefix sprintf('%02d', Nprobe) '\'];
%str= {'i3', 'i4', 'i6', 'i8', 'i9', 'i11', 'i12', 'i13'};
%str= {'i4', 'i7', 'i14', 'i15'};
%str= {'i2','i3','i4','i5','i6','i7'};
str= {'i3','i4','i7','i8'}; %'i06','i11','i13',
FilterDat= struct('Time',struct('NoDataSeparator',1.1*dayHz));
%FilterDat.Time.min= datenum('19.05.2013 17:00:00', 'dd.mm.yyyy HH:MM:SS'); %'19.05.2013 17:00:00'
FilterDat.Time.max= datenum('04.03.2014 17:00:00', 'dd.mm.yyyy HH:MM:SS');
%FilterDat.smooth= [3 33];  FilterDat.decimate= 5;
FilterDat.Gx= struct('max',32768+32639,	'min',129 ,'MaxSpike', 1200);
FilterDat.Gy= FilterDat.Gx; FilterDat.Gz= FilterDat.Gx;
FilterDat.Hx= struct('max',32768+4096-1, 'min',32768-4096+1 ,'MaxSpike',100);
FilterDat.Hy= FilterDat.Hx; FilterDat.Hz= FilterDat.Hx;
FilterDat.Gsum= struct('max',1+0.11,	'min',1-0.11 ,'MaxSpike', 0.05, ...
  'MaxSpikeWidth',9, 'rep2mean', false); %[G]
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
for k= 1:numel(str)
  inclinometr2([SLoad.strDir '/' str{k}], FilterDat, SLoad);
end



%%% Old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
str= 'd:\WorkData\Experiment\naclinometr\2013\0605Lotok\_source\';
for k=11
  SCoef= listFiles(sprintf('%s#*,%d[l,]*og.TXT',str, k));
  if numel(SCoef)~=1; SCoef= listFiles(sprintf('%s#%d[l,]*og.TXT',str, k));
  end
  if numel(SCoef)==1
    [~, ~, aLog]= loadData2DB(struct('strDir', [str SCoef.name], 'col',...
      struct('inDate',1, 'DateForm',{'dd.mm.yyyy HH:MM:SS'}), 'strCmd','N'), ...
      struct('Time', struct('max',now, 'MaxSpike',1)));
  else stopHere();
  end
  Tadd= 1.2e-4; aLog.Time= aLog.Time+Tadd;
  [temp t]= max(aLog.Time);
  FilterDat.Time= struct('NoDataSeparator', 1.1*dayHz, ...
    'min', min(aLog.Time), 'max', temp+(aLog.dT(t)+300)*dayHz);
  inclinometr([str sprintf('#%d.txt',k)], '', FilterDat);
end