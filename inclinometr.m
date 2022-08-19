function [DATA a1D ax SLoad]= inclinometr(DataDirMask, txtFileCoef, FilterDat, SParam)
%"inclinometr(DataDirMask, txtFileCoef)" loads data from DataDirMask, 
%upply coefficients from txtFileCoef and save to DataDir
% tab-delimited file with columns: 
%"Time    	Heading     	Peach     	Roll     Date     Time "
%
%SParam: false => Not to display figure,
%   [] => Not to save time string column (less memory demand)
%   struct with fields Time and raw Gx Gy Gz Hx Hy Hz) => Not to load data
%nargin=1 => use coefficients [0 1]
%If provided 4-th argument saves to DataBase instead of tab-delimited file
%4-th argument is struct: 
%DB.name - system name of database to connect
%DB.tbl  - table name

% Examples
% txtFileCoef= 'd:\WorkData\Cruises\_BalticSea\130510\_source\inclinometr\coef#7.txt';
% inclinometr( ...
%   'd:\WorkData\Cruises\_BalticSea\130510\_source\inclinometr\#7.txt', ...
%   txtFileCoef, struct('Time',struct('NoDataSeparator', 1.1*dayHz), ...
%   'smooth',[3 33], 'decimate',5));

%DataDirMask= 'd:\WorkData\Experiment\naclinometr\COMPAS\130430.txt';
%txtFileCoef= 'd:\WorkData\Experiment\naclinometr\coef130506.txt';
%inclinometr(DataDirMask, txtFileCoef, true);
% Not to display figure: inclinometr(DataDirMask, txtFileCoef, false);
% Not save time string column: DATA= inclinometr(DataDirMask,txtFileCoef,[]);
% Not to load data: DATA= inclinometr(DataDirToSave, txtFileCoef, a1D);
% Saving to DB:
%DB= struct('name','BlS1110', 'tbl','Tchain processed');
% %old: inclinometr(DataDirMask, txtFileCoef, true, DB);

if nargin<3||isempty(FilterDat); FilterDat= struct(); end
if ~isfield(FilterDat, 'Time'); FilterDat.Time= struct(); end
if ~isfield(FilterDat.Time, 'max')
  FilterDat.Time.max= now;
end
if ~isfield(FilterDat.Time, 'min')
  FilterDat.Time.min= datenum( '20.12.1977', 'dd.mm.yyyy');
end
if ~isfield(FilterDat.Time, 'MaxSpike')
  FilterDat.Time.MaxSpike= dayHz;
end
if ~isfield(FilterDat.Time, 'NoDataSeparator')
  FilterDat.Time.NoDataSeparator= [];
end
SLoad= struct('strDir',DataDirMask, 'col',struct('FirstRow',0,  ...
  'inDate','', 'DateForm',{'yyyy mm dd HH MM SS'}, ...
  'header'       ,  'yyyy mm dd HH MM SS P Power Gx Gz Gy Hx Hy Hz' , ...
  'headerReplace', ['yyyy mm dd HH MM SS Gx Gy Gz Hx Hy Hz' , ...
  char(13), 'Time Date1 Date2 Date3 Date4 Date5 Gx Gy Gz Hx Hy Hz'], ...
  'ColSeparator',char(9), 'NotWarningNewFields',{'Time'}, ...
  'valueFormat','u'), 'strCmd','D', ...
  'DB', struct('name', '.h5')); %
str= '';
if (nargin<4)||(~isfield(SParam, 'a1D'))
  if strcmp(DataDirMask((end-2):end), 'mat'); load(DataDirMask); b= false;
  else  
    str= [SLoad.strDir(1:end-4) 'incl_src.mat'];
    b= ~exist(str ,'file');
    if ~b; load(str); end               
  end
  if b
    [SLoad SLog a1D]= loadData2DB(SLoad, struct('Time', FilterDat.Time));
     str= [SLoad.DataDirName bIf(numel(SLoad.FilesInfo)~=1, '', ...
      SLoad.FilesInfo.name(1:end-4)) 'incl_src.mat'];    
    save(str, 'a1D','SLog','SLoad');
  end
  if ~strcmp(SLoad.strDir, DataDirMask) %file moved
    if b
      SLoad.strDir= [DataDirMask(1:(strfind(DataDirMask,'incl_src.mat')-1)) '.txt'];
    else SLoad.strDir= DataDirMask;
    end
    SLoad.FilesInfo.name= SLoad.strDir;
    SLoad= loadPrepare(rmfield(SLoad, {'col','SaveDirName','DataDirName','strFileAdd'}));
    fprintf(' in renamed dir ');
  else
    p= strfind(SLoad.FilesInfo.name, '#');
    if any(p)&&strcmp(SLoad.strFileAdd,'_') %for old mat
      SLoad.strFileAdd= [SLoad.FilesInfo.name(p(end):(end-4)) SLoad.strFileAdd];
    end
    fprintf('\nAccording to loaded mat-file settings result will be saved to\n%s\n', ...
      fileparts( SLoad.SaveDirName(1:end-1)));
  end
  
  a1D= FSFilter(a1D, FilterDat, 'Time',0,     [],[], [],SLoad);           %rem bad in Time
  str= fieldnames(FilterDat); b= strncmp(str, 'H',1);
  a1D= FSFilter(a1D, FilterDat, str(b),0, 'Time',[], [],SLoad);           %set NaNs in H
  str(strcmp(str,'Time')|b)= [];
  a1D= FSFilter(a1D, FilterDat, str,[], 'Time','Time', [],SLoad); %interp Bad in other fields
  %if nargin<2; DATA= []; ax= []; return; end
else
  a1D= SParam.a1D;
end
str= fieldnames(a1D);
k= find(strcmp(str, 'Time'));
a1D= orderfields(a1D, [k 1:(k-1) (k+1):numel(str)]);
%datestr(a1D.iTimeDbl([1 end]));
if ~issorted(a1D.Time)
  warning('Check:LoadedTime', 'Time not sorted, sorting...')
  [a1D.Time, ind]= sort(a1D.Time); ind= int32(ind);
  for k=2:numel(str)
    a1D.(str{k})= a1D.(str{k})(ind);
  end  
  bInd= diff(a1D.Time)==0;
  if any(bInd)
    ind= [];
    for k=1:numel(str)
      if any(a1D.Gx(bInd)~=a1D.Gx([false; bInd]))
        ind= find(a1D.Gx(bInd)~=a1D.Gx([false; bInd]));
        disp(ind)
      end
    end
    if any(ind)
      stopHere('Duplicate Time %d pairs with not equal data... Delete?', sum(bInd))
    end
    for k=1:numel(str)
      a1D.(str{k})(bInd)= [];
    end
  end
end
if (nargin<2)||isempty(txtFileCoef)
  SCoef= struct();
  if(SLoad.strFileAdd(1)=='#')&&exist('d:\Work\MatlabHistory\mat\coef_Electronics.mat', 'file')
    p= str2double(regexp(SLoad.strFileAdd(2:end), '\d+', 'match','once'));       
    SCoef= load('d:\Work\MatlabHistory\mat\coef_Electronics.mat', 'Inclinometer');
    if ~isfield(SCoef.Inclinometer,'i'); stopHere('Bad coefficients file!');
    end
    bInd= [SCoef.Inclinometer.i]==p;
    k= sum(bInd);
    if k==0
      p= str2double(input(sprintf('No coef for probe#%d. if use other # input it:', p),'s'));
      bInd= [SCoef.Inclinometer.i]==p;
      k= sum(bInd);
    end
    if k>1
      for k= find(bInd);
        bInd(k)= ~((a1D.Time(1)< SCoef.Inclinometer(k).TimeRange(1))||...
                   (a1D.Time(1)>=SCoef.Inclinometer(k).TimeRange(2))); %invesion for right NaN proc 
      end
      if sum(bInd)==1;
        k= k(bInd); bInd= false(size(SCoef.Inclinometer)); bInd(k)= true;
      end
      k= sum(bInd);
    end
    if k==1
      fprintf('\nUse coefficients for probe #%d obtained %s ', p, ...
        datestr(SCoef.Inclinometer(bInd).TimeProcessed, 'dd.mm.yyyy HH:MM'));
      SCoef= SCoef.Inclinometer(bInd);
    end
  end
else
  if ~any(txtFileCoef=='\')
    if strcmp(txtFileCoef, '*')&&(SLoad.strFileAdd(1)=='#')
      SCoef.strDir= [fileparts(DataDirMask) '\coef' SLoad.strFileAdd(1:end-1) '*.??t']; %mat or txt
    else
      SCoef.strDir= [fileparts(DataDirMask) '\' txtFileCoef];
    end
    SCoef= loadPrepare(SCoef);
  else
    SCoef.strDir= txtFileCoef;
  end
  if strcmp(txtFileCoef((end-3):end), 'mat')
    SCoef= load(txtFileCoef);
  else
    SCoef.col= struct('FirstRow',1, 'ColSeparator',char(9));
    [~,~,sc]= loadData2DB(SCoef);
    b= strncmp(sc.Probe, {'G'},1);
    if any(b)
      SCoef.G.A= [sc.A1x(b) sc.A2x(b) sc.A3x(b)];
      SCoef.G.C= sc.C(b);
    end
    bInd= strncmp(sc.Probe, {'H'},1);
    if any(bInd)
      SCoef.H.A= [sc.A1x(bInd) sc.A2x(bInd) sc.A3x(bInd)];
      SCoef.H.C= sc.C(bInd);
    end
    b= b|bInd; %other coef
    str= sc.Probe(~b);
    if ~isempty(str)
      sc= cell2mat(struct2cell(rmfield(sc, 'Probe'))'); sc(b,:)= [];
      for k=1:numel(str)
        SCoef.(str{k})= fliplr(sc(k,:));
      end
      f= fittype(sprintf('poly%d', numel(SCoef.Vabs)-find(SCoef.Vabs, 1, 'first')));
      SCoef.Vabs= cfit(f,cell(SCoef.Vabs));
    end
  end
end
%% correct int32 NaN values 
str= {'Gx', 'Gy', 'Gz', 'Hx', 'Hy', 'Hz'};
for k= 1:6
  if ~isfloat(a1D.(str{k}))
    bInd= a1D.(str{k})==int32(0);
    a1D.(str{k})= single(a1D.(str{k}));
    a1D.(str{k})(bInd)= NaN;
  end
end
if(~isfield(SCoef, 'G'))||any(isnan(SCoef.G.A(:)))
  warning('inclinometr:noCoef', 'no acceletometr coeficients. set to [0 1]');
  SCoef.G.A= 1;
  SCoef.G.C= zeros(3,1);
end
if(~isfield(SCoef, 'H'))||any(isnan(SCoef.H.A(:)))
  warning('inclinometr:noCoef', 'no magnetometr coeficients. Set to [0 1]');
  SCoef.H.A= 1;
  SCoef.H.C= zeros(3,1);
else %% Correct bad H component
  str= {'Hx', 'Hy', 'Hz'};
  for k= [3 1 2] %process Z at first
    bInd= isnan(a1D.(str{k}));
    if any(bInd)
      bGood= [true true true]; bGood(k)= false; iStEn= find(bGood); 
      a1D.(str{k})= single(a1D.(str{k})); temp= sum((SCoef.H.A(bGood, bGood)* ...
        ([single(a1D.(str{iStEn(1)})(bInd)), single(a1D.(str{iStEn(2)})(bInd))]' ...
        - repmat(SCoef.H.C(bGood),1, sum(bInd)))).^2);
      temp(temp>1)= 1; temp= sqrt(1 - temp)/SCoef.H.A(k,k);     
      a1D.(str{k})= rep2mean(a1D.(str{k}), bInd); %interpolate sign
      bGood= a1D.(str{k})(bInd)<SCoef.H.C(k);
      temp(bGood)= -temp(bGood);
      a1D.(str{k})(bInd)= SCoef.H.C(k) + temp;
      a1D.(str{k})= rep2mean(a1D.(str{k}));       %interpolate remained NaNs
    end
    
  end
end
if(~isfield(SCoef, 'Vabs'))||(isfloat(SCoef.Vabs)&&any(isnan(SCoef.Vabs)))
  warning('inclinometr:noCoef', 'no Vabs coeficients. Set to [0 1]');
  SCoef.Vabs= @(x) x;
end


%for k=1:numel(str); sc.(str{k})(b)= []; end

N= int32(size(a1D.Time, 1));
%% Work in blocks
iStEn= save2DB_TimeSections(a1D.Time, [1 1]);%, Separator
LastSep= int32(size(iStEn,1));
if LastSep>1; fprintf(' break in %d time intervals:', LastSep)
else          fprintf('...');
end
% a1D.Gxyz= NaN(N,3,'single');
% a1D.Hxyz= NaN(N,3,'single');
str= fieldnames(a1D); str(strcmp(str, 'Time'))= [];
for p= 1:LastSep
  fprintf(' %2d,', p);
  bInd= false(N, 1);
  bInd(iStEn(p,1):iStEn(p,2))= true;
  n= iStEn(p,2) - iStEn(p,1) + 1;
  DATA= struct('Time',[], 'Vabs',[], 'Direction',[], ...
    'Inclination',[], 'Heading',[], 'Pitch',[], 'Roll',[]);  
  % Apply accelerometer coef
  try
    Gxyz= (SCoef.G.A*([single(a1D.Gx(bInd)), single(a1D.Gy(bInd)), ...
                single(a1D.Gz(bInd))]' - repmat(SCoef.G.C, 1, n)))';
  catch ME1
    if strcmp(ME1.identifier, 'MATLAB:nomem')
      %     fprintf(1, '\nNeed to reduce size of block. Ok. Continuing...');
      %     ind= int32(numel(fieldnames(DATA))); ind= [ind 1:(ind-1)];
      %     SLoad.col.NRows= int32(FileSize)/2 + 2;
      %     clear Etemp
      stopHere();
      continue;
    else
      stopHere();
    end
  end
  % Apply compas coef
  Hxyz= (SCoef.H.A*([single(a1D.Hx(bInd)), single(a1D.Hy(bInd)), ...
              single(a1D.Hz(bInd))]' - repmat(SCoef.H.C ,1, n)))';
% DATA= forApplyCoefficients(rmfield(a1D,{'Gx','Gy','Gz','Hx','Hy','Hz'}), ...
%   cell2mat(struct2cell(rmfield(sc, 'Probe'))'), sc.Probe', sc);
%DATA.Time= DATA.iTimeDbl; DATA= rmfield(DATA,'iTimeDbl'); a1D.P fliplr()
%DATA.Direction  = single(wrapToPi(atan2(tan(be.elevation),tan(-be.bank)) - hdg.heading)*180/pi);
 


%% Filter
% for k=1:numel(str) %filter bad
%   a1D.(str{k})(a1D.(str{k})==0)= NaN;
% end
temp= sqrt(sum(Gxyz.^2,2)); %+a1D.Gy.^2+a1D.Gz.^2
% if numel(SCoef.G.A)>1
%   t= nanmean(temp); %1
%   bBad= (t-0.3 > temp)|(temp > t + 0.3);
% else
%   bBad= bSingleSpikes(temp, 300);
% end
%[false(iStEn(p,1),1) bBad false(N-iStEn(p,2),1)] a1D.Gy(bBad)= NaN; a1D.Gz(bBad)= NaN; rep2mean(bBad)

separator= int32(findSeparator(a1D.Time(bInd), FilterDat.Time.NoDataSeparator)); 
for k=1:3
  %Gxyz(:,k)= rep2mean(Gxyz(:,k), bBad, a1D.Time(bInd));
  if isfield(FilterDat, 'smooth')
    Gxyz(:,k)= separateDo(@(x,ind) medfiltSmooth(x(ind), FilterDat.smooth), ...
      Gxyz(:,k), separator);
    Hxyz(:,k)= separateDo(@(x,ind) medfiltSmooth(x(ind), FilterDat.smooth), ...
      Hxyz(:,k), separator);
  end
end
if isfield(FilterDat, 'decimate')
  ind= 1:FilterDat.decimate:n;
  Gxyz= Gxyz(ind,:);
  Hxyz= Hxyz(ind,:);
  DATA.Time= a1D.Time(iStEn(p,1):FilterDat.decimate:iStEn(p,2)); % DATA.Time= DATA.Time(ind);
  separator= separator/FilterDat.decimate;
else DATA.Time= a1D.Time(bInd);
end
%   temp= separateDo(@(x,ind) FilterApplyCorrectly(x(ind), ...
%     'HdLPFIR0p1f1mf10',{'saveMeanTrend','Up&Down'}), a1D.P, ...
%     reshape([iMin; iMax],1, numel(iMin)*2));
%    0.443*fd/w
%0.443*10/33

% Direction calc
[DATA.Pitch DATA.Roll]= clcPitchRoll(Gxyz);
DATA.Heading= -clcHeading(Hxyz, DATA.Pitch, DATA.Roll, false);
DATA.Direction= wrapTo180or360((atan2(tan(DATA.Roll),tan(DATA.Pitch)) + ...
  DATA.Heading)*180/pi);
% DATA.Direction  = wrapToPi(atan((...
% (tan(DATA.Roll)./tan(DATA.Pitch)).*(1+tan(DATA.Roll).^2)+tan(DATA.Pitch).*tan(DATA.Roll)...
% )./sqrt(1+tan(DATA.Pitch).^2+tan(DATA.Roll).^2)))*180/pi; %- DATA.Heading
%DATA.Inclination= atan(sqrt(tan(DATA.Pitch).^2+tan(-DATA.Roll).^2));
DATA.Inclination= atan2(sqrt(sum(Gxyz(:,[true true false]).^2,2)),Gxyz(:,3));

DATA.Vabs= NaN(size(DATA.Inclination), 'single');
bInd= DATA.Inclination<pi/2;
DATA.Vabs(bInd)= SCoef.Vabs(1./sqrt(tan(pi/2 - DATA.Inclination(bInd)))); %polyval();

DATA.Inclination= DATA.Inclination*180/pi; %wrapTo360();
DATA.Heading    = DATA.Heading*180/pi;
DATA.Pitch      = DATA.Pitch  *180/pi;
DATA.Roll       =-DATA.Roll   *180/pi;
%% Speed calc

  if nargin<4||~isfield(SParam, 'bNotShowFigure')
% be = calcBankElevation(Gxyz(:,2), Gxyz(:,3), Gxyz(:,1));
% hdg= calcHeading(DATA.Roll, DATA.Pitch, a1D.Hy,Hxyz(:,3),a1D.Hx, false);
% DATA.Direction= single(-DATA.Heading*180/pi);
% 
% DATA.Inclination= single(DATA.Roll*180/pi);
% DATA.Heading    = single(-DATA.Heading*180/pi);
% DATA.Pitch      = single(DATA.Pitch*180/pi);
% DATA.Roll       = single(-DATA.Roll*180/pi);
% DATA.Direction= single(wrapToPi(atan2(tan((-DATA.Roll)*pi/180),tan((DATA.Pitch)*pi/180)) - ...
%   -DATA.Heading*pi/180)*180/pi);
    strID= [datestr(DATA.Time(1), 'yymmdd_HHMM') SLoad.strFileAdd];
    ax= get(0,'ScreenSize'); ax= ax.*[(ax(3)*0.11)/ax(1) 1 0.89 0.9];
    figure('Name', sprintf('%d %s: source and result', p, strID), ...
      'Position', ax);
    ax(1)= subplot(4,1,1); hold on
    plot(ax(1), DATA.Time,sqrt(sum(Hxyz.^2,2)),'y', DATA.Time,Hxyz(:,1), ...
    'r' , DATA.Time,Hxyz(:,2),'g',  DATA.Time,Hxyz(:,3),'b');
    legend({'|H|','Hx','Hy','Hz'}, 'Location','NorthEastOutside');  legend boxoff
    set(ax(1), 'yLim', [-1.2 1.2]); datetick(ax(1), 'x', 'HH:MM'); grid on;
    
    ax(2)= subplot(4,1,2); hold on; plot(ax(2), DATA.Time, temp, 'y'); 
    plot(ax(2), DATA.Time,Gxyz(:,1),'r', DATA.Time,Gxyz(:,2),'g', ...
    DATA.Time,Gxyz(:,3),'b', DATA.Time,DATA.Vabs,'k');
    legend({'|G|','Gx','Gy','Gz','Vabs'}, 'Location','NorthEastOutside');  legend boxoff
    set(ax(2), 'yLim',[-1.2 1.2]); datetick(ax(2), 'x', 'MM:SS', 'keeplimits'); grid on;
    
    ax(3)= subplot(4,1,3); hold on
    plot(ax(3), DATA.Time, DATA.Inclination, 'r');
    plot(ax(3), DATA.Time, DATA.Direction,   'c');
    plot(ax(3), DATA.Time, DATA.Heading,     'k');
    temp= get(ax(3),'Ylim');
    legend(ax(3), {'Incl','Dir','Heading'}, 'Location','NorthEastOutside');  legend boxoff
    set(ax(3),'yLim', [max(-180,temp(1)) min(180,temp(2))]); 
    datetick(ax(3), 'x', 'MM', 'keeplimits'); grid on;
   
    ax(4)= subplot(4,1,4); hold on
    plot(ax(4), DATA.Time,DATA.Pitch,'r', DATA.Time,DATA.Roll,'b');
    legend(ax(4), {'Pitch','Roll'}, 'Location','NorthEastOutside'); legend boxoff
    datetick(ax(4),'x','HH:MM'); linkaxes(ax, 'x'); 
    set(ax(1),'xLim', DATA.Time([1,end])); grid on;
    temp= cell2mat(get(ax,'Position'));
    temp(:, 3)= max(temp(:, 3)); temp(:, 4)= max(temp(:, 4))*1.3;
    arrayfun(@(p) set(ax(p),'Position',temp(p, :)), 1:4);
    %% Save figure
    try print('-dpng', [SLoad.DataDirName strID '.png']);
    catch ME
      fprintf('can not save fig');
    end
    try if LastSep>1; close(get(ax(1),'Parent')); end
    catch ME
    end
  end

if nargin>3&&isfield(SParam, 'DB')
%   conn= database(DB.name, '', '');
%   save2DB(conn, DB.tbl, a1D);  
%   %str= fieldnames(a1D)'; str= [str {'Update'}]; %(~strcmp(str,'iTimeDbl'))
%   %save2DB(conn, DB.tbl, a1D, {'iTimeDbl' 'Depth'}, str);      
%   close(conn);
  DB.tbl= struct('D1', SParam.DB.tbl, 'log','Tchain sections log');
  
  Track.TimeRange= SLoad.TimeRange;
  Track.processedTime= datestr(now(), '#yyyy/mm/dd HH:MM:SS#');
  t= find(SLoad.strDir=='\');
  for k= numel(t):-1:2
    if ~strcmp(SLoad.strDir((t(k-1)+1):t(k)), {'source\','subproduct\'})
      t(1:(k-2))= []; break;
    end
  end
  Track.File= SLoad.FilesInfo(1).name;
  Track.iName= sprintf('%s%02.f',SLoad.strDir((t(1)+1):t(2)),0);
  Track.Comment= ['P=' mat2str(P)];
  Track.bType= 's';
  save2DBL(struct('DB',DB), DB.tbl, a1D, [], Track, fieldnames(a1D)); 
  %%
  stopHere('Continue ("F5") to save netCDF?')
  if diff(DATA.Time([1, end])) > 30
    strFile= [datestr(DATA.Time(1),'yymmdd-') ...
      datestr(DATA.Time(end),'mmdd')];
  else
    strFile= [datestr(DATA.Time(1),'yymmdd_HHMM-') ...
      datestr(DATA.Time(end)  ,'dd_HHMM')];
  end
  strFile= [SLoad.SaveDirName strFile '.nc'];
  k= 0;
  k=k+1; strCFnames{k}= 'time';        strCFunits{k}= 'days since 1970-01-01 00:00:00';
  k=k+1; strCFnames{k}= 'depth';       strCFunits{k}= 'm';
  k=k+1; strCFnames{k}= 'temperature'; strCFunits{k}= 'degree_Celsius';
  %% Write
  k= netcdf.getConstant('NC_GLOBAL');
  ncid= netcdf.create(strFile,  netcdf.getConstant('CLOBBER') ...
    ); %netcdf.getConstant('NETCDF4')+netcdf.getConstant('CLASSIC_MODEL')|
  netcdf.putAtt(ncid,k, 'title','Black sea inclinometr');
  netcdf.putAtt(ncid,k, 'institution', ['SIO RAS project\n , ' ...
    'Device produced  and data processed at AB SIO RAS, Kaliningrad, Russia']);
  %netcdf.putAtt(ncid,k, 'source' ,'Teledyne RDI DVS + CTD Depth');
  netcdf.putAtt(ncid,k, 'desc' ,'data from chain of termo sensors and pressure sensor');
  netcdf.putAtt(ncid,k, 'featureType' ,'timeSeriesProfile');
  
IDdimStation= netcdf.defDim(ncid, 'N_PROFILE', 1);
IDdimDepth=   netcdf.defDim(ncid, 'NCells', nT);
IDdimRecord=  netcdf.defDim(ncid, 'N_DATE_TIME', netcdf.getConstant('NC_UNLIMITED')); %or L
IDvarStation(1)= netcdf.defVar(ncid, 'latitude' , 'NC_DOUBLE',IDdimStation);
netcdf.putAtt(ncid,IDvarStation(1),  'units'    , 'degrees_north');
netcdf.putAtt(ncid,IDvarStation(1),  'long_name', 'station latitude');
IDvarStation(2)= netcdf.defVar(ncid, 'longitude', 'NC_DOUBLE',IDdimStation);
netcdf.putAtt(ncid,IDvarStation(2),  'units'    , 'degrees_east');
netcdf.putAtt(ncid,IDvarStation(2),  'long_name', 'station longitude');
IDvarRecord= NaN(1,numel(str));

for k=1:numel(str)
  if strcmp(str{k}, {'Time'});      
    strType= 'NC_DOUBLE';
  else
    strType= 'NC_FLOAT';
  end
  if strcmp(str{k}, {'Temp'})
    IDvarRecord(k)= netcdf.defVar(ncid, str{k}, strType, [IDdimDepth, IDdimRecord]);
  else
    IDvarRecord(k)= netcdf.defVar(ncid, str{k}, strType, IDdimRecord);
  end
  if k<numel(strCFunits)&&~isempty(strCFunits{k})
    netcdf.putAtt(ncid,IDvarRecord(k), 'units', strCFunits{k});
    if ~isempty(strCFnames{k})
      netcdf.putAtt(ncid,IDvarRecord(k), 'long_name', strCFnames{k});
    end
  end
end
netcdf.endDef(ncid);  % Leave define mode.
netcdf.putVar(ncid, IDvarStation(1), 1, 44); 
netcdf.putVar(ncid, IDvarStation(2), 1, 38);
netcdf.putVar(ncid, IDvarRecord(1), 1,  L, DATA.Time - ...
  datenum('1970.01.01', 'yyyy.mm.dd'));
  ind= 1:numel(str); ind(strcmp(str,'Time'))= [];
  for k= ind
    netcdf.putVar(ncid, IDvarRecord(k), 1, L, a1D.(str{k}));
  end
  netcdf.close(ncid);
  %% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid)
  ncinfo(strFile)
elseif ((nargin <=3)||((nargin > 3)&&~isfield(SParam, 'noTXT')))&&(nargout < 2)
  if nargin>3&&isfield(SParam, 'bNotSaveTimeString')
    DATA.TimeDbl= DATA.Time; DATA= rmfield(DATA, 'Time');
    t= int32(numel(fieldnames(DATA))); ind= [t 1:(t-1)];
    DATA= orderfields(DATA, ind);
  end
  DATA= insertNaNs(DATA, separator); fprintf('.');
  saveTxt(rmfield(DATA,{'Inclination','Heading','Pitch','Roll'}), ...
    [SLoad.SaveDirName SLoad.strFileAdd]);
  saveTxt(rmfield(DATA,{'Vabs','Direction'}), ...
    [SLoad.SaveDirName SLoad.strFileAdd 'IHPR']);
end
end
