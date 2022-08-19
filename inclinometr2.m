function [DATAout a1D ax SLoad]= inclinometr2(DataDirMask, FilterDat, SParam)
%"inclinometr(DataDirMask)" loads data from DataDirMask, upplies 
%coefficients from SParam.coefFile depended on file name and save to
%DataDir
% tab-delimited file with columns: 
%"Time    	Heading     	Peach     	Roll     Date     Time "
%nargin=1 => use coefficients [0 1]
%SParam:
% prefix - word before probe number
% bNotShowFigure= false => Not to display figure,
%   [] => Not to save time string column (less memory demand)
%   struct with fields Time and raw Gx Gy Gz Hx Hy Hz) => Not to load data
%'noTXT' => not saves text files
%'DB' => saves to DataBase instead of tab-delimited file
%argument DB must be struct: 
%DB.name - system name of database to connect
%DB.tbl  - table name

% Examples
% SParam.coefFile= 'd:\WorkData\Cruises\_BalticSea\130510\_source\inclinometr\coef#7.txt';
% inclinometr( ...
%   'd:\WorkData\Cruises\_BalticSea\130510\_source\inclinometr\#7.txt', ...
%   struct('Time',struct('NoDataSeparator', 1.1*dayHz), ...
%   'smooth',[3 33], 'decimate',5), SParam);

%DataDirMask= 'd:\WorkData\Experiment\naclinometr\COMPAS\130430.txt';
%SParam.coefFile= 'd:\WorkData\Experiment\naclinometr\coef130506.txt';
%inclinometr(DataDirMask, true, SParam);
% Not to display figure: inclinometr(DataDirMask, false, SParam);
% Not save time string column: DATA= inclinometr(DataDirMask);
% Not to load data: DATA= inclinometr(DataDirToSave, [], a1D);
% Saving to DB:
%SParam.DB= struct('name','BlS1110', 'tbl','Tchain processed');
% %old: inclinometr(DataDirMask, true, SParam);
DATAout= [];
if nargin<2||isempty(FilterDat); FilterDat= struct(); end
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
  'headerReplace', ['yyyy mm dd HH MM SS P Power Gx Gy Gz Hx Hy Hz' , ...
  char(13), 'Time Date1 Date2 Date3 Date4 Date5 P Power Gx Gy Gz Hx Hy Hz'], ...
  'ColSeparator',char(9), 'NotWarningNewFields',{'Time'}, ...
  'valueFormat','u'), 'strCmd','D', ...
  'DB', struct('name', '.h5'));
if isfield(SParam, 'SaveDirName'); SLoad.SaveDirName= SParam.SaveDirName;
else SParam= rmfield(SParam, 'SaveDirName');
end
if ~isfield(SParam,'prefix'); SLoad.prefix= '#';
else SLoad.prefix= SParam.prefix; SParam= rmfield(SParam, 'prefix');
end
str= ''; a1D= struct('Time',[]);
if (nargin<3)||(~isfield(SParam, 'a1D'))
  if strcmp(DataDirMask((end-2):end), 'mat'); load(DataDirMask); b= false;
  elseif strfind(DataDirMask, '.h5/')
    n= strfind(DataDirMask, '.h5/');
    if SLoad.prefix=='#'; SLoad.prefix= DataDirMask(n+4: n+2+...
        regexp(DataDirMask((n+4):end), '\d+', 'start'));
    end
    SLoad.strDir= DataDirMask(1:n+2);
    SLoad= rmfield(SLoad, 'col');
    SLoad= loadPrepare(SLoad);
    SLoad.FilesInfo.name= DataDirMask(n+4:end);
    %SLoad.DataDirName= [fileparts(SLoad.strDir) '\']; SLoad.SaveDirName
    DataDirMask= SLoad.strDir; %to not gen.warning of moved file
    b= false;
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
    if (~b)&&~isempty(str); SLoad.strDir= DataDirMask;
    else
      k= strfind(DataDirMask,'incl_src.mat');
      if ~isempty(k); SLoad.strDir= [DataDirMask(1:(k-1)) '.txt']; end
    end
    %SLoad.FilesInfo.name= SLoad.strDir;
    SLoad= loadPrepare(rmfield(SLoad, {'col','SaveDirName','DataDirName','strFileAdd'}));
    fprintf(' in renamed dir ');
  else
    p= strfind(SLoad.FilesInfo.name, SLoad.prefix);
    if strcmp(DataDirMask(end-2:end), '.h5');      
    else
      fprintf('\nAccording to loaded mat-file settings result will be saved to\n%s\n', ...
        fileparts( SLoad.SaveDirName(1:end-1)));
      if any(p)&&strcmp(SLoad.strFileAdd,'_') %for old mat
        SLoad.strFileAdd= [SLoad.FilesInfo.name(p(end):(strfind(SLoad.FilesInfo.name)-4)) SLoad.strFileAdd];
      end
    end
  end
else
  a1D= SParam.a1D;
end
k= find(SLoad.FilesInfo.name=='#', 1, 'last');
if isempty(k)
  if ~strncmp(SLoad.FilesInfo.name, SLoad.prefix, numel(SLoad.prefix))
    if SLoad.FilesInfo.name(1)=='i'&&strcmp(SLoad.prefix,'Inclinometer')
      k= 1;
    else k=0;
    end
  else k= numel(SLoad.prefix);
  end
else k= k+1;
end
if k>0
  Nprobe= str2double(regexp(SLoad.FilesInfo.name(k:end), '\d+', 'match', ...
    'once'));
  if isnan(Nprobe); Nprobe= 0; end
  SLoad.strFileAdd= sprintf('%c%02d\\', SLoad.prefix(1), Nprobe);
else Nprobe= -1;
end
if (nargin<3)||~isfield(SParam,'coefFile')||isempty(SParam.coefFile)
  if isempty(a1D.Time)
    if strcmp(DataDirMask(end-2:end), '.h5')
      a1D= h5read(DataDirMask, sprintf('/%s/table', SLoad.FilesInfo.name),1,1);
      a1D.Time= datenum('01.01.1970','dd.mm.yyyy')+double(a1D.index*1e-8)./(24*36000);
    else
      stopHere('data exists?')
    end
  end
  SCoef= coefLoad(a1D.Time, SLoad.prefix, Nprobe);
else
  if ~any(SParam.coefFile=='\')
    if strcmp(SParam.coefFile, '*')&&(SLoad.strFileAdd(1)=='#')
      SCoef.strDir= [fileparts(DataDirMask) '\coef' SLoad.strFileAdd(1:end-1) '*.??t']; %mat or txt
    else
      SCoef.strDir= [fileparts(DataDirMask) '\' SParam.coefFile];
    end
    SCoef= loadPrepare(SCoef);
  else
    SCoef.strDir= SParam.coefFile;
  end
  if strcmp(SParam.coefFile((end-3):end), 'mat')
    SCoef= load(SParam.coefFile);
    if isfloat(SCoef.Vabs); str= 'Vabs'; else str= '';
    end
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
    end   
  end
  if ~isempty(str)
    f= fittype(sprintf('poly%d', numel(SCoef.Vabs)-find(SCoef.Vabs, 1, 'first')));
    P= num2cell(SCoef.Vabs);
    SCoef.Vabs= cfit(f, P{:});
  end
end
if(~isfield(SCoef, 'Vabs'))||(isfloat(SCoef.Vabs)&&any(isnan(SCoef.Vabs)))
  warning('inclinometr:noCoef', 'no Vabs coeficients. Set to [0 1]');
  SCoef.Vabs= @(x) x;
end
if     ~isfield(a1D  , 'P'); SCoef.P= NaN;
elseif ~isfield(SCoef, 'P')
  if strcmp(SLoad.prefix,'WaveGage')
    warning(sprintf('inclinometr(%s):noCoef', SLoad.prefix), ...
      'no P coeficients. Set to [0 1]');
    SCoef.P= @(x) x;
  else SCoef.P= NaN;
  end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(DataDirMask(end-2:end), '.h5')
  strWhere= '';
  if isfield(FilterDat.Time, 'max')
    if isempty(strWhere); strWhere= {'-where'}; end
    strWhere= [strWhere {['index<' datestr(FilterDat.Time.max(), ...
      '"yyyy.mm.dd HH:MM:SS"')]}];
  end
  if isfield(FilterDat.Time, 'min')
    if isempty(strWhere); strWhere= {'-where'}; end
    strWhere= [strWhere {['index>' datestr(FilterDat.Time.min(), ...
      '"yyyy.mm.dd HH:MM:SS"')]}];
  end
  if isfield(SParam, 'chunkDays')
    if isempty(SParam.chunkDays); SParam.chunkDays= 10000; end
    strWhere= [strWhere {'-chunkDays'}, {sprintf('%g',SParam.chunkDays)}];
  end
  if ~isempty(strWhere)
    str= python('d:\Work\_Python3\_fromMat\h5toh5.py', DataDirMask, ...
      sprintf('/%s', SLoad.FilesInfo.name), strWhere{:}, '-columns','index');
    iStEn= sscanf(str(2:(end-1)),'%d');
    iStEn= [iStEn(1:(end-1))+1, iStEn(2:end)];
    iStEn(end)= iStEn(end)+1;
  else
    N= h5info(DataDirMask, sprintf('/%s', SLoad.FilesInfo.name)); %/table
    N= int32(N.Datasets.Dataspace.Size);
    iStEn= [1 N];
%     str= python('d:\Work\_Python\fromMat\h5toh5.py', DataDirMask, ...
%     sprintf('/%s', SLoad.FilesInfo.name));
  end

else
  %N= int32(size(a1D.Time, 1));
  iStEn= save2DB_TimeSections(a1D.Time, [1 1]);%, Separator
end
LastSep= int32(size(iStEn,1));
if LastSep>1; fprintf(' break in %d time intervals:', LastSep)
elseif ~LastSep
  stopHere('no data');
  1;
else          fprintf('...');
end
% a1D.Gxyz= NaN(N,3,'double');
% a1D.Hxyz= NaN(N,3,'double');
%% Process on chunks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p= 1:LastSep
  fprintf(' %2d,', p);
  n= iStEn(p,2) - iStEn(p,1) + 1;
  if strcmp(DataDirMask(end-2:end), '.h5')
     a1D= h5read(DataDirMask, sprintf('/%s/table', SLoad.FilesInfo.name), ...
       iStEn(p,1),n); %1,2,N);
     a1D.Time= datenum('01.01.1970','dd.mm.yyyy')+double(a1D.index)/(24*36000*1e8); %datestr(a1D.Time([1,end]));
     a1D= rmfield(a1D, {'index','U'});
  end
  %% Filter
  a1D= FSFilter(a1D, FilterDat, 'Time',0,[],[], [],SLoad);  %rem bad in Time
  n= numel(a1D.Time);
  str= fieldnames(FilterDat); b= strncmp(str, 'H',1);
  if any(b);        a1D= FSFilter(a1D, FilterDat, str(b),0, 'Time'); end %set NaNs in H   
  str(strcmp(str,'Time')|b)= [];
  if ~isempty(str); a1D= FSFilter(a1D, FilterDat, str,  [], 'Time'); end %not interp Bad in other fields
  %if nargin<2; DATA= []; ax= []; return; end
  %a1D= FSFilter(a1D, FilterDat, str); 
%% correct int32 NaN values 
str= {'Gx', 'Gy', 'Gz', 'Hx', 'Hy', 'Hz'};
for k= 1:6
  if ~isfloat(a1D.(str{k}))
    bInd= a1D.(str{k})==int32(0);
    a1D.(str{k})= double(a1D.(str{k}));
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
else
  fprintf('correct bad H component: ');
  str= {'Hx', 'Hy', 'Hz'};
  for k= [3 1 2] %process Z at first
    bInd= isnan(a1D.(str{k}));
    if any(bInd)
      fprintf('%s ', str{k});
      bGood= [true true true]; bGood(k)= false; separator= find(bGood);
      a1D.(str{k})= double(a1D.(str{k})); temp= sum((SCoef.H.A(bGood, bGood)* ...
        ([double(a1D.(str{separator(1)})(bInd)), double(a1D.(str{separator(2)})(bInd))]' ...
        - repmat(SCoef.H.C(bGood),1, sum(bInd)))).^2);
      temp(temp>1)= 1; temp= sqrt(1 - temp)/SCoef.H.A(k,k);     
      a1D.(str{k})= rep2mean(a1D.(str{k}), bInd); %interpolate sign
      bGood= a1D.(str{k})(bInd)<SCoef.H.C(k);
      temp(bGood)= -temp(bGood);
      a1D.(str{k})(bInd)= SCoef.H.C(k) + temp;
      a1D.(str{k})= rep2mean(a1D.(str{k}));       %interpolate remained NaNs
    end   
  end
  fprintf('\n');
end

  DATA= struct('Time',[], 'Vabs',[], 'Direction',[], ...
    'Inclination',[], 'Heading',[], 'Pitch',[], 'Roll',[]);  
  % Apply accelerometer coef
  try
    Gxyz= (SCoef.G.A*([double(a1D.Gx), double(a1D.Gy), ...
                double(a1D.Gz)]' - repmat(SCoef.G.C, 1, n)))';
    %Gxyz(imag(Gxyz)~=0)= NaN;
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
  temp= sqrt(sum(Gxyz.^2,2));
  a1D.Gsum= temp; temp(isnan(temp))= 1;
  a1D= FSFilter(a1D, FilterDat, 'Gsum', inf, 'Time'); %set NaNs in G
  a1D.Gsum= isnan(a1D.Gsum); %binary indicator
%   for k=3:3:21
%     a1D.Gsum(smooth(a1D.Gsum,k)>(1/k))= true;
%   end
  a1D.Gsum((smooth(a1D.Gsum,3)>(2/3))|(smooth(a1D.Gsum,5)>(1/5))| ...
    smooth(a1D.Gsum,15)>(1/15))= true;
  % Apply compas coef
  Hxyz= (SCoef.H.A*([double(a1D.Hx), double(a1D.Hy), ...
              double(a1D.Hz)]' - repmat(SCoef.H.C ,1, n)))';
% DATA= forApplyCoefficients(rmfield(a1D,{'Gx','Gy','Gz','Hx','Hy','Hz'}), ...
%   cell2mat(struct2cell(rmfield(sc, 'Probe'))'), sc.Probe', sc);
%DATA.Time= DATA.iTimeDbl; DATA= rmfield(DATA,'iTimeDbl'); a1D.P fliplr()
%DATA.Direction  = double(wrapToPi(atan2(tan(be.elevation),tan(-be.bank)) - hdg.heading)*180/pi);
 


%% Filter
% for k=1:numel(str) %filter bad
%   a1D.(str{k})(a1D.(str{k})==0)= NaN;
% end
%[false(iStEn(p,1),1) bBad false(N-iStEn(p,2),1)] a1D.Gy(bBad)= NaN; a1D.Gz(bBad)= NaN; rep2mean(bBad)

separator= int32(findSeparator(a1D.Time, FilterDat.Time.NoDataSeparator)); 
for k=1:3
  Gxyz(a1D.Gsum,k)= NaN;
  %Gxyz(:,k)= rep2mean(Gxyz(:,k), bBad, a1D.Time);
  if isfield(FilterDat, 'smooth')
    Gxyz(:,k)= separateDo(@(x,ind) medfiltSmooth(x(ind), FilterDat.smooth), ...
      Gxyz(:,k), separator);
    Hxyz(:,k)= separateDo(@(x,ind) medfiltSmooth(x(ind), FilterDat.smooth), ...
      Hxyz(:,k), separator);
  end
end
a1D= rmfield(a1D, 'Gsum');
if isfield(FilterDat, 'decimate')
  ind= 1:FilterDat.decimate:n;
  Gxyz= Gxyz(ind,:);
  Hxyz= Hxyz(ind,:);
  DATA.Time= a1D.Time(iStEn(p,1):FilterDat.decimate:iStEn(p,2)); % DATA.Time= DATA.Time(ind);
  separator= separator/FilterDat.decimate;
else DATA.Time= a1D.Time;
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

DATA.Vabs= NaN(size(DATA.Inclination), 'double');
%bInd= DATA.Inclination<pi/2;
DATA.Vabs= SCoef.Vabs(1./sqrt(tan(pi/2 - DATA.Inclination))); %polyval();
DATA.Vabs(imag(DATA.Vabs)~=0)= NaN;
DATA.Inclination= DATA.Inclination*180/pi; %wrapTo360();
DATA.Heading    = DATA.Heading*180/pi;
DATA.Pitch      = DATA.Pitch  *180/pi;
DATA.Roll       =-DATA.Roll   *180/pi;
if ~((isfloat(SCoef.P)&&any(isnan(SCoef.P))))
  a1D.P= single(a1D.P); bInd= diff(a1D.P)~=0;
  bInd= [bInd; true]&[true; bInd]; bInd= bInd&(a1D.P==0);
  DATA.P= feval(SCoef.P, a1D.P);
  DATA.P= rep2mean(DATA.P, bInd); %intrep only single zeros %a1D.P= a1D.P==0);
end
  if nargin<3||~isfield(SParam, 'bNotShowFigure')
% be = calcBankElevation(Gxyz(:,2), Gxyz(:,3), Gxyz(:,1));
% hdg= calcHeading(DATA.Roll, DATA.Pitch, a1D.Hy,Hxyz(:,3),a1D.Hx, false);
% DATA.Direction= double(-DATA.Heading*180/pi);
% 
% DATA.Inclination= double(DATA.Roll*180/pi);
% DATA.Heading    = double(-DATA.Heading*180/pi);
% DATA.Pitch      = double(DATA.Pitch*180/pi);
% DATA.Roll       = double(-DATA.Roll*180/pi);
% DATA.Direction= double(wrapToPi(atan2(tan((-DATA.Roll)*pi/180),tan((DATA.Pitch)*pi/180)) - ...
%%  -DATA.Heading*pi/180)*180/pi);
    strID= [datestr(DATA.Time(1), 'yymmdd_HHMM') SLoad.strFileAdd(...
      (SLoad.strFileAdd~='\')&(SLoad.strFileAdd~='/'))];
    ax= get(0,'ScreenSize'); ax= ax.*[(ax(3)*0.11)/ax(1) 1 0.89 0.9];
    figure('Name', sprintf('%d %s: source and result', p, strID), ...
      'Position', ax);
    ax(1)= subplot(4,1,1); hold on
    plot(ax(1), DATA.Time,sqrt(sum(Hxyz.^2,2)),'m', DATA.Time,Hxyz(:,1), ...
    'r' , DATA.Time,Hxyz(:,2),'g',  DATA.Time,Hxyz(:,3),'b');
    legend({'|H|','Hx','Hy','Hz'}, 'Location','NorthEastOutside');  legend boxoff
    set(ax(1), 'yLim', [-1.2 1.2]); datetick(ax(1), 'x', 'HH:MM'); grid on;
    if ~((isfloat(SCoef.P)&&any(isnan(SCoef.P)))) %strcmp(SLoad.prefix,'WaveGage')
      plot(ax(1), DATA.Time, (DATA.P/20-1), 'c');
    end
    ax(2)= subplot(4,1,2); hold on; plot(ax(2), DATA.Time, temp, 'm'); 
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
if nargin>2&&isfield(SParam, 'DB')
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
str= fieldnames(a1D);
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
elseif (~(isfield(SParam, 'noTXT')&&(SParam.noTXT)))
  DATA= insertNaNs(DATA, separator+1); fprintf('.');
  fldSecondary= {'Inclination','Heading','Pitch','Roll'};
  strDir= SLoad.SaveDirName;
  SLoad.SaveDirName= [SLoad.SaveDirName SLoad.strFileAdd];
  if strcmp(SLoad.prefix,'WaveGage')
    saveTxt(rmfield(DATA,[{'Vabs','Direction'} fldSecondary]), SLoad);
    b= isfield(SParam, 'bSaveAsInclinometr')&&SParam.bSaveAsInclinometr;
    if b
      DATA= rmfield(DATA,{'P'});
      str= SLoad.strFileAdd;
      if (SLoad.strFileAdd(end)=='\')||(SLoad.strFileAdd(end)=='\')
        SLoad.strFileAdd= [SLoad.strFileAdd 'Incl_data\' ]; %(end-1) SLoad.strFileAdd(end)
      else SLoad.strFileAdd= [SLoad.strFileAdd '_Incl_data'];
      end
    end
  else
    b= true;
  end
  if b
    SLoad.SaveDirName= [strDir SLoad.strFileAdd];
    saveTxt(rmfield(DATA, fldSecondary), SLoad.SaveDirName);
    fldSecondary= {'Vabs','Direction'};
    if(isfield(DATA,'P')); fldSecondary= [fldSecondary {'P'}]; end
    SLoad.SaveDirName= [SLoad.SaveDirName 'IHPR'];
    saveTxt(rmfield(DATA, fldSecondary), SLoad);
    if strcmp(SLoad.prefix,'WaveGage'); SLoad.strFileAdd= str; %return to saved
    end
  end
  SLoad.SaveDirName= strDir;                                 %return to saved
end
if nargout
  if ~isstruct(DATAout); DATAout= DATA; strWhere= fieldnames(DATA);
  else %DATAout= addFieldsFromStruct(DATAout, DATA);
    for nL=1:numel(strWhere)
      DATAout.(strWhere{nL})= [DATAout.(strWhere{nL}); DATA.(strWhere{nL})];
    end
  end
end
end
%if isfield(DATAout,'')
SLoad.coef= SCoef;


end
