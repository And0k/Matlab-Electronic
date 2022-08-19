function DATA= termochain(SLoad, txtFileCoef, P, DB, FilterDat)
%"termochain(DataDirMask, txtFileCoef)" loads data from DataDirMask, 
%upply coefficients from txtFileCoef and save to DataDirMask
% tab-delimited file with columns: 
%"Time     	P     	TP     	Power     	T01     	T02     ...		Date     Time "
%Example:
% DataDirMask= 'd:\WorkData\Cruises\_BlackSea\1110\Split\source\SPLIT000.txt';
% txtFileCoef= 'd:\WorkData\Cruises\_BlackSea\1110\Split\source\coef111110.txt';
% termochain(DataDirMask, txtFileCoef)
%
%"termochain(DataDirMask, txtFileCoef, P)" saves in tab-delimited file
% with columns: "Time     	Temp     	Depth"
%which contains depth of each sensor.
%P should be set as array: [LT LP Bot] (see picture):
%Bot= depth of bottom (load): at BlackSea2011 it is calcululated as 
%BotX-min(|Press|) = 55-3.5
%LT is array of distanses between T sensors in order of its data columns:
%from T0 to T1, from T1 to T2, T2 to T3, ..., T(end-1) to Tend, Tend to Bot.
%, can be NaN if 
%If distanses are equal then you can
%specify LT as L1= distance between Press. and first Temp. probe or between
% any adjacent Temp. probes and rope length from Press. to bottom only
%
%             __Press,T0
%        L(1)/|
%       L(i)/ |_P(i)
%          ...
%         /   |
%        /    |
% L0toB_/     |_Bot
%so
% (P(i)-Press)/L(i) = (Bot-Press)/L0toB => P(i)= L(i)*(Bot-Press)/L0toB + Press
%where i - temperature sensor number,
%L(i)= L1*i for equal length
%L0toB = sum(L(i)) + (Tend to Bot)
% Or set P as structure with fields:
%- Points - indexes of Tsensor with P
%- LT - array of distanses between T sensors in order of its data columns:
%from T0 to T1, from T1 to T2, T2 to T3, ..., T(end-1) to Tend, Tend to Bot
%
%Example:
% DataDirMask= 'd:\WorkData\Cruises\_BlackSea\1110\Split\source\SPLIT000.txt';
% txtFileCoef= 'd:\WorkData\Cruises\_BlackSea\1110\Split\source\coef111110.txt';
%
% termochain(DataDirMask, txtFileCoef, [2 51 55-3.5]);

%If provided 4-th argument saves to DataBase instead of tab-delimited file
%4-th argument is struct: 
%DB.name - system name of database to connect
%DB.tbl  - table name
%Example:
% DB= struct('name','BlS1110', 'tbl','Tchain processed');
% %old: termochain(DataDirMask, txtFileCoef, [2 51 55-3.5], DB);
% DataDirMask= 'd:\WorkData\Cruises\_BlackSea\1206\source\termochain\0715_0924.txt';
% txtFileCoef= 'd:\WorkData\Cruises\_BlackSea\1206\source\termochain\coef.txt';
% termochain(DataDirMask, txtFileCoef, [1 1 repmat(2, 1, 6) 4 -2 4 repmat(2, 1, 6) 1 41.1], DB);
if nargin<5; FilterDat= struct(); end
bInd= ischar(SLoad);
if bInd||(isstruct(SLoad)&&isfield(SLoad,'strDir'))
  if bInd
    SLoad= struct('strDir',SLoad, 'col',struct( 'FirstRow',1, ...
      'inDate','', 'DateForm',{'yyyy mm dd HH MM SS'}, ...
      'header'       , '*yyyy mm dd HH MM SS	P T00 TP          Power T%02d', ...
      'headerReplace', ['yyyy mm dd HH MM SS          Time P TP Power T00 T%02d', ...
      char(13), 'Date1 Date2 Date3 Date4 Date5 Date6 Time P0 TP Power T00 T%02d'], ...
      'ColSeparator',char(9), 'notWarningLostFields',{'Time'}));
  end
  [SLoad SLog a1D]= loadData2DB(SLoad, FilterDat);
else
  a1D= SLoad;
  SLoad= struct('TimeRange',a1D.Time([1 end]), 'strDir', ...
    'C:\Temp\termochain.txt', 'FilesInfo',[], 'SaveDirName','C:\', ...
    'strFileAdd','T-ch');
end
str= fieldnames(a1D);
ind= 1:numel(str);
bInd= strcmp(str, 'Time');  k= find(bInd);
bInd= strncmp(str, 'Power', numel('Power')); L= find(bInd);
ind([k L])= [];
bInd= strncmp(str(ind),'T',1);
a1D= orderfields(a1D, [k ind(bInd) ind(~bInd) L]);
str= fieldnames(a1D);
if ~issorted(a1D.Time) %datestr(a1D.Time([1 end]));
  warning('Check:LoadedTime', 'Time not sorted, sorting...')
  [a1D.Time, ind]= sort(a1D.Time); ind= int32(ind);
  for k=2:numel(str)
    a1D.(str{k})= a1D.(str{k})(ind);
  end  
  bInd= diff(a1D.Time)==0;
  if any(bInd)
    ind= [];
    if isfield(a1D,'P')
        for k=1:numel(str)
            if any(a1D.P(bInd)~=a1D.P([false; bInd]))
                ind= find(a1D.P(bInd)~=a1D.P([false; bInd]));
                disp(ind)
            end
        end
    end
    if any(ind)
      stopHere('duplicate Time (%d), but not equal rows... del?', sum(bInd))
    end
    for k=1:numel(str)
      a1D.(str{k})(bInd)= [];
    end
  end
end
%filter bad
for k=1:numel(str)
  a1D.(str{k})(a1D.(str{k})==0)= NaN;
end
SCoef= struct('strDir',txtFileCoef, 'col',struct(...
 'FirstRow',1, 'header','*Sensor	p%d', 'ColSeparator',char(9)));
[~,~,sc]= loadData2DB(SCoef);
sc.Name= setdiff(str, {'Time', 'TP', 'Power'}); %retain existed coef. only

% nCh= numel(str);
% clr= hsv(nCh);
[DATA I]= forApplyCoefficients(a1D, ...
  cell2mat(struct2cell(rmfield(sc, {'Sensor','Name'}))'), sc.Sensor', sc);
%DATA.T11(DATA.Time>735447.8993287)= NaN;
%DATA.Time= DATA.Time; DATA= rmfield(DATA,'Time'); a1D.P fliplr()
% ind= (find(strncmp(str,'P',1)&~strcmp(str,'Power')))';
% for k=ind; DATA.(str{k})= -DATA.(str{k}); end
%% saveCTD(DATA, struct('Ens',1, 'EnsEnd',numel(DATA.Time)), [], ...
%   SLoad.strDir, {'xyz'});
L= int32(size(DATA.Time,1));
bInd= cellfun(@(x) ~isempty(x)&&x==1, regexp(str, '^P\d', 'once'))';
strP= str(bInd);
%% Filter bad P where bad Power
if isfield(DATA,'Power_V')&&isfield(FilterDat,'Power_V')&& ...
  isfield(FilterDat.Power_V, 'StartFiltBelow') %= 6.3909; %Umin=14032.3015*0.0004175 + 0.5324
  iSt= find(DATA.Power_V < FilterDat.Power_V.StartFiltBelow, 1,'first');
  if ~isempty(iSt)
    fprintf(1, '\nFilter bad P where bad Power:\n');
    figure('Name', 'Filter bad P where bad Power'); hold on; ha= gca;
    for k= 1:numel(strP)
      plot(ha, 1:L, DATA.(strP{k}), 'Color', [.5 .5 .5]);
      plot(ha, iSt, DATA.(strP{k})(iSt), '+k');
    end
    plot(ha, 1:L, DATA.Power_V, 'Color', [1 0.5 0]);
    plot(ha, iSt, DATA.Power_V(iSt), '+k');
    
    ind= iSt:L;
    %Paramerers of top:
    nT= 1;
    m= nanmean(DATA.(strP{nT})(1:iSt));
    t= DATA.(strP{nT})(iSt:end) - m; t= medfiltSmooth(t, 9);
    dm= m-nanmin(DATA.(strP{nT})(1:iSt));
    %Correct others pressure using top
    for k= 2:numel(strP)
      mc= nanmean(DATA.(strP{k})(1:iSt));
      bInd= DATA.(strP{k})(ind)-mc < t; %bottom go down more than top!
      if any(bInd)
        dmc= mc - nanmin(DATA.(strP{k})(1:iSt));
        dmc= dmc/dm;
        if dmc<1; fprintf(1, '%s:%dpoints, k= %g\n', strP{k}, sum(bInd), dmc);
        else stopHere('May be incorrect coeff (%g) of pressure decrease', dmc);
        end
        DATA.(strP{k})(ind(bInd))= mc + t(bInd)*dmc; %correct
        plot(ha, ind(bInd), DATA.(strP{k})(ind(bInd)), 'Color', 'r');
      end
    end
    print(get(ha,'Parent'),'-dpng',[SLoad.SaveDirName SLoad.FilesInfo(1).name(1:end-4) 'correctP.png']);
    close(get(ha,'Parent'));
  end
end
bInd= cellfun(@(x) ~isempty(x)&&x==1, regexp(str, '^P\d', 'once'));
bInd= bInd|strncmp(str, 'Time', numel('Time'))|strncmp(str, 'Power', numel('Power')); %{'Time' 'TP' 'P0' 'P1' 'P2' 'P3' 'Power'}
nT= int32(sum(~bInd));
ind= nT*L;
if nargin>2
  a1D= struct('Temp',NaN(ind,1), 'Depth',NaN(ind,1), 'Time',NaN(ind,1), ...
    'FilterMask',zeros(ind,1,'uint8'));
  a1D.Time(:)= repmat(DATA.Time', nT, 1);
  a1D.Temp(:)= cell2mat(struct2cell(rmfield(DATA, str(bInd)))')';
  m= numel(strP);
  if isstruct(P)
    if isfield(P,'LT')
      if numel(P.LT)==1
        dists= P.LT*(1:double(nT));
%         LT= P.LT;
        %       else
        %         stopHere('Check it');
      else
        if numel(P.LT)>(nT-1)
          dists= cumsum([0 P.LT(1:(nT-1))]);
          L0toB= dists(end)+P.LT(end);
        else
          dists= cumsum([0 P.LT]);
          L0toB= dists(end)+1;
          stopHere('Length of full chain to load is not defined - added 1m to load');
        end
%         LT= dists(end);
      end     
    end
    if isfield(P,'Bot'); Bot= P.Bot;
    else
      Bot= min(DATA.(strP{end})(int32(L*0.2):int32(L*0.8)))+L0toB;      
      stopHere('Length of depth is not defined - think it = %g', Bot);
    end
    if isfield(P,'PPoints'); P= int32(P.PPoints);
    elseif m==1; P= int32(1);
    end
  else
    Bot=-P(end);
    if n==3
      dists= cumsum([0 P(1)*(1:(double(nT)-1))]);
      L0toB= P(2)+dists(end);
    else
      dists= cumsum([0 P(1:(nT-1))]);
      L0toB= P(end)+dists(end); %+sum(P(end-1)); ???
    end
  end
  distsDiv= dists/L0toB; %=(Bot-Press)/L0toB
  n= numel(P); 
  if n>1
    if n==m
      ind= 0:nT:(nT*L-1);
      for k= 1:n %where P with T
        DATA.(strP{k})= rep2mean(DATA.(strP{k}));
        iSt= P(k);
        a1D.Depth(iSt+ind)= DATA.(strP{k});
        a1D.FilterMask(iSt+ind)= uint8(iSt);
      end
      if ~any(P)~=nT %exist edge T sensor without P
        %calc P
        P0toB= Bot - DATA.(strP{end});
        bInd= -P0toB > L0toB;
        ind= nT+ind;
        a1D.Depth(ind( bInd))= DATA.(strP{end})( bInd) - dists(end); %Probe is not set yet
        a1D.Depth(ind(~bInd))= DATA.(strP{end})(~bInd) + distsDiv(end)*P0toB(~bInd); %linspace(DATA.(strP{1})(m), -P, nT);
        ind= 1:nT; ind([P, end])= [];
      else
        ind= 1:nT; ind(P)= [];
      end
      iSt= int32(0);
      for m= 1:L
        a1D.Depth(iSt+ind)= interp1(P, a1D.Depth(P+iSt), ind);
        a1D.FilterMask(iSt+ind)= uint8(ind);
        iSt= iSt + nT;
      end
    else
      stopHere();
    end      
  else
    ind= 1:nT;
    P0toB= Bot - DATA.(strP{1});
    bInd= -P0toB > L0toB;
    for m= 1:L
      if bInd(m)
        a1D.Depth(ind)= DATA.(strP{1})(m) - dists; %Probe is not set yet
      else
        a1D.Depth(ind)= DATA.(strP{1})(m) + distsDiv*P0toB(m); %linspace(DATA.(strP{1})(m), -P, nT);
      end
      ind= ind + nT;
    end
  end
end
if nargin>3
%   conn= database(DB.name, '', '');
%   save2DB(conn, DB.tbl, a1D);  
%   %str= fieldnames(a1D)'; str= [str {'Update'}]; %(~strcmp(str,'Time'))
%   %save2DB(conn, DB.tbl, a1D, {'Time' 'Depth'}, str);      
%   close(conn);
  DB.tbl= struct('D1', DB.tbl, 'log','Tchain sections log');
  
  Track.TimeRange= SLoad.TimeRange';
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
  %% DB
  save2DBL(struct('DB',DB), DB.tbl, a1D, [], Track, fieldnames(a1D));
  %save2DBL(struct('DB',DB), DB.tbl, a1D, [], Track,[fieldnames(a1D);{'Update'}]); %very long!
  t= diff(Track.TimeRange);
  if t > 10
    save2DB_TimeSections(a1D.Time, 10, struct('DB',DB), true); %, a1D, 240
  end
  %%
  stopHere('Continue ("F5") to save netCDF?')
  if t > 30
    strFile= [datestr(a1D.Time(1),'yymmdd-') ...
      datestr(a1D.Time(end),'mmdd')];
  else
    strFile= [datestr(a1D.Time(1),'yymmdd_HHMM-') ...
      datestr(a1D.Time(end)  ,'dd_HHMM')];
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
  netcdf.putAtt(ncid,k, 'title','Black sea termochain');
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
else
  saveTxt(DATA, [fileparts(SLoad.SaveDirName(1:end-1)) '\' SLoad.strFileAdd]); %strFile= 
end

