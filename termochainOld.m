function DATA= termochain(DataDirMask, txtFileCoef, P, DB)
%"termochain(DataDirMask, txtFileCoef)" loads data from DataDirMask, 
%upply coefficients from txtFileCoef and save to DataDirMask
% tab-delimited file with columns: 
%"Time     	P     	TP     	Power     	T01     	T02     ...		Date     Time "
%Example:
% DataDirMask= 'd:\WorkData\Cruises\_BlackSea\1110\Split\source\SPLIT000.txt';
% txtFileCoef= 'd:\WorkData\Cruises\_BlackSea\1110\Split\source\coef111110.txt';
% termochain(DataDirMask, txtFileCoef)
%
%"termochain(DataDirMask, txtFileCoef, P)" saves in tab-delimited file with columns:
%"Time     	Temp     	Depth"
%which contains depth of each sensor.
%P should be set as array: [d1 dB Bot] (see picture):
%d1= distance between Press. and first Temp. probe or between any adjacent
%Temp. probes
%Bot= depth of bottom (load)
%dB= distance between Press. and bottom (=Bot - min(|Press|) =55-3.5 in BlackSea2011)
%
%        __Press
%     d1/|
%  i*d1/ |_P(i)
%     ...
%    /   |
% dB/    |_Bot
%
%so
% d1*i/(P(i)-Press) = dB/(Bot-Press) => P(i)= d1*i*(Bot-Press)/dB + Press
%where i - temperature sensor number
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
% termochain(DataDirMask, txtFileCoef, [2 51 55-3.5], DB);


SLoad= struct('strDir',DataDirMask, 'col',struct('FirstRow',0, ...
 'inDate','', 'DateForm',{'yyyy mm dd HH MM SS'}, ...
 'header'       , '*yyyy mm dd HH MM SS               P TP Power T%02d' , ...
 'headerReplace', ['yyyy mm dd HH MM SS          Time P TP Power T%02d' , ...
  char(13), 'Date1 Date2 Date3 Date4 Date5 Date6 Time P TP Power T%02d'], ...
 'ColSeparator',char(9), 'NotWorningNewFields',{'Time'}));
[SLoad SLog a1D]= loadData2DB(SLoad); a1D.P= -a1D.P;
%datestr(a1D.iTimeDbl([1 end]));

SCoef= struct(...
 'strDir',txtFileCoef, 'FirstRow',1, 'col',struct( ...
 'header', '*Sensor	p%d', 'ColSeparator',char(9)));
[~,~,sc]= loadData2DB(SCoef);

str= fieldnames(a1D);
sc.Name= setdiff(str, {'iTimeDbl', 'TP'});

% nCh= numel(str);
% clr= hsv(nCh);
[DATA I]= forApplyCoefficients(a1D, ...
  fliplr(cell2mat(struct2cell(rmfield(sc, {'Sensor','Name'}))')), sc.Sensor', sc);
%DATA.Time= DATA.iTimeDbl; DATA= rmfield(DATA,'iTimeDbl'); a1D.P
k= find(strcmp(str, 'Time'));
DATA= orderfields(DATA, [k 1:(k-1) (k+1):numel(str)]); DATA.P= -DATA.P;
%% saveCTD(DATA, struct('Ens',1, 'EnsEnd',numel(DATA.iTimeDbl)), [], ...
%   SLoad.strDir, {'xyz'});
nT= int32(numel(setdiff(sc.Name, {'Time', 'TP', 'P', 'Power'})));
if nargin>2
  d1= P(1);
  dB= P(2);
  Bot=-P(3);
  a1D= struct('Temp',   NaN(nT*size(DATA.Time,1),1), ...
              'Depth',  NaN(nT*size(DATA.Time,1),1), ...
              'iTimeDbl',   NaN(nT*size(DATA.Time,1),1));
  ind= 1:nT; dists= d1*double(ind)/dB;
  for m= 1:size(DATA.Time,1)
    a1D.iTimeDbl(ind)= DATA.Time(m);
    a1D.Temp(ind)= structfun(@(x) x(m), rmfield(DATA, {'Time', 'TP', 'P', ...
     'Power'}));
    if -(Bot-DATA.P(m))>dB
      a1D.Depth(ind)= DATA.P(m) - d1*(1:double(nT)); %Probe is not set yet
    else
      a1D.Depth(ind)= DATA.P(m) + dists*(Bot-DATA.P(m)); %linspace(DATA.P(m), -P, nT);
    end
    ind= ind + nT;
  end
end
if nargin>3
  conn= database(DB.name, '', '');
  save2DB(conn, DB.tbl, a1D);  
  %str= fieldnames(a1D)'; str= [str {'Update'}]; %(~strcmp(str,'iTimeDbl'))
  %save2DB(conn, DB.tbl, a1D, {'iTimeDbl' 'Depth'}, str);      
  close(conn);
else
  saveTxt(DATA, SLoad.SaveDirName);
end

