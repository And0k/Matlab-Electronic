function rasp(DataDirMask, txtFileCoef)
%%Example:
%% rasp('d:\WorkData\Cruises\_BlackSea\1110\RASP\source\1\111015_080130.TXT', ...
%% 'd:\WorkData\Cruises\_BlackSea\1110\RASP\source\1\coef111110.txt');
%
%DataDirMask= 'd:\WorkData\Experiment\Raspred\111110rasp,chein\R3\DATA*.TXT';
%T1	T2	T3	T4	T5
%0040	9DAE	A716	D2D0
SLoad= struct('strDir',DataDirMask, 'col',struct(...
 'FirstRow',1, 'inDate','', 'DateForm',{'yyyy mm dd HH MM SS'},  ...
 'header',       ['P	T01	T02	T03' char(13) 'A A A A'], ...
 'headerReplace', 'P	T01	T02	T03	Time', ...
 'ColSeparator',char(9), 'NotWorningNewFields',{'Time'}));
SLoad.col.DateForm= {[] [] 5}; %SLoad.col.DateForm{3} is sampling period, s
[SLoad SLog a1D]= loadData2DB(SLoad);
for k=find(~strcmp(SLoad.col.Name,'Time'))
  str= char(a1D.(SLoad.col.Name{k})); str(str==char(13))= ' ';
  a1D.(SLoad.col.Name{k})= hex2dec(str);
end
a1D.Time= a1D.Time + datenum(SLoad.FilesInfo.name( ...
  (1:(find(SLoad.FilesInfo.name=='.',1)-1))), 'yymmdd_HHMMSS') + ...
  + datenum('01', 'dd') + 10*dayHz;
a1D.P= -a1D.P; %datestr(a1D.Time([1,end]))
% figure('Name', 'All pressure'); plot(a1D.P)

SCoef= struct('strDir',txtFileCoef, 'FirstRow',1, 'col',struct( ...
 'header', '*Sensor	p%d', 'ColSeparator',char(9)));
[~,~,sc]= loadData2DB(SCoef);

str= setdiff(fieldnames(a1D), {'Time'});
sc.Name= str;
% nCh= numel(str);
% clr= hsv(nCh);
[DATA I]= forApplyCoefficients(a1D, ...
  fliplr(cell2mat(struct2cell(rmfield(sc, {'Sensor','Name'}))')), sc.Sensor', sc);
DATA= orderfields(DATA, [numel(str)+1 1:(numel(str))]);
%% saveCTD(DATA, struct('Ens',1, 'EnsEnd',numel(DATA.iTimeDbl)), [], ...
%   SLoad.strDir, {'xyz'}); DATA.P
saveTxt(DATA, [SLoad.SaveDirName 'R']);




