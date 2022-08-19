SLoad= struct('strDir','', 'col',struct( 'FirstRow',1, ...
  'inDate','', 'DateForm',{'yyyy mm dd HH MM SS'}, ...
  'header'       , '*yyyy mm dd HH MM SS	P0 T00 0 Power P1 P2 P3 T%02d', ...
  'headerReplace', ['yyyy mm dd HH MM SS         Time P0 P1 P2 P3 Power T00 T%02d', ...
  char(13), 'Date1 Date2 Date3 Date4 Date5 Date6 Time P0_dBar P1_dBar P2_dBar P3_dBar Power_V T00 T%02d'], ...
  'ColSeparator',char(9), 'notWarningLostFields',{'Time'}, 'valueFormat','u'));
FilterDat= struct();

%%
if true % false % 
  SLoad.strDir= 'd:\WorkData\Cruises\_BalticSea\140719_Nord3\termochain\140719TCh01.txt';
  txtFileCoef= 'd:\WorkData\Cruises\_BalticSea\140719_Nord3\termochain\coef#N1_140403cur_forMatLab.txt';
  DB= struct('name','d:\WorkData\Cruises\_BlackSea\1110\BlS1110.mdb', ...
    'tbl','Tchain processed', 'tblProcessed', struct('log', 'Tchain sections log'));
  P= struct('LT',repmat(2, 1,9)); P.LT(8)= P.LT(8)+2; P.LT(9)= 1.5;
  %DB.name= [SLoad.strDir(1:(end-3)) 'h5'];
  P= struct('LT',2); %, 'PPoints',[1,11,21,31]
  %FilterDat.Power_V.StartFiltBelow= 6.3909;
elseif false
  SLoad.strDir= 'd:\WorkData\Cruises\_BlackSea\1309\_source\Termochain\N3\test\tank\1b.bin';
  txtFileCoef= 'd:\WorkData\Cruises\_BlackSea\1309\_source\Termochain\N3\test\tank\coef130925#TC3.TXT';
  P= struct('LT',2, 'PPoints',[1,11,20,29]);
elseif true
  %TC2
  SLoad.col.header= ...
  ['*Time mm dd HH MM SS P0 T00 0 Power P1 P2 T%02d', char(13), ... 
  '"25/" 01	01 01	00 32	04288	31937	22566	20323	04093	01792	39223	39629	39594	39576	39447	39738	39810	39633	39677	39499	40013	39680	39664	39805	39852	39612	39975'];
                        %P0   T00   0     Power P1    P2    T21   T24   T25   T26   T28   T29   T30   T31   T32   T33   T34   T35   T36   T37   T38   T39   T40
  %P3 remove
  SLoad.col.headerReplace= regexprep(SLoad.col.headerReplace, ...
  {' P3\w*' 'yyyy' 'mm' 'dd' 'HH' 'MM' 'SS' 'Date(?<D>[1-9])'},  '');
  SLoad.col.DateForm= 'xdd/mm/yyyy HH:MM:SSxxxxxx';
  SLoad.col.ColSeparator= char(9);
  %SLoad.col.Delimiter= [char(9) ','];
  P= struct('PPoints',[1,8,18]);
  
  SLoad.strDir= 'd:\WorkData\Cruises\_BlackSea\1306\_source\termochain\#2_onCable\13*.txt';
  txtFileCoef= 'd:\WorkData\Cruises\_BlackSea\1306\_source\termochain\#2_onCable\coef#TC2_130325-T22,23,27.txt';
  DB= struct('name','d:\WorkData\Cruises\_BlackSea\1110\BlS1110.mdb', ...
    'tbl','Tchain processed', 'tblProcessed', struct('log', 'Tchain sections log'));  
  DB.name= [SLoad.strDir(1:(end-3)) 'h5'];
else
  %%
  SLoad.strDir= 'd:\WorkData\Cruises\_BlackSea\1309\_source\Termochain\N1\Kosa_Depth=28m\#Tch1.txt';
  txtFileCoef= 'd:\WorkData\Cruises\_BlackSea\1309\_source\Termochain\N1\Kosa_Depth=28m\coef#TC1_130323-1.txt';
  DB= struct('name','d:\WorkData\Cruises\_BlackSea\1110\BlS1110.mdb', ...
    'tbl','Tchain processed', 'tblProcessed', struct('log', 'Tchain sections log'));
  P= struct('LT',0.95, 'PPoints',[1,10,20]);
  FilterDat.Time.min= 735508.6954;
  %P3 remove
  strRem= {'P1' 'P2' 'P3'};
  k= arrayfun(@(x) strfind(SLoad.col.header, x{:}), strRem);
  ind= []; for t= 1:numel(strRem); ind= [ind k(t):(k(t)+size(strRem{t}))]; end;
  SLoad.col.header(ind)= []; 
  k= arrayfun(@(x) strfind(SLoad.col.headerReplace, x{:}), strRem, 'UniformOutput',false);
  ind= [];
  for t= 1:numel(strRem)
    for tt= 1:numel(k{t})
    ik= k{t}(tt)+length(strRem{t}) - 1;
      ind= [ind k{t}(tt):(ik+strfind(SLoad.col.headerReplace((ik+1):end), ' '))];
    end
  end
  SLoad.col.headerReplace(ind)= [];
  
end
termochain(SLoad, txtFileCoef, P, [], FilterDat);