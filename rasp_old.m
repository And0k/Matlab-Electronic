function rasp(DataDirMask, MFileCoefficients)
%mmddHHMM
[SLoad N_Files]= loadPrepare(struct('strDir', DataDirMask, ...
  'MFileCoefficients', MFileCoefficients, 'strFileAdd','p'));
for i_file=1:numel(N_Files)
    SLoad.strDir= [SLoad.DataDirName SLoad.FilesInfo(i_file).name];
    [a1D SLoad]= loadData(SLoad);
    
    a1D= forApplyCoefficients(a1D, k, Ch, SLoad);
    k= struct2cell(a1D); k= k(cellfun(@isfloat, k));
    DATA= reshape([k{:}],[size(a1D.Time,1), numel(k)]);
    [DATA, SLoad]= calcAddData(DATA, CalcData, SLoad);
    iTime= find(strcmp(SLoad.col.Name, 'Time'));
    DATA(:, iTime)= (DATA(:, iTime) - DATA(1, iTime))/(dayHz*3600);
    fprintf(1,'\nSave runs '); 
    StrSave= cell2mat(strcat(SLoad.col.Name,'\t'));  StrSave(end)='n';
    SaveMask= [repmat('%8.7g\t',1,size(SLoad.col.Name,2)), '\n'];
    StrFileShortName= [SLoad.StrParVal{i_file,1}(1:end-4), 'p'];
    StrFile= strcat(SLoad.DataDirName, StrFileShortName, '.txt');
    BaklanStat(DATA, SLoad.col.Name, StrFile);
    fid= fopen(StrFile, 'wt');
    fprintf(fid, StrSave);
    fprintf(fid, SaveMask, DATA');
    fclose(fid);  
end
   
% saveCTD(a1D,[],[], 'D:\WorkData\_Cruises\_BlackSea\0809\Raspred\source\10091611R1__\', {'Runs'});
% fid= fopen('D:\WorkData\_Cruises\_BlackSea\0809\Raspred\Calibration\Raspred1Contro2.txt', 'wt');
% for k=1:numel(Pout)
%     fprintf(fid,'%g\t%g\t%g\t%g\t%s\n', Pout(k), T2out(k), T3out(k), T4out(k), datestr(Tdat(k), 'yyyy.mm.dd HH:MM:SS'));
% end
% fclose(fid);
