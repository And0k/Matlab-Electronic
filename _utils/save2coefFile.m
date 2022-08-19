function txtFileCoef= save2coefFile(kT, SLoad, strHeader, strNames)
%Saves coefficients to text file named "coef[DATE].txt" by compose some
% SLoad felds
%   kT - coefficients matrix kT
%   SLoad - struct with fields:
%DataDirName - directory where to save. If omited retrieved from SLoad.strDir
%TimeRange - must contain serial date to determine [DATE] part of string
%FilesInfo - struct with field "name". if it has '#' char adds string 
%after "coef" with all chars after '#'. If omited retrieved from SLoad.strDir
%strFileAdd - string to add more before [DATE]
%   strHeader - optional first string replacement,if omitted generates 
%default names: Probe, k0, k1*x, k2*x^2 ...
%   strNames -  optional first column channel names,if omitted generates
%default names: Ch1, Ch2...
%   Returns full file name. If kT is empty not anything saves

%Examples
%Save to txt file:
%SLoad.strFileAdd= '_'; txtFileCoef= ...
%save2coefFile(kT, SLoad, 'Probe     	C	A1x	A2x	A3x', strcat('G', num2cell(strXYZ')));
%Save to d:\Work\MatlabHistory\mat\coef_Electronics.mat:
%save2coefFile(curInc, SLoad, 'Inclinometr');
strk= '';
if isstruct(SLoad)
  if isfield(SLoad,'FilesInfo')
    c= strfind(SLoad.FilesInfo.name, '#')+1;
    if any(c); strk= SLoad.FilesInfo.name(c(end):(end-4)); else strk= ''; end
  elseif isfield(SLoad,'strDir')
    c= strfind(SLoad.strDir, '#')+1;
    if any(c); strk= SLoad.strDir(c(end):(end-4)); else strk= ''; end
  end
  if isfield(SLoad,'strDir')&&~isfield(SLoad,'DataDirName')
    SLoad.DataDirName= [fileparts(SLoad.strDir) '\'];
  end
else
  %istrk= SLoad;
end
[r c]= size(kT); c= c-1;
if isstruct(kT)
  %% Add current coeficients to file with all coeficients %%%%%%%%%%%%%%%%%%
  if exist('d:\Work\MatlabHistory\mat\coef_Electronics.mat', 'file')
    S= load('d:\Work\MatlabHistory\mat\coef_Electronics.mat', strHeader);
  else
    S.(strHeader)(1:15)= struct('i', int32(0), ...
      'TimeRange',[NaN NaN], 'TimeProcessed', NaN, ...
      'H', struct('A',NaN(3,3),'C',NaN(3,1)), ...
      'G', struct('A',NaN(3,3),'C',NaN(3,1)), 'Vabs',NaN);
  end
  k= str2double(strk);
  if isnan(k)
    if isfield(kT,'i'); k= kT.i;
    else stopHere('Item not specified!');
    end
  end
  bInd= [S.(strHeader).i]==k;
  k= int32(k);
  p= sum(bInd);
  if p==1
    p= find(bInd);
    fprintf('Will update item %d by i=#%d: ', p, k);
  elseif p>1;
    p= find(bInd);
    stopHere('write prog');
    bInd= ~((SLoad.TimeRange(1)< [S.(strHeader)(p).TimeRange(1)])| ...
            (SLoad.TimeRange(1)>=[S.(strHeader)(p).TimeRange(2)])); %#ok<NBRAK>
    if sum(bInd)==1
      SCoef= S.(strHeader)(p(bInd));
    end
  else
    if S.(strHeader)(k).i~=0
      stopHere();
    else
      fprintf('Write to empty "%s" item %d: ', strHeader, k);      
    end
  end
  strk= fieldnames(kT);
  for c=1:numel(strk)
    if strncmp(strk{c}, 'Time', 4)
      fprintf('%s= ', strk{c})     
      for t= 1:numel(kT.(strk{c}))
        if isnan(kT.(strk{c})(t)); fprintf('NaN ')
        else
          fprintf('%s ', datestr(kT.(strk{c})(t), 'dd.mm.yyyy HH:MM:SS'));
        end
      end
      fprintf(', ');
    elseif strcmp(strk{c}, 'i'); fprintf('i=%d ', kT.i);
    else fprintf('%s, ',strk{c})
    end
    S.(strHeader)(k).(strk{c})= kT.(strk{c});
  end
  if S.(strHeader)(k).i==0; S.(strHeader)(k).i= k; fprintf('i=%d ',k); end
  if ~isfield(S.(strHeader),'TimeRange')||~isfinite(S.(strHeader)(k).TimeRange(1))
    fprintf('TimeRange(1)= %s', datestr(SLoad.TimeRange(1), ...
              'dd.mm.yyyy HH:MM:SS'))
    S.(strHeader)(k).TimeRange(1)= SLoad.TimeRange(1);
  end
  if ~any(strcmp(strk,'TimeProcessed'))
    S.(strHeader)(k).TimeProcessed= now();
  end
  s= input('.\nPress "Enter" to proceed, Ok? ');
  if isempty(s)
    save('d:\Work\MatlabHistory\mat\coef_Electronics.mat', '-struct', 'S', ...
      strHeader, '-append');
  else
    fprintf('\ndone nothing');
    return;
  end
  fprintf('- Ok');
else
  if ~isfield(SLoad,'strFileAdd'); SLoad.strFileAdd= ''; end
  if ~isempty(strk); strk= ['#' strk]; end
  txtFileCoef= [SLoad.DataDirName 'coef' strk SLoad.strFileAdd datestr( ...
  SLoad.TimeRange(1),'yymmdd') '.txt'];
  if r<1; return; end
  fid= fopen(txtFileCoef, 'wt');
  if(nargin>2)&&ischar(strHeader); fprintf(fid, strHeader);
  else fprintf(fid, ['Probe     \tk0     \tk1*x' repmat('     \tk%d*x^%d', ...
      1, c-1)], [2:c; 2:c]);
  end
  if(nargin>3)
    strk= ['\n%s\t%8.8f' repmat('\t%8.8g', 1, c)];
    for p=1:r; fprintf(fid, strk, strNames{p}, kT(p,:)); end
  else
    strk= ['\nCh%d' repmat('\t%8.8g', 1, c+1)];
    for p=1:r; fprintf(fid, strk, p, kT(p,:)); end
  end
  fclose(fid);
end
end

