function coefDisp(strFile, probeType)
%% Display all fields in struct strFile or coef. saved in file 'strFile'

if (nargin>0)&&isstruct(strFile)
  ind= 1:numel(strFile);
  probeType= 'struct';
  S= struct(probeType, strFile);
else %try load file
  if nargin<1||isempty(strFile)
    strFile= 'd:\Work\MatlabHistory\mat\coef_Electronics.mat';
  elseif ~exist(strFile, 'File')
    stopHere('no file %s', strFile);
  end
  if nargin<2; probeType= 'Inclinometer'; end
  S= load(strFile, probeType);
  ind= find([S.(probeType).i]~=0);
end
str= fieldnames(S.(probeType));
for k= ind
  fprintf(1, '\n%s(%d). ', probeType, k);
  for p=1:numel(str)
    if isstruct(S.(probeType)(k).(str{p}))
      disp(str{p});
      strt= fieldnames(S.(probeType)(k).(str{p}));
      for t= 1:numel(strt)
        fprintf('.%s=\n', strt{t});
        disp(S.(probeType)(k).(str{p}).(strt{t}));
      end
    else
      if (size(str{p},1)==1)&&isnumeric(S.(probeType)(k).(str{p}))
        fprintf('%s= ', str{p});
        if strncmp(str{p}, 'Time', 4)
          for t= 1:numel(S.(probeType)(k).(str{p}))
            if isnan(S.(probeType)(k).(str{p})(t)); fprintf('NaN ')
            else fprintf('%s ', datestr(S.(probeType)(k).(str{p})(t), ...
                'dd.mm.yyyy HH:MM:SS'));
            end
          end
          fprintf('\n');
        else
          fprintf('%g\n', S.(probeType)(k).(str{p}));
        end
      else
        disp(str{p});
        disp(S.(probeType)(k).(str{p}));
      end
    end
  end
end