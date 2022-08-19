function SCoef= coefLoad(dataTime, probeType, Nprobe)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Example

%SCoef= coefLoad(datenum('140617','yymmdd'),'Inclinometr',6);
%coefDisp(SCoef) % Show loaded coeficients
SCoef= struct();
if Nprobe>0&&exist('d:\Work\MatlabHistory\mat\coef_Electronics.mat','file')
  if any(strncmpi(probeType, {'i','#'}, 1)); probeType= 'Inclinometer'; end
  SCoef= load('d:\Work\MatlabHistory\mat\coef_Electronics.mat', probeType);
  if ~isfield(SCoef.(probeType),'i'); stopHere('Bad coefficients file!');
  end
  bInd= [SCoef.(probeType).i]==Nprobe;
  k= sum(bInd);
  while k>1
    if k>1;
      k= find(bInd);
      for t= k
        bInd(t)= ~((dataTime(1)< SCoef.(probeType)(t).TimeRange(1))...
          ||(dataTime(1)>=SCoef.(probeType)(t).TimeRange(2))); %invesion for right NaN proc
      end
      %       if sum(bInd)==1;
      %         k= k(bInd); bInd= false(size(SCoef.(probeType))); bInd(k)= true;
      %       end
      if ~any(bInd)
        p= str2double(input(sprintf(['No coef for probe#%d for data %s.\n', ...
          'If use some of %d existed for this probe then input it sequence number:'], ...
          Nprobe, datestr(dataTime(1)), numel(k)),'s'));
        if 0<p&&p<=numel(k); bInd(:)= false; bInd(k(p))= true; k=1;
        else k= 0;
        end
      else k= sum(bInd);
      end
    end
    if k==0
      p= str2double(input(sprintf(['No coef for probe#%d for data %s.\n', ...
        'If use other # input it:'], Nprobe, datestr(dataTime(1))),'s'));
      bInd= [SCoef.Inclinometer.i]==p;
      k= sum(bInd);
    elseif k>1
      p= str2double(input(sprintf(['More than 1 coef (%d) for probe#%d.\n', ...
        'Input used sequence #:'], k, Nprobe),'s'));
      if 0<p&&p<=k; k= find(bInd); bInd(:)= false; bInd(k(p))= true; k=1;
      else k= 0;
      end
    end
  end
  if k==1
    fprintf('\nUse coefficients for probe %s#%d obtained %s ', ...
      probeType, Nprobe, datestr(SCoef.(probeType)(bInd ...
      ).TimeProcessed, 'dd.mm.yyyy HH:MM'));
    SCoef= SCoef.(probeType)(bInd);
  else fprintf('\nNo coef for probe#%d', Nprobe);   
  end
end

