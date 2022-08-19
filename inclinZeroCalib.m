str= ...
  'd:\WorkData\Cruises\_BalticSea\140301\_source\inclinometr\inclinometr.h5';
  'd:\WorkData\Cruises\_BlackSea\1309\_source\test_Incl#7+Vector\Inclinometr';
SLoad= struct('SaveDirName','', 'strDir',str, ...
  'bNotSaveTimeString',true, 'noTXT',true); % c:\TEMP\130609.h5
%SLoad.noTXT= true; %SLoad.noTXT= false;  SLoad.chunkDays= '';
str= {'i3','i4','i7','i8'}; %i4&i7 with P. str= {'WaveGage1'}; str= {'131003#i7.txt'}; %'i1
%% Period of zero velocity to load:
FilterDat= struct();
FilterDat.Time.min= datenum('05.03.2014 06:00:00', 'dd.mm.yyyy HH:MM:SS');
FilterDat.Time.max= datenum('05.03.2014 10:00:00', 'dd.mm.yyyy HH:MM:SS');
% FilterDat.Time.max= now(); FilterDat.Time= rmfield(FilterDat.Time, 'min'); 
%'07.07.2013 00:00:00' %FilterDat.Time.max= datenum('09.06.2013 23:19:59'
%'26.06.2013 11:00:00' %FilterDat.Time.max= datenum('26.06.2013 23:19:59'
%'10.05.2013 19:00:00' %FilterDat.Time.max= datenum('10.05.2013 23:30:00'
%'27.09.2013 14:43:43'
 for k= 1:numel(str)
  SLoad.prefix= regexp(str{k}, '[a-zA-Z_]*','match','once');
  Nprobe= str2double(regexp(str{k}(numel(SLoad.prefix):end), '\d+', 'match','once'));
  if SLoad.prefix(1)=='i'; SLoad.prefix= 'Inclinometer'; end
  [a1D, DATA, ax, SL]= inclinometr2([SLoad.strDir bIf(SLoad.strDir(end)=='5','/','\') str{k}], FilterDat, SLoad);

    %% Zero calibration
    %   a1D= FSFilter(a1D, FilterDat, {'Gx','Gy','Gz'} ,0, 'Time'); %set NaNs in G
    %   bGood= true(size(DATA.Gx));
    %   n= numel(bGood)
    %   iSt= 55:150:n;
    %   for iIn=0:20; bGood(iSt+ iIn)= false; end
    %   bGood(n+1:end)= [];
    %   plot(ax(2), DATA.Time, bGood, 'y');
    a1D.Vabs(imag(a1D.Vabs)~=0)= NaN;
    a1D.Pitch= medfiltSmooth(a1D.Pitch, [5 55]);
    a1D.Roll = medfiltSmooth(a1D.Roll , [5 55]);
    a1D.Vabs = medfiltSmooth(a1D.Vabs , [5 55]);
    plot(ax(2), DATA.Time, a1D.Vabs,  'c');
    plot(ax(4), DATA.Time, a1D.Pitch, 'k');
    plot(ax(4), DATA.Time, a1D.Roll , 'c');
    bGood= true(size(DATA.Gx));
    if false
      %figure(); hist(a1D.Vabs,100);
      [f,v]= hist(a1D.Vabs,100);
      p= find(f/max(f)>0.1,1);
      bGood= (v(p-1)<=a1D.Vabs)&(a1D.Vabs<=nanmean(a1D.Vabs));
      plot(ax(2), DATA.Time, bGood, 'y');
      bGood(407992:end)= false;

      p=-nanmean(double(a1D.Pitch(bGood)))*pi/180; mode(a1D.Vabs)
      r= nanmean(double(a1D.Roll( bGood)))*pi/180;
    else

      p=-double(mode(a1D.Pitch(bGood)))*pi/180;
      r= double(mode(a1D.Roll( bGood)))*pi/180;
    end
    %Remove founded zero offset and save new coef
    fprintf('Add offset to zero: Pitch= %gdeg, Roll= %gdeg', p*180/pi, r*180/pi);
    Rp= [cos(p) 0 -sin(p); 0 1 0; sin(p) 0 cos(p)];
    Rr= [1 0 0; 0 cos(r) sin(r); 0 -sin(r) cos(r)];
 
    curInc= SL.coef;
    curInc.G.A= Rr*Rp*SL.coef.G.A;
    curInc.TimeRange= [fix(a1D.Time(1)), NaN]; %'10.05.2013 11:00:00'
    %curInc.TimeRange= [datenum('09.06.2013 12:00:00', 'dd.mm.yyyy HH:MM:SS')-dayHz, NaN];
    %curInc.TimeRange= [datenum('01.01.2013 00:00:00', 'dd.mm.yyyy HH:MM:SS'), NaN];
    %[DATA.Time(1), NaN];
    curInc.TimeProcessed= now();

    SCoef= load('d:\Work\MatlabHistory\mat\coef_Electronics.mat');
    bInd= [SCoef.(SLoad.prefix).i]==Nprobe; %dispCoef(SCoef.(SLoad.prefix)(end), 'Inclinometr')
    k= find(bInd);
    for p= k
      b= curInc.TimeRange(1) < SCoef.(SLoad.prefix)(p).TimeRange(1);           %p is next interval
      if b && ~(curInc.TimeRange(2)>SCoef.(SLoad.prefix)(p).TimeRange(1)) %save start of closest interval as end of this
        curInc.TimeRange(2)= SCoef.(SLoad.prefix)(p).TimeRange(1)-dayHz;
      elseif ~(curInc.TimeRange(1)>(SCoef.(SLoad.prefix)(p).TimeRange(2))) %p intersect this
        if ~isfinite(SCoef.(SLoad.prefix)(p).TimeRange(2)); s= 'inf';
        else s= datestr(SCoef.(SLoad.prefix)(p).TimeRange(2), 'dd.mm.yyyy HH:MM:SS');
        end
        stopHere('Change coef%d: %s-\p%s\nend limit (to be before this) to\p%s', ...
          p, datestr(SCoef.(SLoad.prefix)(p).TimeRange(1), 'dd.mm.yyyy HH:MM:SS'), ...
          s, datestr(curInc.TimeRange(1)-dayHz, 'dd.mm.yyyy HH:MM:SS'));
        SCoef.(SLoad.prefix)(p).TimeRange(2)= fix(curInc.TimeRange(1))-1; %curInc.TimeRange(1)= datenum('27.02.2014 00:00:00', 'dd.mm.yyyy HH:MM:SS');
      end
      %intersect remains?
      bInd(p)= ~(b||(a1D.Time(1)>=SCoef.(SLoad.prefix)(p).TimeRange(2))); %invesion for right NaN proc
    end
    fprintf('\nSet coef range:\p%s', datestr(curInc.TimeRange(1),'dd.mm.yyyy HH:MM:SS'));
    if ~isnan(curInc.TimeRange(2));
      fprintf(' -\p%s', datestr(curInc.TimeRange(2),'dd.mm.yyyy HH:MM:SS'));
    else fprintf(' -\tNaN');
    end
    if any(bInd); stopHere('intersect remains?');
    end
    n= find([SCoef.(SLoad.prefix).i]>0,1,'last')+1;
    stopHere('Write to empty item %d by i=#%d: ', n, Nprobe);
    curInc.comment= 'accelerometer zeroing';
%     for p= k  %add new item, change latence of others
%       if curInc.TimeRange(1)<=SCoef.(SLoad.prefix)(p).TimeRange(2)
%       else
%         if ~isfinite(curInc.TimeRange(2)); s= 'inf';
%         else s= datestr(curInc.TimeRange(2));
%         end
%         stopHere('Change End Limit of this coef %d from to %s', p, ...
%           s, datestr(SCoef.(SLoad.prefix)(p).TimeRange(1)-dayHz));
%         curInc.TimeRange(2)= SCoef.(SLoad.prefix)(p).TimeRange(1)-dayHz;
%       end
% 
%     else %replace item? which?
%       %if last item
%       n= find(bInd, 1, 'last');
%       fprintf('Write to existed item %d by i=#%d: ', n, Nprobe);
%       stopHere()
%     end
    SCoef.(SLoad.prefix)(n)= curInc; %save2coefFile(curInc, 'next', 'Inclinometer');
    save('d:\Work\MatlabHistory\mat\coef_Electronics.mat', '-struct', 'SCoef');
       % , '-append');
end
fprintf('\nfinished zeroing\n');
if false %correct NaN time in prev data
  SCoef= load('d:\Work\MatlabHistory\mat\coef_Electronics.mat');
  SLoad.prefix= 'Inclinometer';
  Nprobe= 4;
  bInd= [SCoef.(SLoad.prefix).i]==Nprobe;
  %for k=1:numel(SCoef.(SLoad.prefix)); SCoef.(SLoad.prefix)(k).TimeRange= SCoef.(SLoad.prefix)(k).TimeRange'; end
  k= find(bInd);
  p= [SCoef.(SLoad.prefix).TimeRange];
  iBad= find(isnan(p(2,k))&k<max(k)); iBad= k(iBad);
  if ~isempty(iBad)
  p(2,iBad)= min(p(1,k(p(1,k)>p(1,iBad))));
  SCoef.(SLoad.prefix)(iBad).TimeRange(2)= p(2,iBad)-dayHz;
  datestr(SCoef.(SLoad.prefix)(iBad).TimeRange)
  end
end