function boundsT= findSeparatorConstr(DATA, freq)
  bInd= isnan(DATA.Direction)|isnan(DATA.Vabs);
  sumBad= sum(bInd);
  zmin= min(max(5,sum(bInd))*dayHz/freq, 1/24); %min detected interval range = [5/ID.fs 1Hour]
  boundsT= [findSeparator(DATA.Time, zmin, sumBad); numel(DATA.Time)+1];
end