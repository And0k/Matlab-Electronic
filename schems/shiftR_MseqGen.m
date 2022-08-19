pol= [1 1 2];
%pol= [1 0 2 1];
L=40; qPower= 3;
Mdegree= numel(pol); T=qPower.^(Mdegree-1)-1; fprintf('\nT= %d:', T);
a= zeros(1, Mdegree-1);  pol2Toend= pol(2:end); %a(1)= 1;
out= NaN(1, L); ain= 1;
for k=1:L;
  a= [ain a(1:(end-1))];
  ain= mod(-sum(pol2Toend.*a), qPower);
  out(k)= a(end); 
  fprintf(1, '%d', out(k));
  if rem(k,T)==0; fprintf(1, ';'); else fprintf(1, ','); end
end
out= out-1;


out1= [ 1, 1, 1,-1, 0, -1, 1, 1, 1, 0]; sum(out1)
out2= [ 1,-1, 1, 1, 1, -1,-1, 1,-1, 0]; sum(out2)

out1= [ 1, 1, 1,-1];
out2= [ 1, 1,-1, 1];
out1= [ 1, 1, 0, 1,-1];
out2= [-1,-1, 0, 1,-1];

out1= [-1, 0, 1, 1,-1];
out2= [ 0, 0, 1,-1, 1];

out1= [-1, 0, 1, 1,-1];
out2= [ 1, 0,-1,-1, 1];

T= numel(out1); %out1= (out1+1)/2; out2= (out2+1)/2;
% out1= repmat(out1,1,5);
% out2= repmat(out2,1,5);
T= numel(out1);
figure(); hold on; plot(xcorr(out1(1:T))+xcorr(out2(1:T)), 'k');
plot(xcorr(out1(1:T)),'b'); plot(xcorr(out2(1:T)),'g');
plot(xcorr(out1(1:T), out2(1:T)), 'r'); 
legend({'AKF1+AKF2','AKF1','AKF2','VKF'}); legend('boxoff');
figure(); plot(conv(out2, out2(1:numel(out1)), 'valid'))
  %  sum
%polyval(pol,1)

c0= 1; %c= 0:(qPower-1)
qPower= 3;
for v=1:qPower
  for jj=1:qPower
    a(jj,v)= mod((jj-1)*(v-1)+c0, qPower);
  end
end