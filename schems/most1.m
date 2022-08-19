R1= 10e3;
Rt= 1e3; %:1.1e3]';
dT= -5:0.1:30; %dT= -dT;
ppm2all= [0 5 10  5   ];   R3_0= 26e3; %10e3;
ppm3all= [0 10 100 100 ];   R4_0= 2.5e3; %910;
%ppmpall= [0  0  0   0     100 100 100 100];  Rp_0= 10e3;
Up= 2.5; 
%------------------------
clr= 'grbkcm'; strLegend= cell(1, numel(ppm2all)-1);
figure('Name', sprintf('Vout(dT) for Rt=%gkOm', Rt/1e3)); hold on;
xlabel(gca, 'dT, degree'); ylabel(gca, 'U, V');
title('corves for resistor''s [R3, R4] different heat tolerance');
Rt= repmat(Rt, size(dT));
R1= repmat(R1, size(Rt));
for k=1:numel(ppm2all)
  ppm2= ppm2all(k);
  ppm3= ppm3all(k);
  %ppmp= ppmpall(k);
  R3= R3_0*(1 + dT*ppm2*1e-6);
  R4= R4_0*(1 + dT*ppm3*1e-6);
%   Rp= Rp_0*(1 + dT*ppmp*1e-6);
% 
%   R_bridge_p= II((R1+Rt)',(R3+R4)')'; % bridge resistance:
%   U= Up*R_bridge_p./(R_bridge_p+Rp);
  Vo= ((Rt./(R1+Rt)) - R4./(R3+R4)).*Up;
  if k==1; Vo_= 0;
  else
    if k<=numel(clr); plot(dT,Vo-Vo_,      clr(k));
    else              plot(dT,Vo-Vo_, ['.' clr(rem(k,numel(clr)))]);
    end
    strLegend{k-1}= sprintf('[%g,%g]ppm', ppm2, ppm3);
  end
end
grid on; legend(strLegend, 'Location' ,'SW'); legend('boxoff');
% figure('Name', 'all dT fo of all values Rt'); hold on;
% plot(Vo');

