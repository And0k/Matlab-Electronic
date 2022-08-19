R1= 10e3;
Rt= 1e3; %:1.1e3]';
dT= 0:0.1:30; %dT= -dT;
ppm2all= [ 5 10 10  5     5 10 10  5];   R2_0= 10e3;
ppm3all= [10 10 100 100   10 10 100 100];   R3_0= 910;
ppmpall= [0  0  0   0     100 100 100 100];  Rp_0= 10e3;
Up= 10; 
%------------------------
clr= 'rgbkcm'; strLegend= cell(1, numel(ppm2all));
figure('Name', sprintf('Vout(dT) for Rt=%gkOm', Rt/1e3)); hold on;
xlabel(gca, 'dT, degree'); ylabel(gca, 'U, V');
title('corves for resistor''s different heat tolerance');
Rt= repmat(Rt, size(dT));
R1= repmat(R1, size(Rt));
for k=1:numel(ppm2all)
  ppm2= ppm2all(k);
  ppm3= ppm3all(k);
  ppmp= ppmpall(k);
  R2= R2_0*(1 + dT*ppm2*1e-6);
  R3= R3_0*(1 + dT*ppm3*1e-6);
  Rp= Rp_0*(1 + dT*ppmp*1e-6);

  R_bridge_p= II((R1+Rt)',(R2+R3)')'; % bridge resistance:
  U= Up*R_bridge_p./(R_bridge_p+Rp);
  Vo= ((Rt./(R1+Rt)) - R3./(R2+R3)).*U;
  if k<=numel(clr); plot(dT,Vo,      clr(k));
  else              plot(dT,Vo, ['.' clr(rem(k,numel(clr)))]);
  end
  strLegend{k}= sprintf('R2,R3,Rp: [%g %g %g]ppm', ppm2, ppm3, ppmp);
end
grid on; legend(strLegend, 'Location' ,'SW'); legend('boxoff');
% figure('Name', 'all dT fo of all values Rt'); hold on;
% plot(Vo');

