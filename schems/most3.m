R1_0= 1e3; R3_0= 1e3; R4_0= 1e3;
Rp= 10e3;
Rt= 1115.9; %e3; %:1.1e3]';
dT= -5:0.1:30; %dT= -dT;
ppm2all= [0  0  0   0    100 10  10  10];   
ppm3all= [0  10 100 100   10  10  100 100];   
ppm1all= [0  5  10  5    100 100 100 10];  
Up= 2.5; 
%------------------------
clr= 'grbkmc'; strLegend= cell(1, numel(ppm2all)-1);
figure('Name', sprintf('Vout(dT) for Rt=%gkOm', Rt/1e3)); hold on;
xlabel(gca, 'dT, degree'); ylabel(gca, 'U, V');
title('corves for resistor''s  [R1, R3, R4]  different heat tolerance');
Rt= repmat(Rt, size(dT));
%R1_0= repmat(R1_0, size(Rt));
for k=1:numel(ppm2all)
  ppm2= ppm2all(k);
  ppm3= ppm3all(k);
  ppm1= ppm1all(k);
  R3= R3_0*(1 + dT*ppm2*1e-6);
  R4= R4_0*(1 + dT*ppm3*1e-6);
  R1= R1_0*(1 + dT*ppm1*1e-6);

  R_bridge_p= II((R1+Rt)',(R3+R4)')'; % bridge resistance:
  U= Up*R_bridge_p./(R_bridge_p+Rp);
  Vo= ((Rt./(R1+Rt)) - R4./(R3+R4)).*U;
  if k==1; Vo_= Vo;
  else
    if k<=numel(clr); plot(dT,Vo, ['.' clr(k)]);
    else              plot(dT,Vo, clr(rem(k,numel(clr))));
    end
    strLegend{k-1}= sprintf('[%g,%g,%g]ppm', ppm1, ppm2, ppm3);
  end
end
grid on; legend(strLegend, 'Location' ,'SW'); legend('boxoff');
% figure('Name', 'all dT fo of all values Rt'); hold on;
% plot(Vo');

