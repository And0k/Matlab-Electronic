%%
U= 5;
R1= 30;
R2= 10;
Umid= U*R2/(R1+R2)

%R1/R2= 5/1.25 - 1;

%%
R1= 51000;
R2= 1000 + II(100,100000);
K= 2+2*R1/R2

%%
R1= 5600;
R2= 560;
R7= 2000;
Up= 5;

I1= Up/(R1+R7)
U7= R7*I1
I=  U7/R2
%%
PresUnit_kg2m= 1.01972e-004/(0.0122980*6.4516e-004)
%% ADUC831 TIMER3 UART
fUART= 4800; fcore= 12000000;
DIV= floor(log(fcore/(32*fUART))/log(2))
T3FD= round(2*fcore/((2^DIV) * fUART)) - 64