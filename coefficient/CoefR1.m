%Coefficients will be used in this order:
%1) DATA= DATA*kFile2Counts;
%2) DATA(n)= k(n,1)+DATA(n)^2*k(n,2)+DATA(n)^2*k(n,3)+ ...

%Used In Cols (in comfortable order):
SLoad.tline= 'T1 T3 T4 T2';
%Out Cols (in same order):
SLoad.tlineReplase= 'P_(dBar) Ttop_(°C) Tmid_(°C) Tbot_(°C) Time';
fd= 0.1;

%Coefficient for convert to counts
k= zeros(5,4);

%Table of coefficients for convert counts to phisic:
%Ch{[Ind]}= '[ChName]';	k([Ind],:)=	[k0	k1	k2];
%Ind 	- Any Indexes, different for different rows
%ChName 	- Consist of some first characters of tline_ChNameReplase, to find  column for process.
%it mast be identical to same number of first characters of only one tline_ChNameReplase word
 
Ch{1 }= 'Ttop';	k(1 ,:)=	[-1.775 0.0004946   -1.768e-010 1.328e-014];
Ch{2 }= 'Tmid';	k(2 ,:)=	[-26.64 0.0008962	1.671e-009  7.579e-016];
Ch{3 }= 'Tbot';	k(3 ,:)=	[-4.795 0.0004953   -1.5e-010   1.567e-014];
Ch{4 }= 'P_';	k(4 ,:)=	[-0.62  0.03858	0	0];	%[-32.219	0.076991];        

CalcData(1).tline= 'DepTermCalc_(m)'; %sw_salt(C_/sw_c3515,T_,P_);
L0= 20; %probe length
%L2= (L0 - (23 - P_))
CalcData(1).function= @(Ttop_, Tmid_, Tbot_, P_) P_ + L0*(Ttop_-Tmid_)./(Ttop_-Tbot_);

%Ch{1 }= 'Ttop';	k(1 ,:)=	[-2.8	0.000537	0];
%Ch{2 }= 'Tmid';	k(2 ,:)=	[-28.18	0.001041	0];
%Ch{3 }= 'Tbot';	k(3 ,:)=	[-6.527	0.0005606	0];
%Ch{4 }= 'P_';	k(4 ,:)=	[-0.62	0.03858	0];	%[-32.219	0.076991]; 

SLoad.inDate= []; SLoad.DateForm= {[] 'mmddHHMM' 1/fd}; %format of file name date