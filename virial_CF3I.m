function virial_CF3I
% check out coefficenit of B and C
% input: T arrary
% 该方程的所有系数均采用自段远源教授博士学位论文-PVTX性质部分
data=xlsread('CF3I_3.xlsx');
R=8.3144;
Rg=0.0424402;
T=data(:,1);
P=data(:,3);%kPa
RHOm=data(:,4);%kg/m^3
RHO=data(:,5);%mol/dm^3
Tc=396.44;%critical point of trifluoroiodomethane
M=195.91; %molar mass of trifluoiodomethane
%{
unit;
B:
C:
%}
Bm0=0.1017882;
Bm1=-0.3547075;
Bm2=0.4369783;
Bm3=-0.1979620;
Bm4=0.01492476;
Bm5=-2.31700e-3;
Cm0=-1.14023e-6;
Cm1=9.84870e-7;
D0=8.98114e-9;
D1=-8.64263e-9;

Tr=T./Tc;
% mass unit
Bm=Bm0+Bm1*Tr.^-1+Bm2*Tr.^-2+Bm3*Tr.^-3+Bm4*Tr.^-6+Bm5*Tr.^-8;
Cm=Cm0*Tr.^-5+Cm1*Tr.^-6;
% Dm=D0+D1*Tr;
% using three order truncated form
Dm=0;
% mole unit
B0=Bm0*M
B1=Bm1*M
B2=Bm2*M
B3=Bm3*M
B4=Bm4*M
B5=Bm5*M
C0=Cm0*(M^2)
C1=Cm1*(M^2)
B=Bm*M;
C=Cm*(M^2);
D=Dm*(M^3);
D=0;
Pcal=Rg*RHOm.*T.*(1+Bm.*RHOm+Cm.*RHOm.^2+Dm.*RHOm.^3);
DevP=100*(Pcal-P)./P;
%{
Pcal=R*RHO.*T.*(1+B.*RHO+C.*RHO.^2+D.*RHO.^3);
DevP=100*(Pcal-P)./P;
%}
figure(1)
title('Bm')
plot(T,Bm,'+')

figure(2)
title('Cm')
plot(T,Cm,'*')

figure(3)
title('压力偏差')
plot(T,DevP,'+')
CF3IFILE=table(T,P,RHO,Pcal,DevP,B,C);
writetable(CF3IFILE,'CF3I_3_virial.csv');
end