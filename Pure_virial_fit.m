%function Pure_virial_fit

%propane
PROPANE=xlsread('PROPANE_CURVE.xlsx');
CF3I=xlsread('CF3I_CURVE.xlsx');
%set propane data
T_propane=PROPANE(:,1);
VIRB_propane=PROPANE(:,2);
VIRC_propane=PROPANE(:,3);
%set CF3I data
T_CF3I=CF3I(:,1);
VIRB_CF3I=CF3I(:,2);
VIRC_CF3I=CF3I(:,3);
%critical point
Tcrit_propane=369.825;%Propane
Tcrit_CF3I=396.44;%CF3I
%B11=b1+b2.*(Tr1.^-1)+b3.*(Tr1.^-2)+b4.*(Tr1.^-3)+b5.*(Tr1.^-6)+b6.*(Tr1.^-8);
%C111=c1.*(Tr1.^-5)+c2.*(Tr1.^-6);

XB_propane=[ones(size(T_propane)),(T_propane/Tcrit_propane).^-1,(T_propane/Tcrit_propane).^-2,(T_propane/Tcrit_propane).^-3, ...
    (T_propane/Tcrit_propane).^-6,(T_propane/Tcrit_propane).^-8];
XC_propane=[(T_propane/Tcrit_propane).^-5,(T_propane/Tcrit_propane).^-6];
VIRB_coeff_pr=XB_propane\VIRB_propane
VIRC_coeff_pr=XC_propane\VIRC_propane

XB_CF=[ones(size(T_CF3I)),(T_CF3I/Tcrit_CF3I).^-1,(T_CF3I/Tcrit_CF3I).^-2,(T_CF3I/Tcrit_CF3I).^-3,(T_CF3I/Tcrit_CF3I).^-6, ...
    (T_CF3I/Tcrit_CF3I).^-8];
XC_CF=[(T_CF3I/Tcrit_CF3I).^-5,(T_CF3I/Tcrit_CF3I).^-6];
VIRB_coeff_CF=XB_CF\VIRB_CF3I
VIRC_coeff_CF=XC_CF\VIRC_CF3I
xlswrite('virial_propane.xlsx',[VIRB_coeff_pr;VIRC_coeff_pr])
xlswrite('virial_CF3I.xlsx',[VIRB_coeff_CF;VIRC_coeff_CF])
%end