function Pcal=viral_TRHO2P(VirCoef,InputData)
% data: input data matrix
% L: coefficients
format long;
global Tc1 Tc2
%Tc1=375.31;% HFC-161
%Tc2=374.25; %HFC-134a
%Tc2=351.26 ;%HFC-32
R=8.3144;
T=InputData(:,1);
RHO=InputData(:,4);
X1=InputData(:,3);
%{
b1=virial_coeff(1);
b2=virial_coeff(2);
b3=virial_coeff(3);
b4=virial_coeff(4);
b5=virial_coeff(5);
b6=virial_coeff(6);
b7=virial_coeff(7);
b8=virial_coeff(8);
b9=virial_coeff(9);
c1=virial_coeff(10);
c2=virial_coeff(11);
c3=virial_coeff(12);
c4=virial_coeff(13);
c5=virial_coeff(14);
c6=virial_coeff(15);
c7=virial_coeff(16);
c8=virial_coeff(17);
c9=virial_coeff(18);
c10=virial_coeff(19);
c11=virial_coeff(20);
c12=virial_coeff(21);
%}

b1=VirCoef(1);
b2=VirCoef(2);
b3=VirCoef(3);
b4=VirCoef(4);
b5=VirCoef(5);
b6=VirCoef(6);
b7=VirCoef(7);
b8=VirCoef(8);
b9=VirCoef(9);
b10=VirCoef(10);
b11=VirCoef(11);
b12=VirCoef(12);
b13=VirCoef(13);
b14=VirCoef(14);
b15=VirCoef(15);
b16=VirCoef(16);
b17=VirCoef(17);
b18=VirCoef(18);
c1=VirCoef(19);
c2=VirCoef(20);
c3=VirCoef(21);
c4=VirCoef(22);
c5=VirCoef(23);
c6=VirCoef(24);
c7=VirCoef(25);
c8=VirCoef(26);
c9=VirCoef(27);
c10=VirCoef(28);

Tc12=sqrt(Tc1*Tc2);
Tc112=(Tc1*Tc1*Tc2)^(1/3);
Tc221=(Tc1*Tc2*Tc2)^(1/3);
Tr1=T./Tc1;
Tr2=T./Tc2;
Tr12=T./Tc12;
Tr112=T./Tc112;
Tr221=T./Tc221;
%{
B11=b1+b2.*(Tr1.^-1)+b3.*exp(Tr1.^-1);
B22=b4+b5.*(Tr2.^-1)+b6.*exp(Tr2.^-1);
B12=b7+b8.*(Tr12.^-1)+b9.*exp(Tr12.^-1);
B21=b7+b8.*(Tr12.^-1)+b9.*exp(Tr12.^-1);
C111=c1+c2.*(Tr1.^-3)+c3.*(Tr1.^-13);%to normalize change -3,-13 to -5 -12
C222=c4+c5.*(Tr2.^-5)+c6.*(Tr2.^-12);
C112=c7+c8.*(Tr112.^-5)+c9.*(Tr112.^-7);
C211=C112;
C121=C112;
C221=c10+c11.*(Tr221.^-5)+c12.*(Tr221.^-6);%change form 6 to 7
C122=C221;
C212=C221;
%}

B11=b1+b2.*(Tr1.^-1)+b3.*(Tr1.^-2)+b4.*(Tr1.^-3)+b5.*(Tr1.^-6)+b6.*(Tr1.^-8);
B22=b7+b8.*(Tr2.^-1)+b9.*(Tr2.^-2)+b10.*(Tr2.^-3)+b11.*(Tr2.^-6)+b12.*(Tr2.^-8);
B12=b13+b14.*(Tr12.^-1)+b15.*(Tr12.^-2)+b16*(Tr12.^-3)+b17.*(Tr12.^-6)+b18.*(Tr12.^-8);
B21=B12;
C111=c1.*(Tr1.^-5)+c2.*(Tr1.^-6);%to normalize change -3,-13 to -5 -12
C222=c3.*(Tr2.^-5)+c4.*(Tr2.^-6);
C112=c5+c6.*(Tr112.^-5)+c7.*(Tr112.^-6);
C211=C112;
C121=C112;
C221=c8+c9.*(Tr221.^-5)+c10.*(Tr221.^-6);%change form 6 to 7
C122=C221;
C212=C221;

X2=1-X1;
Bm=X1.*X1.*B11+X1.*X2.*B12+X2.*X1.*B21+X2.*X2.*B22;
Cm=X1.*X1.*X1.*C111+X1.*X1.*X2.*C112+X1.*X2.*X1.*C121+X1.*X2.*X2.*C122+X2.*X1.*X1.*C211+X2.*X1.*X2.*C212+X2.*X2.*X1.*C221+X2.*X2.*X2.*C222;
Z=1+Bm.*RHO+Cm.*(RHO.^2);
Pcal=Z.*RHO.*T*R;
end