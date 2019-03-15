function curvefit_viral
%% �ֲ��Ż���Ϻ���
%NP=140; Number of Points
clear
%����Ϊ�������ͣ�
format long;
%�������峣��
R=8.3144;
% �����ٽ�ֵ
global Tc1 Tc2
% Tc1=304.1282;% CO2
% Tc2=369.825; %Propane
Tc1=369.825;%Propane
Tc2=396.44;%trifluoroiodomethane
%Tc1=375.31;% HFC-161
%Tc2=374.25; %HFC-134a
%Tc2=351.26 ;%HFC-32
%% ȫ����0 �����������
coeff=0;%��0
RHOcal=0;
rho=0;
RHOC=0;
Dev=0;
Dev_P=0;
Dev_RHO=0;

%% �˻��������
Readfile=input('which file do you want to use:\n ������������Ҫʹ���ļ����ļ���������.csv):\n','s');
Readedfile=[Readfile,'.csv'];
%LP=21;
try
datafile=csvread(Readedfile);
catch
datafile=xlsread([Readfile,'.xlsx']);
end
[NP,N]=size(datafile);
M=['��������ļ�����',num2str(NP),'�����ݵ㣬ÿ�����ݵ㹲',num2str(N),'��ά��'];
disp(M);

%��5��ά�ȵ�������ȫ����ȡ������á�
Zexp=datafile(:,5);
T=datafile(:,1);
P=datafile(:,2);
rho=datafile(:,4);
X1=datafile(:,3);
X2=1-X1;
%����virial���̵��������
data=datafile(:,1:4);
%X1=data(:,3);
%X2=1-X1;
askoperation=input('what do you want to do? Optimize(Y)/caculate(C or other random symbol):\n��������ϻ��Ǽ��㣬�������Y����������C)\n','s');
%% read coeffcient from file viralco.xls
coeff=csvread('viralco.csv');
coeff=coeff';
%disp('read done\n(ά��ϵ����ȡ��ϣ�\n');
%% ȫ���Ż�
if askoperation=='Y'
% ʹ����ϳ��򣬿�ʼ���
disp('��ʼ��ϣ��԰�����......');
%coeff0=coeff;
options=optimset('TolFun',1e-10);
lb=[];
ub=[];
%L=lsqcurvefit(@viral3,L0,data,Zexp,lb,ub,options);
% �鿴���-Pcal
%Zcal=viral3(L,data)
%Dev=(Zcal-Zexp)./Zexp
model=@viral_TRHO2P;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',data,'ydata',P,'x0',coeff,'lb',lb,'ub',ub);
%���ģ���� model �ж��壬������data���壬���P,�����ĳ�ʼֵ��coeff0
ms=MultiStart;
%coeff_fin�����յĴ�ȡֵ
[coeff_fin,fval,exitflag,output,solutions]=run(ms,problem,50);
coeff=coeff_fin;
%disp('optimize done');
disp('������');
else
    %����ϣ���������һ�ε��Ż�
    disp('��������ϣ�ֱ��ʹ��ϵ������');
end
%% caculate Pressure/kPa
disp('���㿪ʼ,�Ե�Ƭ��......')
P_cal=viral_TRHO2P(coeff,data);
%Z_cal=P_cal./(rho.*R*T);
Dev_P=(P_cal-P)./P;
Dev_P=100*Dev_P;
%L=L1;
%% ά��ϵ��
% ��������η��̵ĽǶ���P X T�����RHO
b1=coeff(1);
b2=coeff(2);
b3=coeff(3);
b4=coeff(4);
b5=coeff(5);
b6=coeff(6);
b7=coeff(7);
b8=coeff(8);
b9=coeff(9);
c1=coeff(10);
c2=coeff(11);
c3=coeff(12);
c4=coeff(13);
c5=coeff(14);
c6=coeff(15);
c7=coeff(16);
c8=coeff(17);
c9=coeff(18);
c10=coeff(19);
c11=coeff(20);
c12=coeff(21);
Tc12=sqrt(Tc1*Tc2);
Tc112=(Tc1*Tc1*Tc2)^(1/3);
Tc221=(Tc1*Tc2*Tc2)^(1/3);
Tr1=T./Tc1;
Tr2=T./Tc2;
Tr12=T./Tc12;
Tr112=T./Tc112;
Tr221=T./Tc221;
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
Bm=X1.*X1.*B11+X1.*X2.*B12+X2.*X1.*B21+X2.*X2.*B22;
Cm=X1.*X1.*X1.*C111+X1.*X1.*X2.*C112+X1.*X2.*X1.*C121+X1.*X2.*X2.*C122+X2.*X1.*X1.*C211+X2.*X1.*X2.*C212+X2.*X2.*X1.*C221+X2.*X2.*X2.*C222;
%% solve root of density
PRT=P./(R*T); % P/RT
% ���������η��̵�ϵ������
I=ones(1,NP); 
I=I';
tri_coeff=[Cm,Bm,I,-PRT];
% 2018-11-07-����
% ����ԭ�������
% ��������������η��̴����������������ϵ�ʵ�ʸ�����ʱ�Ա�����һ��������
% ʵ��Ϊ��ȷֵ�ǷǺ����
% 
% ������Ҫ��
%     1������жԳ�ʼ��ʵ��ֵ����������
%     2�����ɹ������������Ӷ�
%     3�����˼·���£�
%          �ڶ��ʵ���б�Ȼ��һ��ֵ����ӽ�ʵ��ֵ�ģ����ǲ�����ʵ��ֵ��Ϊ�Ƚ϶����ɸѡ��
%          ���������ʧȥͨ���ԣ������Ҫ����PR���̵ļ���ֵΪԤ��ֵ�����ֵ�Ƿǳ��ӽ�ʵ
%          �ʵĶ�����ʹ��ѡ��viral�������Ķ��ʵ������ӽ�PR����Ԥ��ֵ������
%     4) ����PR������Ҫ���±�дPR�ĳ��򣬲����ã�ֱ����ȡ���ʵ�ʸ��нӽ�ʵ��ֵ��
rhoj=[];    %����һ���յ��м�������������洢���η��̵Ķ���ʵ�� 
deltaRHOj=[];
% ��ÿһ�����̿���ɸѡʵ�ʸ�
for i=1:NP
    j=0;
    s=roots(tri_coeff(i,:));
    for k=1:length(s);
        ss=isreal(s(k));
        if (ss==1)&&(s(k)>0) 
            j=j+1;
            rhoj(j)=s(k);
            deltaRHOj(j)=rhoj(j)-rho(i);   
           % rhoc=s(k);
        else
            %pass
        end
    end
[deltamin,index]=min(abs(deltaRHOj(1:j)));
delta=deltaRHOj(index);
rhoc=rho(i)+delta;
RHOC(i)=rhoc;
end
RHOcal=RHOC';
Dev_RHO=100*(RHOcal-rho)./rho;
%Dev3=-100*Dev3;
disp('�������󣬿�ʼ��ͼ......');
figure(1);
plot(T,Dev_P,'o');
title('Pressure devitation');
figure(2);
plot(T,Dev_RHO,'o');

title('density devitation');
%% plot
figure(3);
plot(T,1000*Bm,'o');
title('Bm*1000');
%% caculate p and density
%%  use NIST to caculate critical point
%disp('Caculate done now save coefficents')
disp('(�����������ȡά��ϵ��)');
%% save coefficents
LT=coeff';
csvwrite('viralco.csv',LT);
SaveYN=input('�Ƿ��ȡ����ѹ��ƫ��ܶ�ƫ������(Y|N):','s');
if SaveYN=='Y'
    csvwrite('devtation.csv',[T P X1 rho RHOcal P_cal Dev_P Dev_RHO]);
end
end


function P=viral_TRHO2P(L,data)
% data: input data matrix
% L: coefficients
global Tc1 Tc2
%Tc1=375.31;% HFC-161
%Tc2=374.25; %HFC-134a
%Tc2=351.26 ;%HFC-32
R=8.3144;
T=data(:,1);
RHO=data(:,4);
X1=data(:,3);
b1=L(1);
b2=L(2);
b3=L(3);
b4=L(4);
b5=L(5);
b6=L(6);
b7=L(7);
b8=L(8);
b9=L(9);
c1=L(10);
c2=L(11);
c3=L(12);
c4=L(13);
c5=L(14);
c6=L(15);
c7=L(16);
c8=L(17);
c9=L(18);
c10=L(19);
c11=L(20);
c12=L(21);
Tc12=sqrt(Tc1*Tc2);
Tc112=(Tc1*Tc1*Tc2)^(1/3);
Tc221=(Tc1*Tc2*Tc2)^(1/3);
Tr1=T./Tc1;
Tr2=T./Tc2;
Tr12=T./Tc12;
Tr112=T./Tc112;
Tr221=T./Tc221;
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
X2=1-X1;
Bm=X1.*X1.*B11+X1.*X2.*B12+X2.*X1.*B21+X2.*X2.*B22;
Cm=X1.*X1.*X1.*C111+X1.*X1.*X2.*C112+X1.*X2.*X1.*C121+X1.*X2.*X2.*C122+X2.*X1.*X1.*C211+X2.*X1.*X2.*C212+X2.*X2.*X1.*C221+X2.*X2.*X2.*C222;
Z=1+Bm.*RHO+Cm.*(RHO.^2);
P=Z.*RHO.*T*R;
end