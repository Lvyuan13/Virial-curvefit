%% �ֲ��Ż���Ϻ���
%NP=140;
clc 
clear
%% ȫ����0 �����������
L=0;%��0
RHOcal=0;
rho=0;
RHOC=0;
Dev=0;
Dev2=0;
Dev3=0;
global Tc1 Tc2
%% �˻��������
R=input('which file do you want to use:\n ������������Ҫʹ���ļ����ļ���������.xlsx):\n','s');
RE=[R,'.xlsx'];
LP=21;
A=xlsread(RE);
[NP,N]=size(A);
M=['��������ļ�����',num2str(NP),'�����ݵ㣬ÿ�����ݵ㹲',num2str(N),'��ά��'];
disp(M);
% Tc1=304.1282;% CO2
% Tc2=369.825; %Propane
Tc1=369.825;
Tc2=396.44;
%Tc1=375.31;% HFC-161
%Tc2=374.25; %HFC-134a
%Tc2=351.26 ;%HFC-32
Zexp=A(:,5);
T=A(:,1);
P=A(:,2);
rho=A(:,4);
data=A(:,1:4);
R=8.3144;
R2=input('what do you want to do? Optimize(Y)/caculate(C or other random symbol):\n��������ϻ��Ǽ��㣬�������Y����������C)\n','s');
%% read coeffcient from file viralco.xls
L=xlsread('viralco.xls');
L=L';
%disp('read done\n(ά��ϵ����ȡ��ϣ�\n');
%% ȫ���Ż�
if R2=='Y'
%L=0; % to ensure the origin L didn't influence
%L0=ones(1,LP).*0.3;
disp('��ʼ��ϣ��԰�����');
L0=L;
options=optimset('TolFun',1e-10);
lb=[];
ub=[];
%L=lsqcurvefit(@viral3,L0,data,Zexp,lb,ub,options);
% �鿴���-Pcal
%Zcal=viral3(L,data)
%Dev=(Zcal-Zexp)./Zexp
model=@viral3;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',data,'ydata',Zexp,'x0',L,'lb',lb,'ub',ub);
ms=MultiStart;
[L1,fval,exitflag,output,solutions]=run(ms,problem,50);
L=L1;
disp('optimize done');
disp('������');
else
    disp('��������ϣ�ֱ��ʹ��ϵ������');
end
%% caculate Pressure/kPa
disp('���㿪ʼ,�Ե�Ƭ��')
Zcal2=viral3(L,data);
Pcal2=Zcal2.*rho.*R.*T;
Dev2=(Pcal2-P)./P;
Dev2=100*Dev2;

%L=L1;
%% ά��ϵ��

X1=data(:,3);
X2=1-X1;
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
Bm=X1.*X1.*B11+X1.*X2.*B12+X2.*X1.*B21+X2.*X2.*B22;
Cm=X1.*X1.*X1.*C111+X1.*X1.*X2.*C112+X1.*X2.*X1.*C121+X1.*X2.*X2.*C122+X2.*X1.*X1.*C211+X2.*X1.*X2.*C212+X2.*X2.*X1.*C221+X2.*X2.*X2.*C222;
%% solve root of density
PRT=P./(R*T); % P/RT
D=ones(1,NP); % ���η���ϵ��
D=D';
coeff=[Cm,Bm,D,-PRT];
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
for i=1:NP
    j=0;
    s=roots(coeff(i,:));
    for k=1:length(s);
        ss=isreal(s(k));
        if (ss==1)&&(s(k)>0) 
            j=j+1;
            rhoj(j)=s(k);
            deltaRHOj(j)=rhoj(j)-rho(i);   
           % rhoc=s(k);
        else
            
        end
    end
[deltamin,index]=min(abs(deltaRHOj(1:j)));
delta=deltaRHOj(index);
rhoc=rho(i)+delta;
RHOC(i)=rhoc;
end
RHOcal=RHOC';
Dev3=100*(RHOcal-rho)./rho;
%Dev3=-100*Dev3;
disp('�������󣬿�ʼ��ͼ');
figure(1);
plot(T,Dev2,'o');
title('Pressure devitation');
figure(2);
plot(T,Dev3,'o');

title('density devitation');
%% plot
figure(3);
plot(T,1000*Bm,'o');
title('Bm*1000');
%% caculate p and density
%%  use NIST to caculate critical point

disp('Caculate done now save coefficents')
disp('(�����������ȡά��ϵ��)');
%% save coefficents
LT=L';
xlswrite('viralco.xls',LT);