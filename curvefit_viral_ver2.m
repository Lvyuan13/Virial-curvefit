function curvefit_viral_ver2
%% �ֲ��Ż���Ϻ���
%NP=140; Number of Points
%����Ϊ�������ͣ�
% version 2
% 2019-03-20 ����
% ���������ά��ϵ�������ٲ��ó¹������ڿ�����Ķ���ʽ��ϣ����Ǹ��ݻ���ԭ��
% �������µ����Զ���ʽ���������幫ʽ�����ģ��Ѿ�������ϵ����21��������28��
% ��������Ҫ��һ�㣬�����ô���������ϵķ������ά��ϵ��������������֤����ǰʹ�õ�
% ���벿�ִ�����������ϵķ�����Ȼ�����൱�ĺ����ԣ���ϵͳ����������ȫ����Ӧ�õģ�Ҳ����һ����
% �����ԣ����ǣ��������㷢�֣���ǰ�ķ������������ά��ϵ��Bmʱ�������ߣ������ڼ�������ά��ϵ��
% ʱ���֣��������Cm�������Ƶ������Ĳ�಻С,������С��Χ�ڻ����Cm������Ϊ��ֵ��������������ǲ���ȷ��
% ���ڲ��õ��·�����
% 1. propane��ά��ϵ��Bm Cm ʹ��refprop �����ʹ�� Helmholtzģ�͵ļ���ֵ�����������������Ѿ����ϣ��ٸ���Bm Cm 
% ����ֵ��Bm�� Cm ��ϳ�Tr �Ķ���ʽ�� ������ʱ���Ӻ��������
% 2. CF3I ��ά��ϵ�������ö�ԶԴ���ڵķ��̼������Ҳ���ȸ��ݼ������Bm Cm���ٽ�Bm��Cm ��ϳ�Tr�Ķ���ʽ��
% ͨ�����㷢�֣��÷���ϳ���Сϵ�� B1 B2 B3 etc ���������н����ȫһ�£���˿����ж�����ʱ��ȫ��ȷ��
% 3. ����������������������Ĵ�����ϵ��ֱ��Ƕ�뷽�̣�������Ҫ��ϵ�����ֻ��cross coefficient������������ϾͿ�����
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


%% �˻��������
Readfile=input('����������Ҫʹ���ļ����ļ���.csv .xlsx .xls �ļ����ɣ�����.csv .xlsx .xls:\n:','s');
Readedfile=[Readfile,'.csv'];
%LP=21;
try
datafile=csvread(Readedfile);
catch
datafile=xlsread([Readfile,'.xlsx']);
end
[NP,N]=size(datafile);
M=[Readedfile,'�ļ�����',num2str(NP),'�����ݵ㣬ÿ�����ݵ㹲',num2str(N),'��ά��'];
disp(M);
%Ԥ�����ڴ�
RHOcal=ones(NP,1);
RHO_temp=0;
%��5��ά�ȵ�������ȫ����ȡ������á�
Z=datafile(:,5);
T=datafile(:,1);
P=datafile(:,2);
rho=datafile(:,4);
X1=datafile(:,3);
X2=1-X1;
%����virial���̵��������
data=datafile(:,1:4);
%X1=data(:,3);
%X2=1-X1;
askoperation=input('what do you want to do? Optimize(Y)/caculate(C or other random symbol):\n��������ϻ��Ǽ��㣬�������Y����������C)\n ������:','s');
%% read coeffcient from file viralco.xls
virial_coeff=csvread('viralco_ver2.csv');
virial_coeff=virial_coeff';
disp('read done(viralcover1.csvϵ���ļ���ȡ���)');
%% ȫ���Ż�
if askoperation=='Y'
% ʹ����ϳ��򣬿�ʼ���
disp('��ʼ��ϣ��԰�����......');
%coeff0=coeff;
lb=[];
ub=[];
%L=lsqcurvefit(@viral3,L0,data,Zexp,lb,ub,options);
% �鿴���-Pcal
%Zcal=viral3(L,data)
%Dev=(Zcal-Zexp)./Zexp
% ��ѹ��ֵ�Ա�׼������С�������

model=@viral_TRHO2P;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',data,'ydata',P,'x0',virial_coeff,'lb',lb,'ub',ub);

%{
% ��ѹ������ֵΪ��׼������С�������
model=@viral_TRHO2Z;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',data,'ydata',Z,'x0',virial_coeff,'lb',lb,'ub',ub);
%���ģ���� model �ж��壬������data���壬���P,�����ĳ�ʼֵ��coeff0
%}
ms=MultiStart;
%coeff_fin�����յĴ�ȡֵ
[coeff_fin,~,~,~,~]=run(ms,problem,50);
virial_coeff=coeff_fin;
%disp('optimize done');
disp('������');
else
    %����ϣ���������һ�ε��Ż�
    disp('��������ϣ�ֱ��ʹ��ϵ������');
end
%% ������β����κθĶ����Ķ������ϲ�����ϲ��֣���ȡ28��ϵ������������ϵ�ֻ��28-12-4=12 ��ϵ��
%% caculate Pressure/kPa
disp('���㿪ʼ,�Ե�Ƭ��......')
P_cal=viral_TRHO2P(virial_coeff,data);
%Z_cal=P_cal./(rho.*R*T);
Dev_P=(P_cal-P)./P;
Dev_P=100*Dev_P;
%L=L1;
%% ά��ϵ��
% ��������η��̵ĽǶ���P X T�����RHO
b1=virial_coeff(1);
b2=virial_coeff(2);
b3=virial_coeff(3);
b4=virial_coeff(4);
b5=virial_coeff(5);
b6=virial_coeff(6);
b7=virial_coeff(7);
b8=virial_coeff(8);
b9=virial_coeff(9);
b10=virial_coeff(10);
b11=virial_coeff(11);
b12=virial_coeff(12);
b13=virial_coeff(13);
b14=virial_coeff(14);
b15=virial_coeff(15);
b16=virial_coeff(16);
b17=virial_coeff(17);
b18=virial_coeff(18);
c1=virial_coeff(19);
c2=virial_coeff(20);
c3=virial_coeff(21);
c4=virial_coeff(22);
c5=virial_coeff(23);
c6=virial_coeff(24);
c7=virial_coeff(25);
c8=virial_coeff(26);
c9=virial_coeff(27);
c10=virial_coeff(28);

NB=18;
NC=10;

Tc12=sqrt(Tc1*Tc2);
Tc112=(Tc1*Tc1*Tc2)^(1/3);
Tc221=(Tc1*Tc2*Tc2)^(1/3);
Tr1=T./Tc1;
Tr2=T./Tc2;
Tr12=T./Tc12;
Tr112=T./Tc112;
Tr221=T./Tc221;

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

% 2019-3-16-����
% Ѱ�Ҹ��Ĺ���Ӧ���ǵ���ϵģ���Ӧ�ô��ں�ʵ�����ݵ��κι�ϵ
% �˴�����������״̬�������Ԥ��ֵ����Ѱ�ҵ���������ʵ�ʸ�ʱ����λ����ͬʹ����������״̬����
% ������������С��ʵ��
% ��ˣ��ÿ�ܿ�����ϵͳ������ʹ��
rhoj=zeros(3,1);    %����һ���յ��м�������������洢���η��̵Ķ���ʵ�� 
deltaRHOj=zeros(3,1);

% ��ÿһ�����̿���ɸѡʵ�ʸ�
for i=1:NP
    % j����ͳ�����η���ʵ�����ĸ�����һ����1����Ҳ�п�����2��������ʱ��Ҫ�ҳ���Ҫ��ֵ
    j=0;
    % s Ϊһ����Ԫ����
    sroot=roots(tri_coeff(i,:));
    for k=1:length(sroot);
        real_S=isreal(sroot(k));
        if (real_S==1)&&(sroot(k)>0)
            j=j+1;
            rhoj(j)=sroot(k);%����
            deltaRHOj(j)=rhoj(j)-P(i)/(R*T(i));   
           %deltaRHOj(j)=rhoj(j)-rho(i);   
           % rhoc=s(k);
        else
            %pass
        end
    end
 % ��λ��ʵֵ����λ�ã�indexΪ���±�
[deltamin,index]=min(abs(deltaRHOj(1:j)));
%delta=deltaRHOj(index);
%rhoc=rho(i)+delta;
%RHO_temp�ݴ����ĸ�
RHO_temp(i)=rhoj(index);
end
RHOcal=RHO_temp';
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
figure(4);
plot(T,1000000*Cm,'o');
title('Cm*1000');
%% caculate p and density
%%  use NIST to caculate critical point
%disp('Caculate done now save coefficents')
disp('(�����������ȡά��ϵ��, ��ȡ������)');
disp('')
disp('')
%% save coefficents
csvwrite('viralco_ver1.csv',virial_coeff');
for i=1:NB
    disp(['b',num2str(i),' = ',num2str(virial_coeff(i))])
end
for i=1:NC
    disp(['c',num2str(i),' = ',num2str(virial_coeff(i+NB))])
end
%SaveYN=input('�Ƿ��ȡ����ѹ��ƫ��ܶ�ƫ������(Y|N):','s');
%if SaveYN=='Y'
DevFile=[Readfile,'_CalResults.csv'];
DevFiledata=table(T,P,X1,rho,RHOcal,P_cal,Dev_P,Dev_RHO,Bm,Cm);
writetable(DevFiledata,DevFile);
%end
end

function P=viral_TRHO2P_fit(virial_coeff,inputdata)
% �ĺ�������Ϊ�������
% ������ϵ�ֻ��12��ϵ��
% data: input data matrix
% L: coefficients
global Tc1 Tc2
%Tc1=375.31;% HFC-161
%Tc2=374.25; %HFC-134a
%Tc2=351.26 ;%HFC-32
R=8.3144;
T=inputdata(:,1);
RHO=inputdata(:,4);
X1=inputdata(:,3);
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

b1=virial_coeff(1);
b2=virial_coeff(2);
b3=virial_coeff(3);
b4=virial_coeff(4);
b5=virial_coeff(5);
b6=virial_coeff(6);
b7=virial_coeff(7);
b8=virial_coeff(8);
b9=virial_coeff(9);
b10=virial_coeff(10);
b11=virial_coeff(11);
b12=virial_coeff(12);
b13=virial_coeff(13);
b14=virial_coeff(14);
b15=virial_coeff(15);
b16=virial_coeff(16);
b17=virial_coeff(17);
b18=virial_coeff(18);
c1=virial_coeff(19);
c2=virial_coeff(20);
c3=virial_coeff(21);
c4=virial_coeff(22);
c5=virial_coeff(23);
c6=virial_coeff(24);
c7=virial_coeff(25);
c8=virial_coeff(26);
c9=virial_coeff(27);
c10=virial_coeff(28);

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
P=Z.*RHO.*T*R;
end

function P=viral_TRHO2P(virial_coeff,inputdata)
% data: input data matrix
% L: coefficients
global Tc1 Tc2
%Tc1=375.31;% HFC-161
%Tc2=374.25; %HFC-134a
%Tc2=351.26 ;%HFC-32
R=8.3144;
T=inputdata(:,1);
RHO=inputdata(:,4);
X1=inputdata(:,3);
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

b1=virial_coeff(1);
b2=virial_coeff(2);
b3=virial_coeff(3);
b4=virial_coeff(4);
b5=virial_coeff(5);
b6=virial_coeff(6);
b7=virial_coeff(7);
b8=virial_coeff(8);
b9=virial_coeff(9);
b10=virial_coeff(10);
b11=virial_coeff(11);
b12=virial_coeff(12);
b13=virial_coeff(13);
b14=virial_coeff(14);
b15=virial_coeff(15);
b16=virial_coeff(16);
b17=virial_coeff(17);
b18=virial_coeff(18);
c1=virial_coeff(19);
c2=virial_coeff(20);
c3=virial_coeff(21);
c4=virial_coeff(22);
c5=virial_coeff(23);
c6=virial_coeff(24);
c7=virial_coeff(25);
c8=virial_coeff(26);
c9=virial_coeff(27);
c10=virial_coeff(28);

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
P=Z.*RHO.*T*R;
end
