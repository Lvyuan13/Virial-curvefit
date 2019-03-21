% function curvefit_viral_ver2
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
% 
format long;
warning('off')
% 2019-03-20������һ���ǳ���������⣬��ʹ����ϳ���ʱ������Ľ��˿����ˬ�����ȼ��ѡ����ǵ���Ϻ��ȡϵ��֮��
% �ٴ�ʹ����Ϻõ�ϵ�����м���ʱ�����ȶ�ʱ���˵�0.8���ϡ������Ǵ�ȡ.xlsx .csv �ļ�������ʧ������
%�������峣��
R=8.3144;
% �����ٽ�ֵ
global Tc1 Tc2
% Tc1=304.1282;% CO2
% Tc2=369.825; %Propane
Tc1=369.825;%Propane
Tc2=396.44;%trifluoroiodomethane
NB=18;
NC=10;
% tolrance
tolr=1e-7;

%Tc1=375.31;% HFC-161
%Tc2=374.25; %HFC-134a
%Tc2=351.26 ;%HFC-32
%% ȫ����0 �����������


%% �˻��������
ReadFile=input('����������Ҫʹ���ļ����ļ���.csv .xlsx .xls �ļ����ɣ�����.csv .xlsx .xls:\n:','s');
ReadedFile=[ReadFile,'.csv'];
%LP=21;
try
DataFile=csvread(ReadedFile);
catch
DataFile=xlsread([ReadFile,'.xlsx']);
end
[NP,N]=size(DataFile);
M=[ReadedFile,'�ļ�����',num2str(NP),'�����ݵ㣬ÿ�����ݵ㹲',num2str(N),'��ά��'];
disp(M);

%��5��ά�ȵ�������ȫ����ȡ������á�

T=DataFile(:,1);
P=DataFile(:,2);
X1=DataFile(:,3);
X2=1-X1;
rho=DataFile(:,4);
%Z=datafile(:,5);

%����virial���̵��������
Data=DataFile(:,1:4);
%X1=data(:,3);
%X2=1-X1;
AskOpe=input('what do you want to do? Optimize(Y)/caculate(C or other random symbol):\n��������ϻ��Ǽ��㣬�������Y����������C)\n ������:','s');
%% read coeffcient from file viralco.xls
VirCoef=csvread('viralco_ver2.csv');
VirCoef=VirCoef';
VirCoefTheo=csvread('viralco_ver2_theo.csv');
VirCoefTheo=VirCoefTheo';
%disp('read done(viralcover1.csvϵ���ļ���ȡ���)');
for i=1:NB
    disp(['b',num2str(i),' = ',num2str(VirCoef(i))])
end
for i=1:NC
    disp(['c',num2str(i),' = ',num2str(VirCoef(i+NB))])
end
%% ȫ���Ż�
if AskOpe=='Y'
% ʹ����ϳ��򣬿�ʼ���
disp('��ʼ��ϣ��԰�����......');
%coeff0=coeff;
lb=VirCoefTheo-tolr;
ub=VirCoefTheo+tolr;
for i=13:18
    lb(i)=-Inf;  
    ub(i)=Inf;
end
for i=23:28
    lb(i)=-Inf;  
    ub(i)=Inf;
end
%L=lsqcurvefit(@viral3,L0,data,Zexp,lb,ub,options);
% �鿴���-Pcal
%Zcal=viral3(L,data)
%Dev=(Zcal-Zexp)./Zexp
% ��ѹ��ֵ�Ա�׼������С�������
% ���Ż������趨�ĵ�ʽԼ��
% b=VirCoefTheo
% A=diag(ones(28,1));
% for i=13:18
%     A(i,i)=0;  
% end
% for i=23:28
%     A(i,i)=0;
% end


model=@viral_TRHO2P_fit;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',Data,'ydata',P,'x0',VirCoef,'lb',lb,'ub',ub);

%{
% ��ѹ������ֵΪ��׼������С�������
model=@viral_TRHO2Z;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',data,'ydata',Z,'x0',virial_coeff,'lb',lb,'ub',ub);
%���ģ���� model �ж��壬������data���壬���P,�����ĳ�ʼֵ��coeff0
%}
ms=MultiStart;
%coeff_fin�����յĴ�ȡֵ

[coeff_fin,~,~,~,~]=run(ms,problem,5);
VirCoef=coeff_fin;
%disp('optimize done');

disp('������,�����Ѿ����浽.mat�ļ�');
else
    %����ϣ���������һ�ε��Ż�
    disp('��������ϣ�ֱ��ʹ��ϵ������');
end
%% ������β����κθĶ����Ķ������ϲ�����ϲ��֣���ȡ28��ϵ������������ϵ�ֻ��28-12-4=12 ��ϵ��
%% caculate Pressure/kPa
disp('���㿪ʼ,�Ե�Ƭ��......')
P_cal=viral_TRHO2P(VirCoef,Data);
%%
RHOcal=viral_TP2RHO(VirCoef,Data);
%%
Bm=VIRB(VirCoef,Data);
Cm=VIRC(VirCoef,Data);
%Z_cal=P_cal./(rho.*R*T);
Dev_P=(P_cal-P)./P;
Dev_P=100*Dev_P;
%L=L1;

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
DevFile=[ReadFile,'_CalResults.csv'];
DevFiledata=table(T,P,X1,rho,RHOcal,P_cal,Dev_P,Dev_RHO,1000*Bm,1000000*Cm);
writetable(DevFiledata,DevFile);
disp('')
% disp('������Ϻ�Ľ��Ϊ��')
%% save coefficents
csvwrite('viralco_ver2.csv',VirCoef');
for i=1:NB
    disp(['b',num2str(i),' = ',num2str(VirCoef(i))])
end
for i=1:NC
    disp(['c',num2str(i),' = ',num2str(VirCoef(i+NB))])
end
%SaveYN=input('�Ƿ��ȡ����ѹ��ƫ��ܶ�ƫ������(Y|N):','s');
%if SaveYN=='Y'

%end
% end
