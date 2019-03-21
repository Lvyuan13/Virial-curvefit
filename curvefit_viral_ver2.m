% function curvefit_viral_ver2
%% 局部优化拟合函数
%NP=140; Number of Points
%设置为长数据型；
% version 2
% 2019-03-20 修正
% 重新拟合了维里系数，不再采用陈光明教授课题组的多项式拟合，而是根据基本原理
% 设置了新的线性多项式方法，具体公式见论文，已经改正；系数由21个调整成28个
% 此外最重要的一点，不再用纯粹数据拟合的方法拟合维里系数，经过周密验证，以前使用的
% 加入部分纯物质数据拟合的方法仍然具有相当的合理性，在系统仿真中是完全可以应用的，也具有一定的
% 外推性，但是，经过核算发现，从前的方法，计算二阶维里系数Bm时精度甚高，但是在计算三阶维里系数
% 时发现，计算出的Cm和理论推导出来的差距不小,甚至在小范围内会出现Cm计算结果为负值的情况，这无疑是不正确的
% 现在采用的新方法：
% 1. propane的维里系数Bm Cm 使用refprop 软件中使用 Helmholtz模型的计算值，该理论在论文中已经补上，再根据Bm Cm 
% 计算值将Bm和 Cm 拟合成Tr 的多项式， 这无疑时更加合理的做法
% 2. CF3I 的维里系数，采用段远源教授的方程计算出，也是先根据计算出的Bm Cm，再将Bm和Cm 拟合成Tr的多项式，
% 通过核算发现，该法拟合出的小系数 B1 B2 B3 etc 与其论文中结果完全一致，因此可以判定其结果时完全正确的
% 3. 第三步，将上述两个步骤的纯物质系数直接嵌入方程，这样需要拟合的数据只有cross coefficient，加以数据拟合就可以了
% 
format long;
warning('off')
% 2019-03-20：存在一个非常不解的问题，当使用拟合程序时计算出的结果丝毫不爽，精度极佳。但是当拟合后存取系数之后
% 再次使用拟合好的系数进行计算时，精度顿时会退到0.8以上。怀疑是存取.xlsx .csv 文件精度损失的问题
%理想气体常数
R=8.3144;
% 两个临界值
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
%% 全部归0 方便二次运行


%% 人机交互设计
ReadFile=input('请输入你需要使用文件的文件名.csv .xlsx .xls 文件均可，不加.csv .xlsx .xls:\n:','s');
ReadedFile=[ReadFile,'.csv'];
%LP=21;
try
DataFile=csvread(ReadedFile);
catch
DataFile=xlsread([ReadFile,'.xlsx']);
end
[NP,N]=size(DataFile);
M=[ReadedFile,'文件共有',num2str(NP),'个数据点，每个数据点共',num2str(N),'个维度'];
disp(M);

%将5个维度的数据先全部读取并定义好。

T=DataFile(:,1);
P=DataFile(:,2);
X1=DataFile(:,3);
X2=1-X1;
rho=DataFile(:,4);
%Z=datafile(:,5);

%定义virial方程的输入矩阵
Data=DataFile(:,1:4);
%X1=data(:,3);
%X2=1-X1;
AskOpe=input('what do you want to do? Optimize(Y)/caculate(C or other random symbol):\n（你想拟合还是计算，拟合输入Y，计算输入C)\n 请输入:','s');
%% read coeffcient from file viralco.xls
VirCoef=csvread('viralco_ver2.csv');
VirCoef=VirCoef';
VirCoefTheo=csvread('viralco_ver2_theo.csv');
VirCoefTheo=VirCoefTheo';
%disp('read done(viralcover1.csv系数文件读取完毕)');
for i=1:NB
    disp(['b',num2str(i),' = ',num2str(VirCoef(i))])
end
for i=1:NC
    disp(['c',num2str(i),' = ',num2str(VirCoef(i+NB))])
end
%% 全局优化
if AskOpe=='Y'
% 使用拟合程序，开始拟合
disp('开始拟合，稍安勿躁......');
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
% 查看结果-Pcal
%Zcal=viral3(L,data)
%Dev=(Zcal-Zexp)./Zexp
% 以压力值以标准进行最小二乘拟合
% 给优化问题设定的等式约束
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
% 以压缩因子值为标准进行最小二乘拟合
model=@viral_TRHO2Z;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',data,'ydata',Z,'x0',virial_coeff,'lb',lb,'ub',ub);
%拟合模型是 model 中定义，输入是data定义，输出P,给定的初始值是coeff0
%}
ms=MultiStart;
%coeff_fin是最终的存取值

[coeff_fin,~,~,~,~]=run(ms,problem,5);
VirCoef=coeff_fin;
%disp('optimize done');

disp('拟合完毕,变量已经保存到.mat文件');
else
    %不拟合，跳过上以一段的优化
    disp('不进行拟合，直接使用系数计算');
end
%% 计算这段不加任何改动，改动的是上部的拟合部分，读取28个系数，而参与拟合的只有28-12-4=12 个系数
%% caculate Pressure/kPa
disp('计算开始,稍等片刻......')
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
disp('计算无误，开始出图......');
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
disp('(程序结束，存取维里系数, 存取计算结果)');
DevFile=[ReadFile,'_CalResults.csv'];
DevFiledata=table(T,P,X1,rho,RHOcal,P_cal,Dev_P,Dev_RHO,1000*Bm,1000000*Cm);
writetable(DevFiledata,DevFile);
disp('')
% disp('经过拟合后的结果为：')
%% save coefficents
csvwrite('viralco_ver2.csv',VirCoef');
for i=1:NB
    disp(['b',num2str(i),' = ',num2str(VirCoef(i))])
end
for i=1:NC
    disp(['c',num2str(i),' = ',num2str(VirCoef(i+NB))])
end
%SaveYN=input('是否存取计算压力偏差，密度偏差数据(Y|N):','s');
%if SaveYN=='Y'

%end
% end
