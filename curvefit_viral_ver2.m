function curvefit_viral_ver2
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
format long;
%理想气体常数
R=8.3144;
% 两个临界值
global Tc1 Tc2
% Tc1=304.1282;% CO2
% Tc2=369.825; %Propane
Tc1=369.825;%Propane
Tc2=396.44;%trifluoroiodomethane
%Tc1=375.31;% HFC-161
%Tc2=374.25; %HFC-134a
%Tc2=351.26 ;%HFC-32
%% 全部归0 方便二次运行


%% 人机交互设计
Readfile=input('请输入你需要使用文件的文件名.csv .xlsx .xls 文件均可，不加.csv .xlsx .xls:\n:','s');
Readedfile=[Readfile,'.csv'];
%LP=21;
try
datafile=csvread(Readedfile);
catch
datafile=xlsread([Readfile,'.xlsx']);
end
[NP,N]=size(datafile);
M=[Readedfile,'文件共有',num2str(NP),'个数据点，每个数据点共',num2str(N),'个维度'];
disp(M);
%预分配内存
RHOcal=ones(NP,1);
RHO_temp=0;
%将5个维度的数据先全部读取并定义好。
Z=datafile(:,5);
T=datafile(:,1);
P=datafile(:,2);
rho=datafile(:,4);
X1=datafile(:,3);
X2=1-X1;
%定义virial方程的输入矩阵
data=datafile(:,1:4);
%X1=data(:,3);
%X2=1-X1;
askoperation=input('what do you want to do? Optimize(Y)/caculate(C or other random symbol):\n（你想拟合还是计算，拟合输入Y，计算输入C)\n 请输入:','s');
%% read coeffcient from file viralco.xls
virial_coeff=csvread('viralco_ver2.csv');
virial_coeff=virial_coeff';
disp('read done(viralcover1.csv系数文件读取完毕)');
%% 全局优化
if askoperation=='Y'
% 使用拟合程序，开始拟合
disp('开始拟合，稍安勿躁......');
%coeff0=coeff;
lb=[];
ub=[];
%L=lsqcurvefit(@viral3,L0,data,Zexp,lb,ub,options);
% 查看结果-Pcal
%Zcal=viral3(L,data)
%Dev=(Zcal-Zexp)./Zexp
% 以压力值以标准进行最小二乘拟合

model=@viral_TRHO2P;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',data,'ydata',P,'x0',virial_coeff,'lb',lb,'ub',ub);

%{
% 以压缩因子值为标准进行最小二乘拟合
model=@viral_TRHO2Z;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',data,'ydata',Z,'x0',virial_coeff,'lb',lb,'ub',ub);
%拟合模型是 model 中定义，输入是data定义，输出P,给定的初始值是coeff0
%}
ms=MultiStart;
%coeff_fin是最终的存取值
[coeff_fin,~,~,~,~]=run(ms,problem,50);
virial_coeff=coeff_fin;
%disp('optimize done');
disp('拟合完毕');
else
    %不拟合，跳过上以一段的优化
    disp('不进行拟合，直接使用系数计算');
end
%% 计算这段不加任何改动，改动的是上部的拟合部分，读取28个系数，而参与拟合的只有28-12-4=12 个系数
%% caculate Pressure/kPa
disp('计算开始,稍等片刻......')
P_cal=viral_TRHO2P(virial_coeff,data);
%Z_cal=P_cal./(rho.*R*T);
Dev_P=(P_cal-P)./P;
Dev_P=100*Dev_P;
%L=L1;
%% 维里系数
% 下面从三次方程的角度由P X T反求出RHO
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
% 下设置三次方程的系数矩阵
I=ones(1,NP); 
I=I';
tri_coeff=[Cm,Bm,I,-PRT];
% 2018-11-07-修正
% 出错原因分析：
% 在少许情况，三次方程存在两个或两个以上的实际根，此时以遍历第一次遇到的
% 实跟为正确值是非合理的
% 
% 处理方法要求：
%     1）求解中对初始的实验值不能有依赖
%     2）不可过分提升程序复杂度
%     3）解决思路如下：
%          在多个实跟中必然有一个值是最接近实际值的，但是不能以实验值作为比较对象而筛选，
%          这样程序会失去通用性，因而需要先以PR方程的计算值为预测值，这个值是非常接近实
%          际的而后再使用选择viral方程求解的多个实根中最接近PR方程预测值的量。
%     4) 采用PR计算需要重新编写PR的程序，不经济，直接用取多个实际根中接近实验值的

% 2019-3-16-修正
% 寻找根的过程应该是低耦合的，不应该存在和实验数据的任何关系
% 此处用理想气体状态方程替代预测值，当寻找到两个以上实际根时，定位根到同使用理想气体状态方程
% 计算结果差异最小的实根
% 如此，该框架可以在系统仿真中使用
rhoj=zeros(3,1);    %定义一个空的中间变量数组用来存储三次方程的多重实根 
deltaRHOj=zeros(3,1);

% 对每一个方程开启筛选实际根
for i=1:NP
    % j用以统计三次方程实数根的个数，一般是1个，也有可能是2个，两个时需要找出需要的值
    j=0;
    % s 为一个三元向量
    sroot=roots(tri_coeff(i,:));
    for k=1:length(sroot);
        real_S=isreal(sroot(k));
        if (real_S==1)&&(sroot(k)>0)
            j=j+1;
            rhoj(j)=sroot(k);%测试
            deltaRHOj(j)=rhoj(j)-P(i)/(R*T(i));   
           %deltaRHOj(j)=rhoj(j)-rho(i);   
           % rhoc=s(k);
        else
            %pass
        end
    end
 % 定位真实值所在位置，index为其下标
[deltamin,index]=min(abs(deltaRHOj(1:j)));
%delta=deltaRHOj(index);
%rhoc=rho(i)+delta;
%RHO_temp暂存合理的根
RHO_temp(i)=rhoj(index);
end
RHOcal=RHO_temp';
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
%SaveYN=input('是否存取计算压力偏差，密度偏差数据(Y|N):','s');
%if SaveYN=='Y'
DevFile=[Readfile,'_CalResults.csv'];
DevFiledata=table(T,P,X1,rho,RHOcal,P_cal,Dev_P,Dev_RHO,Bm,Cm);
writetable(DevFiledata,DevFile);
%end
end

function P=viral_TRHO2P_fit(virial_coeff,inputdata)
% 改函数纯粹为了拟合用
% 参与拟合的只有12个系数
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
