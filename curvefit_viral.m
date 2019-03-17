function curvefit_viral
%% 局部优化拟合函数
%NP=140; Number of Points
%设置为长数据型；
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
Zexp=datafile(:,5);
T=datafile(:,1);
P=datafile(:,2);
rho=datafile(:,4);
X1=datafile(:,3);
X2=1-X1;
%定义virial方程的输入矩阵
data=datafile(:,1:4);
%X1=data(:,3);
%X2=1-X1;
askoperation=input('what do you want to do? Optimize(Y)/caculate(C or other random symbol):\n（你想拟合还是计算，拟合输入Y，计算输入C)\n','s');
%% read coeffcient from file viralco.xls
virial_coeff=csvread('viralco.csv');
virial_coeff=virial_coeff';
disp('read done(viralco.csv系数文件读取完毕)');
%% 全局优化
if askoperation=='Y'
% 使用拟合程序，开始拟合
disp('开始拟合，稍安勿躁......');
%coeff0=coeff;
options=optimset('TolFun',1e-10);
lb=[];
ub=[];
%L=lsqcurvefit(@viral3,L0,data,Zexp,lb,ub,options);
% 查看结果-Pcal
%Zcal=viral3(L,data)
%Dev=(Zcal-Zexp)./Zexp
model=@viral_TRHO2P;
problem=createOptimProblem('lsqcurvefit','objective', ...
    model,'xdata',data,'ydata',P,'x0',virial_coeff,'lb',lb,'ub',ub);
%拟合模型是 model 中定义，输入是data定义，输出P,给定的初始值是coeff0
ms=MultiStart;
%coeff_fin是最终的存取值
[coeff_fin,fval,exitflag,output,solutions]=run(ms,problem,50);
virial_coeff=coeff_fin;
%disp('optimize done');
disp('拟合完毕');
else
    %不拟合，跳过上以一段的优化
    disp('不进行拟合，直接使用系数计算');
end
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
%% save coefficents
csvwrite('viralco.csv',virial_coeff');
%SaveYN=input('是否存取计算压力偏差，密度偏差数据(Y|N):','s');
%if SaveYN=='Y'
DevFile=[Readfile,'_CalResults.csv'];
DevFiledata=table(T,P,X1,rho,RHOcal,P_cal,Dev_P,Dev_RHO,Bm,Cm);
writetable(DevFiledata,DevFile);
%end
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