
function RHOcal=viral_TP2RHO(VirCoef,InputData)
% data: input data matrix
% L: coefficients
format long;
global Tc1 Tc2
%Tc1=375.31;% HFC-161
%Tc2=374.25; %HFC-134a
%Tc2=351.26 ;%HFC-32
R=8.3144;
T=InputData(:,1);
P=InputData(:,2);
X1=InputData(:,3);
X2=1-X1;
%预分配内存
[NP,~]=size(InputData);
% RHOcal=ones(NP,1);
RHO_temp=0;
%RHO=inputdata(:,4);

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
%% 维里系数
% 下面从三次方程的角度由P X T反求出RHO
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

end