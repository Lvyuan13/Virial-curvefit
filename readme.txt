1.运行globalrun.m
2.程序会询问使用什么文件（which file do you want to use:）
填你想用的数据文件名，如1112.xlsx文件，填1112。enter,程序会将1112.xlsx中的数据加载到矩阵中
3.程序再次询问，你想用这个数据干什么？（what do you want to do? Optimize(Y)/caculate(C or other random symbol):）
       填Y，程序会拟合出数据，并且将拟合好的系数放到文件viralco.xls文件中。自动生成计算效果图
       不填Y，填C或者其他任何字母，程序会直接用已经存放在viralco.xls中的系数计算该数据。自动生成计算效果图
4.如果想用A数据拟合，用A数据拟合的结果去处理B数据计算，那么先按1到3步骤运行，3步骤填Y，然后再次运行程序，
安照1-3步骤运行，三步骤填C。此时出的图计算用A数据拟合后，再用A拟合的系数去计算B的效果了