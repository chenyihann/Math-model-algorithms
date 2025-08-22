%{


%}


clc,clear;
%% 初始化参数
b = 1.7/(2*pi);


%% 计算一系列点坐标
% 计算圆与螺旋线相交的店A和B
[x_a,y_a] = cal_point(b,0);
[x_b,y_b] = cal_point(b,pi);

% 计算两圆弧相交的点C
x_c = (1/3)* (2*x_b+x_a);
y_c = (1/3)* (2*y_b+y_a);

