function [sin_alpha_s, cos_alpha_s, sin_beta_s,cos_beta_s, DNI] = cal_alpha_and_beta(time, month)
% 非数组
% 初始化参数
varphi = 39.4*pi/180;  % 纬度(弧度)
omega = (time - 12)*(pi/12);   % 为太阳时角
D = [-59, -28, 0, 31, 61, 92, 122, 153, 184, 214, 254, 275];
sin_delta =  sin(2*pi*D(month)/365)*sin(deg2rad(23.45));  % 太阳赤纬角
delta = asin(sin_delta);
cos_delta = cos(delta);

% 计算高度角和方位角

sin_alpha_s = cos_delta*cos(varphi)*cos(omega) + sin_delta*sin(varphi);
alpha_s = asin(sin_alpha_s);
cos_alpha_s = cos(alpha_s);
cos_beta_s = (sin_delta - sin_alpha_s*sin(varphi))/(cos_alpha_s*cos(varphi));
beta_s = acos(cos_beta_s);
sin_beta_s = sin(beta_s);

% 计算DNI
G0 = 1.366;
H = 3;  % 单位千米
a = 0.4237-0.00821*(6-H)^2;
b = 0.5055+0.00595*(6.5-H)^2;
c = 0.2711+0.01858*(2.5-H)^2;

DNI = G0*(a+b*exp(-c/sin_alpha_s));


end

