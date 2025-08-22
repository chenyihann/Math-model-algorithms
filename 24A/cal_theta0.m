function [theta0_initial, theta0_values] = cal_theta0(t_range, start_theta0, b, C, options)
% 初始化存储数组theta0
theta0_values = zeros(1, length(t_range));
L = 2.86;
for i = 1:length(t_range)
    
    if i==1
        theta0_val = start_theta0;
    else
        theta0_prev = theta0_values(i-1);  % 防止初值跳变
        r_prev = b * theta0_prev;
        theta0_guess = theta0_prev - L/r_prev;
        t = t_range(i);
        fun1 = @(theta) b*0.5*(theta*sqrt(theta^2 + 1) + log(theta + sqrt(theta^2 + 1))) + t + C;
        theta0_val = fsolve(fun1,theta0_guess,options);
    end
    theta0_values(i) = theta0_val;
end
theta0_initial = theta0_values(length(t_range));  % 获取最后一秒时的数值theta0
end