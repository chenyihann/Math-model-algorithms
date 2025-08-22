%% 问题1：由递推公式来求得每秒的速度和位置
%% 初始化参数
clc,clear;

total_num = 224;  % 把手数目
options = optimoptions('fsolve', 'Display','off', 'FunctionTolerance',1e-9);

tic;
%% 求解theta与t的关系式

b = 0.55/(2*pi);
t_val = 0:1:300;  % 总时长300s
C = -b*0.5*(32*pi*sqrt((32*pi)^2 + 1) + log(32*pi + sqrt((32*pi)^2 + 1)));

theta0_values = zeros(1, length(t_val));

for i = 1:length(t_val)
    if i==1
        theta0_val = 32*pi;
    else
        theta0_prev = theta0_values(i-1) - 1;
        t = t_val(i);
        fun1 = @(theta) b*0.5*(theta*sqrt(theta^2 + 1) + log(theta + sqrt(theta^2 + 1))) + t + C;
        theta0_val = fsolve(fun1,theta0_prev,options);
    end
    theta0_values(i) = theta0_val;
end

theta_sol = zeros(length(t_val), total_num); % 初始化数组，存储所有解

%% 求解各个把手theta的递推关系式

for i = 1:length(t_val)
    for j = 1:total_num
        if j == 1
            theta_val = theta0_values(i);
        else
            if j == 2
                L = 2.86;
            else
                L = 1.65;
            end
            
            theta_prev = theta_sol(i, j-1);  % 使用同一个时间，前一个把手的解作为初值
            
            fun2 = @(theta_k) b^2*(theta_prev^2 + theta_k^2) - L^2 - 2*b^2*theta_k*theta_prev*cos(theta_k - theta_prev);
            theta_val = fsolve(fun2, theta_prev, options);
        end 
        theta_sol(i, j) = theta_val; % 索引从1开始
    end
end

%% 求解所有把手在0-300s时的位置
x_pos = b*theta_sol.*cos(theta_sol);
y_pos = b*theta_sol.*sin(theta_sol);

draw_full_dragon(x_pos(end,:),y_pos(end,:), b);

%% 速度的求解
v_sol = zeros(length(t_val), total_num); % 初始化数组，存储所有解
for i = 1:length(t_val)
    for j = 1:total_num
        if j==1
            v_val = 1;
            v_sol(i,j) = v_val;
        else
            theta1 = theta_sol(i,j-1);
            theta2 = theta_sol(i,j);
            k1 = (b*sin(theta1) + b*theta1*cos(theta1))/(b*cos(theta1)-b*theta1*sin(theta1));
            k2 = (b*sin(theta2) + b*theta2*cos(theta2))/(b*cos(theta2)-b*theta2*sin(theta2));
            k = (y_pos(i,j)-y_pos(i,j-1))/(x_pos(i,j)-x_pos(i,j-1));
            apha = atan(abs((k1-k) / (1 + k1*k)));
            beta = atan(abs((k2-k) / (1 + k2*k)));  
            v_prev = v_sol(i,j-1);
            fun3 = @(v) v*cos(beta) - v_prev*cos(apha);
            v_val = fsolve(fun3,v_prev, options); 
            v_sol(i,j) = v_val;
        end
    end
end

toc;           

