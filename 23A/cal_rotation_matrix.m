function R = cal_rotation_matrix(n)
% 输出：单个定日镜的旋转矩阵，从地面坐标系到定日镜
% 输入：单个定日镜的法向量

cos_alpha = n(3);
alpha = acos(cos_alpha);
tan_gamma = n(1)/n(2);

gamma = atan(tan_gamma);

R_z = [cos(gamma), - sin(gamma), 0; sin(gamma), cos(gamma), 0; 0, 0, 1];
R_x = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];

R = R_x*R_z;  % 矩阵相乘


end