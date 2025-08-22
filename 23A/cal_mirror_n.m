function [every_mirr_n_unit_vector_array, reflect_light_unit_vector_array]= cal_mirror_n(S,x,y,z)
% 输入：
%   S: 1×3 单位向量，太阳光方向向量，某一个时刻下
%   x, y, z: N×1 向量，镜子中心坐标
% 输出：
%   every_mirr_n_unit_vector_array: N×3，镜子法向量
%   reflect_light_unit_vector_array: N×3，反射光单位向量
%   mirror_cneter: N×3，镜子中心坐标
S = S(:).';

x = x(:); 
y = y(:); 
z = z(:); 

N = length(x);

% 构建 mirror_cneter 为 N×3 数组
mirror_cneter = [x, y, z]; % N×3

% 将 T_center 复制为 N×3
T_center = repmat([0, 0, 84], N,1); % N×3



% 反射光直线方程
reflect_light_vector_array = mirror_cneter - T_center;

% 计算每个反射向量的模
magnitudes = sqrt(sum(reflect_light_vector_array.^2, 2));

if any(magnitudes < 1e-6)
    error('Mirror(s) at T_center ([0, 0, 84]) cause zero magnitude');
end

reflect_light_unit_vector_array = reflect_light_vector_array./magnitudes;  % 反射光向量的单位向量

every_mirr_n_unit_vector_array = (reflect_light_unit_vector_array - S);

% 归一化法向量
norms = sqrt(sum(every_mirr_n_unit_vector_array.^2, 2)); % 1×N, 计算模
every_mirr_n_unit_vector_array = every_mirr_n_unit_vector_array ./ norms;

end