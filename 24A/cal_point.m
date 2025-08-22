function [point_x, point_y] = cal_point(b, theta0)

R = 4.5;
theta = R/b - theta0;

point_x = b*theta*cos(theta);
point_y = b*theta*sin(theta);

end

