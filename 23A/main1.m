clear,clc;
tic;
month_array = 1:1:12;

month_eta_sb_array = zeros(12,1);
month_eta_sb_mirr_reciever_array = zeros(12,1);
month_eta_sb_mirr_reflect_array = zeros(12,1);
month_eta_sb_tow2mirr_array = zeros(12,1);
for i = 1:12
    month = month_array(i);
    [eta_sb_array, eta_sb_mirr_reciever_array, eta_sb_mirr_reflect_array, eta_sb_tow2mirr_array] = cal_eta_sb(month);  
    month_avg_eta_sb = mean(eta_sb_array(:));
    month_avg_eta_sb_mirr_reciever = mean(eta_sb_mirr_reciever_array(:));
    month_avg_eta_sb_mirr_reflect = mean(eta_sb_mirr_reflect_array(:));
    month_avg_eta_sb_tow2mirr = mean(eta_sb_tow2mirr_array(:));
    % 储存每个月的数据
    month_eta_sb_array(i) = month_avg_eta_sb;
    month_eta_sb_mirr_reciever_array(i) = month_avg_eta_sb_mirr_reciever;
    month_eta_sb_mirr_reflect_array(i) = month_avg_eta_sb_mirr_reflect;
    month_eta_sb_tow2mirr_array(i) = month_avg_eta_sb_tow2mirr;
end

eta_sb_mirr_reciever = mean(month_avg_eta_sb_mirr_reciever(:));
eta_sb_mirr_reflect = mean(month_avg_eta_sb_mirr_reflect(:));
eta_sb_tow2mirr = mean(month_avg_eta_sb_tow2mirr(:));
eta_sb = mean(month_eta_sb_array(:));

fprintf('eta_sb_mirr_reciever:%f\n', eta_sb_mirr_reciever);
fprintf('eta_sb_mirr_reflect:%f\n', eta_sb_mirr_reflect);
fprintf('eta_sb_tow2mirr:%f\n', eta_sb_tow2mirr);
fprintf('eta_sb:%f\n', eta_sb);

toc;