close all; clear
save_dir = 'E:\Research\Work\HCS_multifold_flapping_by_PSP_SDO\Encounter 7\';
e_LMN_file = [save_dir, 'e_LMN.csv'];
e_LMN_mean_file = [save_dir, 'e_LMN_mean.csv'];

e_LMN_arr = readmatrix(e_LMN_file);
e_L_arr = e_LMN_arr(:, 1:5);
e_M_arr = e_LMN_arr(:, 6:10);
e_N_arr = e_LMN_arr(:, 11:15);

e_LMN_mean = readmatrix(e_LMN_mean_file);
e_L_mean = e_LMN_mean(:,1);
e_M_mean = e_LMN_mean(:,2);
e_N_mean = e_LMN_mean(:,3);

n_M_arr = e_N_arr' * e_M_mean;
n_M_angle = acosd(n_M_arr);
n_N_arr = e_N_arr' * e_N_mean;
n_N_angle = acosd(n_N_arr);
delta_BL_arr = [-1;1;-1;1;-1];

k_arr = sign(n_M_arr .* n_N_arr .* delta_BL_arr);