clear; close all;
encounter = 7;
dTime = 0.02;
save_or_not = 0;
mean_or_local = 0;
which_cross = 0;
data_dir = ['E:\Research\Data\PSP\Encounter ',num2str(encounter),'\'];
save_dir  = ['E:\Research\Work\HCS_multifold_flapping_by_PSP_SDO\Encounter ',num2str(encounter),'\'];

%% which encounter
year = 2021;
cross_num = 5;
spc_or_spi = 'i';
fld_lst = '12';

%% begin time and end time
plot_beg_time = '2021-01-17 12:00:00';
plot_end_time = '2021-01-17 16:00:00';
% calc_beg_time = '2021-01-17 12:00:00';
% calc_end_time = '2021-01-17 14:00:00';
calc_beg_time = '2021-01-17 13:00:00';
calc_end_time = '2021-01-17 15:00:00';
plot_beg_epoch = datenum(plot_beg_time);
plot_end_epoch = datenum(plot_end_time);

%% fld_data: Brtn
fld_file = ['psp_fld_l2_mag_rtn_',plot_beg_time(1:4),plot_beg_time(6:7),plot_beg_time(9:10),fld_lst,'_v02.cdf'];
fld_dir = [data_dir,fld_file];
fld_info = spdfcdfinfo(fld_dir);

fld_Epoch = spdfcdfread(fld_dir,'Variables','epoch_mag_RTN');
Brtn = spdfcdfread(fld_dir,'Variables','psp_fld_l2_mag_RTN');
Brtn(abs(Brtn)>1e3) = nan;
Br = Brtn(:,1); 
Bt = Brtn(:,2); 
Bn = Brtn(:,3); 
B_mod = sqrt(Br.^2 + Bt.^2 + Bn.^2);

fld_plot_index = find(fld_Epoch >= plot_beg_epoch & fld_Epoch <= plot_end_epoch);
fld_Epoch_plot = fld_Epoch(fld_plot_index);
Br_plot = Br(fld_plot_index); 
Bt_plot = Bt(fld_plot_index); 
Bn_plot = Bn(fld_plot_index); 
B_mod_plot = B_mod(fld_plot_index);

%% smooth the time series of magnetic field
window = 1000001;
Br_plot_sm = smooth1d(Br_plot, 'movmean', window);
Bt_plot_sm = smooth1d(Bt_plot, 'movmean', window);
Bn_plot_sm = smooth1d(Bn_plot, 'movmean', window);
B_mod_plot_sm = smooth1d(B_mod_plot, 'movmean', window);

%% show the field before and after smoothing
figure
LineWidth = 1.5;
FontSize = 16;

subplot(3,1,1)
plot(fld_Epoch_plot, Br_plot, 'r', 'LineWidth', LineWidth); hold on
plot(fld_Epoch_plot, Bt_plot, 'g', 'LineWidth', LineWidth); hold on
plot(fld_Epoch_plot, Bn_plot, 'b', 'LineWidth', LineWidth); hold on
plot(fld_Epoch_plot, B_mod_plot, 'k', 'LineWidth', LineWidth); 
legend('Br','Bt','Bn','|B|','Location','northeast')
datetick('x','HH:MM');
set(gca, 'FontSize', FontSize)

subplot(3,1,2)
plot(fld_Epoch_plot, Br_plot_sm, 'r', 'LineWidth', LineWidth); hold on
plot(fld_Epoch_plot, Bt_plot_sm, 'g', 'LineWidth', LineWidth); hold on
plot(fld_Epoch_plot, Bn_plot_sm, 'b', 'LineWidth', LineWidth); hold on
plot(fld_Epoch_plot, B_mod_plot_sm, 'k', 'LineWidth', LineWidth); hold on
legend('Br smooth','Bt smooth','Bn smooth','|B| smooth','Location','northeast')
datetick('x','HH:MM');
set(gca, 'FontSize', FontSize)

%% determine mean LMN coordinate    
calc_beg_epoch = datenum(calc_beg_time);
calc_end_epoch = datenum(calc_end_time);
[V_mean, D_mean,ratio21_mean,ratio32_mean] = get_LMN(Br,Bt,Bn,fld_Epoch,calc_beg_epoch,calc_end_epoch); % ratio21 >= 3 is required
e_L_mean = V_mean(:,3)/(sqrt(dot(V_mean(:,3), V_mean(:,3))));
e_M_mean = V_mean(:,2)/(sqrt(dot(V_mean(:,2), V_mean(:,2))));
e_N_mean = V_mean(:,1)/(sqrt(dot(V_mean(:,1), V_mean(:,1))));

e_LMN_mean = cat(2,e_L_mean,e_M_mean,e_N_mean);
% csvwrite([save_dir, 'e_LMN_mean.csv'], e_LMN_mean); % e_L, e_M, e_N

%% interpolate
nTime = floor((plot_end_epoch - plot_beg_epoch)*86400/dTime);
std_time = linspace(plot_beg_epoch, plot_end_epoch, nTime);
std_Epoch = interp1(fld_Epoch_plot,fld_Epoch_plot,std_time,'pchip');

Br_interp = interp1(fld_Epoch_plot,Br_plot_sm,std_time,'linear'); 
Br_interp(isnan(Br_interp)) = mean(Br_interp,'omitnan');
Bt_interp = interp1(fld_Epoch_plot,Bt_plot_sm,std_time,'linear'); 
Bt_interp(isnan(Bt_interp)) = mean(Bt_interp,'omitnan');
Bn_interp = interp1(fld_Epoch_plot,Bn_plot_sm,std_time,'linear'); 
Bn_interp(isnan(Bn_interp)) = mean(Bn_interp,'omitnan');
Brtn_interp = cat(2,Br_interp',Bt_interp',Bn_interp');
B_mod_interp = sqrt(Br_interp.^2 + Bt_interp.^2 + Bn_interp.^2);

%% calculate B_LMN
Bnml_interp = Brtn_interp * V_mean; 
BL_interp = Bnml_interp(:,3);
BM_interp = Bnml_interp(:,2);
BN_interp = Bnml_interp(:,1);

subplot(3,1,3)
plot(std_Epoch, BL_interp, 'r', 'LineWidth', LineWidth); hold on
plot(std_Epoch, BM_interp, 'g', 'LineWidth', LineWidth); hold on
plot(std_Epoch, Bn_interp, 'b', 'LineWidth', LineWidth); hold on
legend('B_L smooth','B_M smooth','B_N smooth','Location','northeast')
datetick('x','HH:MM');
set(gca, 'FontSize', FontSize)

% close

%% functions
function [V,D,ratio21,ratio32] = get_LMN(Br,Bt,Bn,fld_Epoch,calc_beg_epoch,calc_end_epoch)
    Br_calc = Br(find(fld_Epoch >= calc_beg_epoch & fld_Epoch <= calc_end_epoch)); 
    Br_calc(isnan(Br_calc)) = mean(Br_calc,'omitnan');
    Bt_calc = Bt(find(fld_Epoch >= calc_beg_epoch & fld_Epoch <= calc_end_epoch)); 
    Bt_calc(isnan(Bt_calc)) = mean(Bt_calc,'omitnan');
    Bn_calc = Bn(find(fld_Epoch >= calc_beg_epoch & fld_Epoch <= calc_end_epoch)); 
    Bn_calc(isnan(Bn_calc)) = mean(Bn_calc,'omitnan');
    M = [mean(Br_calc.*Br_calc) - mean(Br_calc)*mean(Br_calc), mean(Br_calc.*Bt_calc) - mean(Br_calc)*mean(Bt_calc), mean(Br_calc.*Bn_calc) - mean(Br_calc)*mean(Bn_calc);
         mean(Bt_calc.*Br_calc) - mean(Bt_calc)*mean(Br_calc), mean(Bt_calc.*Bt_calc) - mean(Bt_calc)*mean(Bt_calc), mean(Bt_calc.*Bn_calc) - mean(Bt_calc)*mean(Bn_calc);
         mean(Bn_calc.*Br_calc) - mean(Bn_calc)*mean(Br_calc), mean(Bn_calc.*Bt_calc) - mean(Bn_calc)*mean(Bt_calc), mean(Bn_calc.*Bn_calc) - mean(Bn_calc)*mean(Bn_calc)];
    [V, D] = eig(M);
    lambda1 = D(1,1);lambda2 = D(2,2); lambda3 = D(3,3); % lambda1 < lambda2 < lambda3; order: N, M, L
    ratio21 = lambda2 / lambda1; ratio32 = lambda3 / lambda2; % ratio21 usually should be limited
end

function [smoothed_data] = smooth1d(data, method, window_size, sigma)
    data = reshape(data, [], 1);
    switch lower(method)
        case 'movmean'
            pad_width = floor(window_size/2);
            data_padded = padarray(data, [pad_width, 0], 'symmetric', 'both');
            smoothed_data = movmean(data_padded, window_size);
            smoothed_data = smoothed_data(pad_width+1 : end-pad_width);
        case 'gaussian'
            half_win = floor(window_size/2);
            x = -half_win : half_win;
            gauss_kernel = exp(-x.^2/(2*sigma^2));
            gauss_kernel = gauss_kernel / sum(gauss_kernel);
            smoothed_data = conv(data, gauss_kernel, 'same');
    end
    if isrow(data)
        smoothed_data = smoothed_data';
    end
end