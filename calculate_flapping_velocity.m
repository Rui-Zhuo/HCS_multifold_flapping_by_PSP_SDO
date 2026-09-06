clear; close all;
cfg = hcs_config();
encounter = 7;
dTime = 0.02;
save_or_not = 0;
data_store = [fullfile(cfg.psp_data_root, ['Encounter ', num2str(encounter)]), filesep];
data_save = [fullfile(cfg.work_root, ['Encounter ', num2str(encounter)], ...
    'flapping_events'), filesep];
%% which encounter
if encounter == 5 % 2020-06-08 spi
    year = 2020;
    flap_num = 2;
    cross_num = [4,4];
    spc_or_spi = ['i';'i'];
    plot_beg_lst = ['2020-06-08 00:00:00';'2020-06-08 00:00:00']; % i_flap
    plot_end_lst = ['2020-06-08 04:00:00';'2020-06-08 04:00:00']; % i_flap
    fld_lst = ['00';'00'];
end
if encounter == 6 % 2020-09-19 spc; 2020-09-21 spc
    year = 2020;
    flap_num = 2;
    cross_num = [7,7];
    spc_or_spi = ['c';'c'];
    plot_beg_lst = ['2020-09-21 01:00:00';'2020-09-21 03:50:00']; % i_flap
    plot_end_lst = ['2020-09-21 05:00:00';'2020-09-21 04:10:00']; % i_flap
    fld_lst = ['00';'00'];
end
if encounter == 7 % 2021-01-17 spi
    year = 2021;
    flap_num = 2;
    cross_num = [5,5];
    spc_or_spi = ['i';'i'];
    plot_beg_lst  = ['2021-01-17 13:00:00';'2021-01-17 13:00:00']; % i_flap
    plot_end_lst  = ['2021-01-17 15:00:00';'2021-01-17 15:00:00']; % i_flap
    fld_lst = ['12';'12'];
end
if encounter == 8 % 2021-04-29 spi
    year = 2021;
    flap_num = 2;
    cross_num = [5,5];
    spc_or_spi = ['i';'i'];
    plot_beg_lst  = ['2021-04-29 07:30:00';'2021-04-29 09:24:00']; % i_flap
    plot_end_lst  = ['2021-04-29 10:30:00';'2021-04-29 09:26:30']; % i_flap
    fld_lst = ['06';'06'];
end
%% calculate flapping_velocity
for i_flap = 1:1%flap_num
    time_plot_beg = datenum(plot_beg_lst(i_flap,:));
    time_plot_end = datenum(plot_end_lst(i_flap,:));
    %% spc_data: Vrtn_spc, sc_pos, sc_vel, carr_lat, carr_lon
    spc_file = ['psp_swp_spc_l3i_',plot_beg_lst(i_flap,1:4),plot_beg_lst(i_flap,6:7),plot_beg_lst(i_flap,9:10),'_v02.cdf'];
    spc_dir = [data_store,spc_file];
    spc_info = spdfcdfinfo(spc_dir);
    
    spc_Epoch = spdfcdfread(spc_dir,'Variables','Epoch');
    general_flag = spdfcdfread(spc_dir,'Variables','general_flag');
    Vrtn_spc = spdfcdfread(spc_dir,'Variables','vp_moment_RTN');
    Vr_spc = Vrtn_spc(:,1)'; Vr_spc(abs(Vr_spc)>1e3) = nan;
    Vt_spc = Vrtn_spc(:,2)'; Vt_spc(abs(Vt_spc)>1e3) = nan;
    Vn_spc = Vrtn_spc(:,3)'; Vn_spc(abs(Vn_spc)>1e3) = nan;
    carr_lat = spdfcdfread(spc_dir,'Variables','carr_latitude');
    carr_lon = spdfcdfread(spc_dir,'Variables','carr_longitude');
    sc_pos_HCI = spdfcdfread(spc_dir,'Variables','sc_pos_HCI');
    sc_vel_HCI = spdfcdfread(spc_dir,'Variables','sc_vel_HCI');
    sc_pos_HCIx = sc_pos_HCI(:,1); sc_vel_HCIx = sc_vel_HCI(:,1);
    sc_pos_HCIy = sc_pos_HCI(:,2); sc_vel_HCIy = sc_vel_HCI(:,2);
    sc_pos_HCIz = sc_pos_HCI(:,3); sc_vel_HCIz = sc_vel_HCI(:,3);
    sc_distan = sqrt(sc_pos_HCIx.^2 + sc_pos_HCIy.^2 + sc_pos_HCIz.^2);
    [sc_vel_RTNr,sc_vel_RTNt,sc_vel_RTNn] = calc_HCI2SCRTN(sc_vel_HCIx,sc_vel_HCIy,sc_vel_HCIz,sc_pos_HCIx,sc_pos_HCIy,sc_pos_HCIz);
    sc_vel_RTN = [sc_vel_RTNr,sc_vel_RTNt,sc_vel_RTNn];
    
    spc_plot_index = find(spc_Epoch >= time_plot_beg & spc_Epoch <= time_plot_end);
    spc_Epoch_plot = spc_Epoch(spc_plot_index);
    Vr_spc_plot = Vr_spc(spc_plot_index); Vt_spc_plot = Vt_spc(spc_plot_index); Vn_spc_plot = Vn_spc(spc_plot_index);
    carr_lat_plot = carr_lat(spc_plot_index); carr_lon_plot = carr_lon(spc_plot_index);
    sc_distan_plot = sc_distan(spc_plot_index);
    %% spi_data: Vrtn_spi, Np, Tp
    spi_file = ['psp_swp_spi_sf00_l3_mom_inst_',plot_beg_lst(i_flap,1:4),plot_beg_lst(i_flap,6:7),plot_beg_lst(i_flap,9:10),'_v03.cdf'];
    spi_dir = [data_store,spi_file];
    spi_info = spdfcdfinfo(spi_dir);
    
    spi_Epoch = spdfcdfread(spi_dir,'Variables','Epoch');
    Vsc_spi = spdfcdfread(spi_dir,'Variables','VEL');
    Vsc_spi(abs(Vsc_spi)>1e3) = nan;
    SC2RTN = spdfcdfread(spi_dir,'Variables','ROTMAT_SC_INST');
    Vrtn2sc_spi = (SC2RTN) * Vsc_spi';
        %% switch Vrtn to inertial RTN frame
    sc_vel_RTNr_interp = interp1(spc_Epoch,sc_vel_RTNr,spi_Epoch,'pchip');
    sc_vel_RTNt_interp = interp1(spc_Epoch,sc_vel_RTNt,spi_Epoch,'pchip');
    sc_vel_RTNn_interp = interp1(spc_Epoch,sc_vel_RTNn,spi_Epoch,'pchip');
    Vr2sc_spi = -Vrtn2sc_spi(3,:); Vr_spi = Vr2sc_spi + sc_vel_RTNr_interp';
    Vt2sc_spi =  Vrtn2sc_spi(1,:); Vt_spi = Vt2sc_spi + sc_vel_RTNt_interp';
    Vn2sc_Spi = -Vrtn2sc_spi(2,:); Vn_spi = Vn2sc_Spi + sc_vel_RTNn_interp';
    Vrtn_spi = [Vr_spi;Vt_spi;Vn_spi];
    
    spi_plot_index = find(spi_Epoch >= time_plot_beg & spi_Epoch <= time_plot_end);
    spi_Epoch_plot = spi_Epoch(spi_plot_index);
    Vr_spi_plot = Vr_spi(spi_plot_index); Vt_spi_plot = Vt_spi(spi_plot_index); Vn_spi_plot = Vn_spi(spi_plot_index);
    %% interpolation
    nTime = floor((time_plot_end - time_plot_beg)*86400/dTime);
    std_time = linspace(time_plot_beg, time_plot_end, nTime);
    
    std_Epoch = interp1(spi_Epoch_plot,spi_Epoch_plot,std_time,'pchip');
    Vr_spi_interp = interp1(spi_Epoch_plot,Vr_spi_plot,std_time,'linear'); Vr_spi_interp(isnan(Vr_spi_interp)) = mean(Vr_spi_interp,'omitnan');
    Vt_spi_interp = interp1(spi_Epoch_plot,Vt_spi_plot,std_time,'linear'); Vt_spi_interp(isnan(Vt_spi_interp)) = mean(Vt_spi_interp,'omitnan');
    Vn_spi_interp = interp1(spi_Epoch_plot,Vn_spi_plot,std_time,'linear'); Vn_spi_interp(isnan(Vn_spi_interp)) = mean(Vn_spi_interp,'omitnan');
    Vr_spc_interp = interp1(spc_Epoch_plot,Vr_spc_plot,std_time,'linear'); Vr_spc_interp(isnan(Vr_spc_interp)) = mean(Vr_spc_interp,'omitnan');
    Vt_spc_interp = interp1(spc_Epoch_plot,Vt_spc_plot,std_time,'linear'); Vt_spc_interp(isnan(Vt_spc_interp)) = mean(Vt_spc_interp,'omitnan');
    Vn_spc_interp = interp1(spc_Epoch_plot,Vn_spc_plot,std_time,'linear'); Vn_spc_interp(isnan(Vn_spc_interp)) = mean(Vn_spc_interp,'omitnan');
    carr_lat_interp = interp1(spc_Epoch_plot,carr_lat_plot,std_time,'linear');
    carr_lon_interp = interp1(spc_Epoch_plot,carr_lon_plot,std_time,'linear');
    sc_distan_interp = interp1(spc_Epoch_plot,sc_distan_plot,std_time,'linear');
    Vrtn_spi_interp = cat(2,Vr_spi_interp',Vt_spi_interp',Vn_spi_interp');
    Vrtn_spc_interp = cat(2,Vr_spc_interp',Vt_spc_interp',Vn_spc_interp');
    %% import e_LMN
    e_LMN_mat = importdata([data_save,'e_LMN.csv']);
    e_L_mat = e_LMN_mat(:,1:cross_num(i_flap));
    e_M_mat = e_LMN_mat(:,cross_num(i_flap)+1:2*cross_num(i_flap));
    e_N_mat = e_LMN_mat(:,2*cross_num(i_flap)+1:3*cross_num(i_flap));
    %% crossing list
    if encounter == 5
        if i_flap == 1
            jump_beg_lst = ['2020-06-08 00:00:55';'2020-06-08 00:01:25';'2020-06-08 00:01:40';'2020-06-08 00:02:40']; % i_cross
            jump_end_lst = ['2020-06-08 00:01:05';'2020-06-08 00:01:35';'2020-06-08 00:01:50';'2020-06-08 00:02:50']; % i_cross
        end
        if i_flap == 2
            jump_beg_lst = ['2020-06-08 00:00:55';'2020-06-08 00:01:25';'2020-06-08 00:01:40';'2020-06-08 00:02:40']; % i_cross
            jump_end_lst = ['2020-06-08 00:01:05';'2020-06-08 00:01:35';'2020-06-08 00:01:50';'2020-06-08 00:02:50']; % i_cross
        end
    end
    if encounter == 6
        if i_flap == 1
            jump_beg_lst = ['2020-09-21 01:27:00';'2020-09-21 01:54:00';'2020-09-21 02:12:00';'2020-09-21 02:27:00';'2020-09-21 02:29:50';'2020-09-21 03:08:30';'2020-09-21 03:58:50']; % i_cross
            jump_end_lst = ['2020-09-21 01:30:00';'2020-09-21 01:59:00';'2020-09-21 02:21:00';'2020-09-21 02:29:30';'2020-09-21 02:31:10';'2020-09-21 03:10:00';'2020-09-21 03:59:15']; % i_cross
        end
        if i_flap == 2
            jump_beg_lst = ['2020-09-21 01:27:00';'2020-09-21 01:54:00';'2020-09-21 02:12:00';'2020-09-21 02:27:00';'2020-09-21 02:29:50';'2020-09-21 03:08:30';'2020-09-21 03:58:50']; % i_cross
            jump_end_lst = ['2020-09-21 01:30:00';'2020-09-21 01:59:00';'2020-09-21 02:21:00';'2020-09-21 02:29:30';'2020-09-21 02:31:10';'2020-09-21 03:10:00';'2020-09-21 03:59:15']; % i_cross
        end
    end
    if encounter == 7
        if i_flap == 1
            jump_beg_lst = ['2021-01-17 13:14:10';'2021-01-17 14:02:55';'2021-01-17 14:05:28';'2021-01-17 14:14:00';'2021-01-17 14:40:34']; % i_cross
            jump_end_lst = ['2021-01-17 13:30:55';'2021-01-17 14:03:10';'2021-01-17 14:05:34';'2021-01-17 14:20:00';'2021-01-17 14:41:04']; % i_cross
        end
        if i_flap == 2
            jump_beg_lst = ['2021-01-17 13:14:10';'2021-01-17 14:02:55';'2021-01-17 14:05:28';'2021-01-17 14:14:00';'2021-01-17 14:40:34']; % i_cross
            jump_end_lst = ['2021-01-17 13:30:55';'2021-01-17 14:03:10';'2021-01-17 14:05:34';'2021-01-17 14:20:00';'2021-01-17 14:41:04']; % i_cross
        end
    end
    if encounter == 8
        if i_flap == 1
            jump_beg_lst = ['2021-04-29 08:14:30';'2021-04-29 09:13:05';'2021-04-29 09:15:00';'2021-04-29 09:24:30';'2021-04-29 09:34:00']; % i_cross
            jump_end_lst = ['2021-04-29 08:28:00';'2021-04-29 09:14:15';'2021-04-29 09:15:20';'2021-04-29 09:25:15';'2021-04-29 10:11:00']; % i_cross
        end
        if i_flap == 2
            jump_beg_lst = ['2021-04-29 08:14:30';'2021-04-29 09:13:05';'2021-04-29 09:15:00';'2021-04-29 09:24:30';'2021-04-29 09:34:00']; % i_cross
            jump_end_lst = ['2021-04-29 08:28:00';'2021-04-29 09:14:15';'2021-04-29 09:15:20';'2021-04-29 09:25:15';'2021-04-29 10:11:00']; % i_cross
        end
    end
    %% calculate jump velocity
    if spc_or_spi(i_flap,:) == 'c'
        VN_spc_mean = zeros(1,cross_num(i_flap));
    end
    if spc_or_spi(i_flap,:) == 'i'
        VN_spi_mean = zeros(1,cross_num(i_flap));
    end
    VN_rotate_mean = zeros(1,cross_num(i_flap));
    VN_flap_mean =  zeros(1,cross_num(i_flap));
    
    for i_cross = 1:cross_num(i_flap)
        time_jump_beg = datenum(jump_beg_lst(i_flap,:));
        time_jump_end = datenum(jump_end_lst(i_flap,:));
        
        jump_index = find(std_Epoch >= time_jump_beg & std_Epoch <= time_jump_end);
        jump_Epoch = std_Epoch(jump_index);
        Vrtn_spi_jump = Vrtn_spi_interp(jump_index,:);
        Vrtn_spc_jump = Vrtn_spc_interp(jump_index,:);
        carr_lat_jump = carr_lat_interp(jump_index);
        sc_distan_jump = sc_distan_interp(jump_index);
        
        V = cat(2,e_N_mat(:,i_cross),e_M_mat(:,i_cross),e_L_mat(:,i_cross));
        Vnml_spi_jump = Vrtn_spi_jump * V;
        VL_spi_jump = Vnml_spi_jump(:,3);
        VM_spi_jump = Vnml_spi_jump(:,2);
        VN_spi_jump = Vnml_spi_jump(:,1);
        Vnml_spc_jump = Vrtn_spc_jump * V;
        VL_spc_jump = Vnml_spc_jump(:,3);
        VM_spc_jump = Vnml_spc_jump(:,2);
        VN_spc_jump = Vnml_spc_jump(:,1);
        
        Vt_rotate = sun_rotation(carr_lat_jump,sc_distan_jump); % rotation velocity (in rtn't' direction)
        Vr_rotate = zeros(1,length(jump_index)); Vn_rotate = zeros(1,length(jump_index));
        Vrtn_rotate = cat(2,Vr_rotate',Vt_rotate',Vn_rotate');
        Vnml_rotate = Vrtn_rotate * V;
        VL_rotate = Vnml_rotate(:,3);
        VM_rotate = Vnml_rotate(:,2);
        VN_rotate = Vnml_rotate(:,1);
        
        if spc_or_spi(i_flap,:) == 'c'
            VN_flap = VN_spc_jump - VN_rotate;
            VN_spc_mean(i_cross) = mean(VN_spc_jump);
        end
        if spc_or_spi(i_flap,:) == 'i'
            VN_flap = VN_spi_jump - VN_rotate;
            VN_spi_mean(i_cross) = mean(VN_spi_jump);
        end
        VN_rotate_mean(i_cross) = mean(VN_rotate);
        VN_flap_mean(i_cross) = mean(VN_flap);
    end
    %% orbit
    rr = zeros(1,cross_num(i_flap));
    tt = 200 * (1:cross_num(i_flap));
    nn = zeros(1,cross_num(i_flap));
    rr_min = -100; rr_max = 100;
    tt_min =  0; tt_max = 200 * (cross_num(i_flap) + 1);
    nn_min = -100; nn_max = 100;
    %% plot VN_flap vectors
    VN_flap_vec = VN_flap_mean .* e_N_mat;
    
    LineWidth = 2;

    subplot(2,2,1)
    quiver3(rr(1),tt(1),nn(1),rr(cross_num(i_flap))-rr(1),tt(cross_num(i_flap))-tt(1),nn(cross_num(i_flap))-nn(1),'off',':k','LineWidth',LineWidth); hold on;
    quiver3(rr,tt,nn,VN_flap_vec(1,:),VN_flap_vec(2,:),VN_flap_vec(3,:),'off','b','LineWidth',LineWidth); grid on
    xlabel('R'); ylabel('T'); zlabel('N');
    axis equal
%     xlim([rr_min rr_max]); ylim([tt_min tt_max]); zlim([nn_min nn_max]);
    subtitle('3-D vector')
    
    subplot(2,2,2)
    quiver(tt(1),rr(1),tt(cross_num(i_flap))-tt(1),rr(cross_num(i_flap))-rr(1),'off',':k','LineWidth',LineWidth); hold on;
    quiver(tt,rr,VN_flap_vec(2,:),VN_flap_vec(1,:),'off','b','LineWidth',LineWidth); grid on
    xlabel('T'); ylabel('R');
    axis equal
%     xlim([tt_min tt_max]); ylim([rr_min rr_max]);
    subtitle('T-R projection')
    
    subplot(2,2,3)
    quiver(tt(1),nn(1),tt(cross_num(i_flap))-tt(1),nn(cross_num(i_flap))-nn(1),'off',':k','LineWidth',LineWidth); hold on;
    quiver(tt,nn,VN_flap_vec(2,:),VN_flap_vec(3,:),'off','b','LineWidth',LineWidth); grid on
    xlabel('T'); ylabel('N');
    axis equal
%     xlim([tt_min tt_max]); ylim([nn_min nn_max]);
    subtitle('T-N projection')
    
    subplot(2,2,4)
    quiver(rr(1),nn(1),rr(cross_num(i_flap))-rr(1),nn(cross_num(i_flap))-nn(1),'off',':k','LineWidth',LineWidth); hold on;
    quiver(rr,nn,VN_flap_vec(1,:),VN_flap_vec(3,:),'off','b','LineWidth',LineWidth); grid on
    xlabel('R'); ylabel('N');
    axis equal
%     xlim([rr_min rr_max]); ylim([nn_min nn_max]);
    subtitle('R-N projection')
    
    sgtitle('V_{flapping} in RTN coordinate');
end
%% functions
function [x_RTN,y_RTN,z_RTN] = calc_HCI2SCRTN(x_HCI,y_HCI,z_HCI,SC_HCIx,SC_HCIy,SC_HCIz)
% Change the coordiantes xyz in HCI frame to xyz in spacecraft RTN frame
%   input: x_HCI, y_HCI, z_HCI, the velocity in HCI frame (km/s)
%          SC_HCIx, SC_HCIy, SC_HCIz, the sapcecraft position in HCI frame (km)
%   output: x_RTN, y_RTN, z_RTN,the velocity in SC RTN frame (km/s)
%   This function does not consider the move of the origin (using for velocity conversion)
    num = length(x_HCI);
    xyz_RTN = zeros(num,3);
    for i = 1:num
        Q = zeros(3,3);
        x1 = [1 0 0];
        y1 = [0 1 0];
        z1 = [0 0 1];
        x2 = [SC_HCIx(i),SC_HCIy(i),SC_HCIz(i)];
        if norm(x2)~= 0
            x2 = x2/norm(x2);
        end
        y2 = cross(z1,x2);
        if norm(y2)~= 0
            y2 = y2/norm(y2);
        end
        z2 = cross(x2,y2);
        if norm(z2)~= 0
            z2 = z2/norm(z2);
        end
        Q(1,1) = dot(x2,x1); Q(1,2) = dot(x2,y1); Q(1,3) = dot(x2,z1);
        Q(2,1) = dot(y2,x1); Q(2,2) = dot(y2,y1); Q(2,3) = dot(y2,z1);
        Q(3,1) = dot(z2,x1); Q(3,2) = dot(z2,y1); Q(3,3) = dot(z2,z1);
        xyz_RTN(i,:) = Q*[x_HCI(i);y_HCI(i);z_HCI(i)];
    end
    x_RTN = xyz_RTN(:,1); y_RTN = xyz_RTN(:,2); z_RTN = xyz_RTN(:,3);
end
function Vt_rotate = sun_rotation(carr_lat,sc_distan)
    deg2rad = 3.1415926535 / 180; % from degree to radian
    day2sec = 24*60*60;
    omega1 = 14.713; omega2 = -2.396; omega3 = -1.787; % unit deg./day
    omega = (omega1 + omega2 * sind(carr_lat).^2 + omega3 * sind(carr_lat).^4) * deg2rad / day2sec; % sun rotation
        
    Vt_rotate = omega .* sc_distan;
end
