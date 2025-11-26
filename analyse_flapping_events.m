clear; close all;
encounter = 7;
dTime = 0.02;
save_or_not = 0;
mean_or_local = 0;
which_cross = 0;
data_store = ['E:\Research\Data\PSP\Encounter ',num2str(encounter),'\'];
data_save  = ['E:\Research\Work\HCS_multifold_flapping_by_PSP_SDO\Encounter ',num2str(encounter),'\'];
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
    plot_beg_lst = ['2021-01-17 13:00:01';'2021-01-17 14:30:01']; % i_flap
    plot_end_lst = ['2021-01-17 15:00:00';'2021-01-17 14:50:00']; % i_flap
    cros_beg_lst = ['2021-01-17 13:10:00';'2021-01-17 13:55:00';'2021-01-17 13:55:00';'2021-01-17 14:10:00';'2021-01-17 14:30:00'];
    cros_end_lst = ['2021-01-17 13:35:00';'2021-01-17 14:15:00';'2021-01-17 14:15:00';'2021-01-17 14:30:00';'2021-01-17 14:50:00'];
    mark_beg_lst = ['2021-01-17 13:14:30';'2021-01-17 14:02:30';'2021-01-17 14:05:00';'2021-01-17 14:14:30';'2021-01-17 14:40:00'];
    mark_end_lst = ['2021-01-17 13:30:00';'2021-01-17 14:04:00';'2021-01-17 14:06:00';'2021-01-17 14:20:00';'2021-01-17 14:41:30'];
    if mean_or_local == 1
        plot_beg_lst(2,:) = cros_beg_lst(which_cross,:);
        plot_end_lst(2,:) = cros_end_lst(which_cross,:);
    end
    fld_lst = ['12';'12'];
end
if encounter == 8 % 2021-04-29 spi
    year = 2021;
    flap_num = 2;
    cross_num = [5,4,4];
    spc_or_spi = ['i';'i';'i'];
    plot_beg_lst = ['2021-04-29 07:30:00';'2021-04-28 19:00:00';'2021-04-28 19:40:01']; % i_flap
    plot_end_lst = ['2021-04-29 10:30:00';'2021-04-28 21:00:00';'2021-04-28 19:49:59']; % i_flap
    cros_beg_lst = ['2021-04-29 08:00:00';'2021-04-29 09:10:00';'2021-04-29 09:10:00';'2021-04-29 09:20:00';'2021-04-29 09:26:00'];
    cros_end_lst = ['2021-04-29 08:40:00';'2021-04-29 09:20:00';'2021-04-29 09:20:00';'2021-04-29 09:30:00';'2021-04-29 10:40:00'];
%     cros_beg_lst = ['2021-04-28 19:40:01';'2021-04-28 19:50:01';'2021-04-28 20:05:01';'2021-04-28 20:17:01'];
%     cros_end_lst = ['2021-04-28 19:49:59';'2021-04-28 19:59:59';'2021-04-28 20:14:59';'2021-04-28 20:26:59'];
    if mean_or_local == 1
        plot_beg_lst(3,:) = cros_beg_lst(which_cross,:);
        plot_end_lst(3,:) = cros_end_lst(which_cross,:);
    end
    fld_lst = ['06';'18';'18'];
end
for i_flap = 1:1%flap_num
    close all;
    %% which flapping
    if encounter == 5
        if i_flap == 1
            calc_beg_lst = ['2020-06-08 00:00:55';'2020-06-08 00:01:25';'2020-06-08 00:01:40';'2020-06-08 00:02:40']; % i_cross
            calc_end_lst = ['2020-06-08 00:01:05';'2020-06-08 00:01:35';'2020-06-08 00:01:50';'2020-06-08 00:02:50']; % i_cross
        end
        if i_flap == 2
            calc_beg_lst = ['2020-06-08 00:00:55';'2020-06-08 00:01:25';'2020-06-08 00:01:40';'2020-06-08 00:02:40']; % i_cross
            calc_end_lst = ['2020-06-08 00:01:05';'2020-06-08 00:01:35';'2020-06-08 00:01:50';'2020-06-08 00:02:50']; % i_cross
        end
    end
    if encounter == 6
        if i_flap == 1
            calc_beg_lst = ['2020-09-21 01:27:00';'2020-09-21 01:54:00';'2020-09-21 02:12:00';'2020-09-21 02:27:00';'2020-09-21 02:29:50';'2020-09-21 03:08:30';'2020-09-21 03:58:50']; % i_cross
            calc_end_lst = ['2020-09-21 01:30:00';'2020-09-21 01:59:00';'2020-09-21 02:21:00';'2020-09-21 02:29:30';'2020-09-21 02:31:10';'2020-09-21 03:10:00';'2020-09-21 03:59:15']; % i_cross
        end
        if i_flap == 2
            calc_beg_lst = ['2020-09-21 01:27:00';'2020-09-21 01:54:00';'2020-09-21 02:12:00';'2020-09-21 02:27:00';'2020-09-21 02:29:50';'2020-09-21 03:08:30';'2020-09-21 03:58:50']; % i_cross
            calc_end_lst = ['2020-09-21 01:30:00';'2020-09-21 01:59:00';'2020-09-21 02:21:00';'2020-09-21 02:29:30';'2020-09-21 02:31:10';'2020-09-21 03:10:00';'2020-09-21 03:59:15']; % i_cross
        end
    end
    if encounter == 7
        if i_flap == 1
            calc_beg_lst = ['2021-01-17 13:29:40';'2021-01-17 14:02:00';'2021-01-17 14:05:20';'2021-01-17 14:14:05';'2021-01-17 14:40:20']; % i_cross
            calc_end_lst = ['2021-01-17 13:31:00';'2021-01-17 14:04:00';'2021-01-17 14:06:10';'2021-01-17 14:18:50';'2021-01-17 14:40:50']; % i_cross
        end
        if i_flap == 2
            calc_beg_lst = ['2021-01-17 13:13:00';'2021-01-17 14:02:00';'2021-01-17 14:05:20';'2021-01-17 14:14:05';'2021-01-17 14:40:20']; % i_cross
            calc_end_lst = ['2021-01-17 13:15:00';'2021-01-17 14:04:00';'2021-01-17 14:06:10';'2021-01-17 14:18:50';'2021-01-17 14:40:50']; % i_cross
            lead_beg_lst = ['2021-01-17 13:13:00';'2021-01-17 13:13:00';'2021-01-17 13:13:00';'2021-01-17 14:10:00';'2021-01-17 14:35:00'];
            lead_end_lst = ['2021-01-17 13:14:00';'2021-01-17 13:13:00';'2021-01-17 13:13:00';'2021-01-17 14:15:00';'2021-01-17 14:40:00'];
            tral_beg_lst = ['2021-01-17 13:15:00';'2021-01-17 13:13:00';'2021-01-17 13:13:00';'2021-01-17 14:23:00';'2021-01-17 14:45:00'];
            tral_end_lst = ['2021-01-17 13:16:00';'2021-01-17 13:13:00';'2021-01-17 13:13:00';'2021-01-17 14:28:00';'2021-01-17 14:50:00'];
        end
        % if choose the leading-edge of crossing 1, it should be calc_beg('2021-01-17 13:13:00') and calc_end('2021-01-17 13:15:00'), and the normal points to (-0.13,0.36,0.92)
    end
    if encounter == 8
        if i_flap == 1
            calc_beg_lst = ['2021-04-29 08:26:00';'2021-04-29 09:13:00';'2021-04-29 09:14:30';'2021-04-29 09:23:00';'2021-04-29 10:10:00']; % i_cross
            calc_end_lst = ['2021-04-29 08:29:00';'2021-04-29 09:14:30';'2021-04-29 09:15:30';'2021-04-29 09:27:00';'2021-04-29 10:13:00']; % i_cross
        end
        if i_flap == 2
            calc_beg_lst = ['2021-04-29 08:13:00';'2021-04-29 09:13:00';'2021-04-29 09:14:30';'2021-04-29 09:23:00';'2021-04-29 10:10:00']; % i_cross
            calc_end_lst = ['2021-04-29 08:16:00';'2021-04-29 09:14:30';'2021-04-29 09:15:30';'2021-04-29 09:27:00';'2021-04-29 10:13:00']; % i_cross
            lead_beg_lst = ['2021-04-29 08:09:00';'2021-04-29 08:29:00';'2021-04-29 08:29:00';'2021-04-29 09:20:00';'2021-04-29 09:26:00'];
            lead_end_lst = ['2021-04-29 08:14:00';'2021-04-29 08:29:00';'2021-04-29 08:29:00';'2021-04-29 09:24:00';'2021-04-29 09:35:00'];
            tral_beg_lst = ['2021-04-29 08:29:00';'2021-04-29 08:29:00';'2021-04-29 08:29:00';'2021-04-29 09:26:00';'2021-04-29 10:31:00'];
            tral_end_lst = ['2021-04-29 08:34:00';'2021-04-29 08:29:00';'2021-04-29 08:29:00';'2021-04-29 09:30:00';'2021-04-29 10:40:00'];
        end
        if i_flap == 3
            calc_beg_lst = ['2021-04-28 19:44:01';'2021-04-28 19:53:31';'2021-04-28 20:09:51';'2021-04-28 20:17:01']; % i_cross
            calc_end_lst = ['2021-04-28 19:48:59';'2021-04-28 19:59:59';'2021-04-28 20:12:59';'2021-04-28 20:26:59']; % i_cross
%             lead_beg_lst = ['2021-04-28 19:00:00';'2021-04-28 19:00:00';'2021-04-28 19:00:00';'2021-04-28 19:00:00'];
%             lead_end_lst = ['2021-04-28 21:00:00';'2021-04-28 21:00:00';'2021-04-28 21:00:00';'2021-04-28 21:00:00'];
%             tral_beg_lst = ['2021-04-28 19:00:00';'2021-04-28 19:00:00';'2021-04-28 19:00:00';'2021-04-28 19:00:00'];
%             tral_end_lst = ['2021-04-28 21:00:00';'2021-04-28 21:00:00';'2021-04-28 21:00:00';'2021-04-28 21:00:00'];
        end
    end
    %% begin time and end time
    time_plot_beg = datenum(plot_beg_lst(i_flap,:));
    time_plot_end = datenum(plot_end_lst(i_flap,:));
    %% fld_data: Brtn
    fld_file = ['psp_fld_l2_mag_rtn_',plot_beg_lst(i_flap,1:4),plot_beg_lst(i_flap,6:7),plot_beg_lst(i_flap,9:10),fld_lst(i_flap,:),'_v02.cdf'];
    fld_dir = [data_store,fld_file];
    fld_info = spdfcdfinfo(fld_dir);
    
    fld_Epoch = spdfcdfread(fld_dir,'Variables','epoch_mag_RTN');
    Brtn = spdfcdfread(fld_dir,'Variables','psp_fld_l2_mag_RTN');
    Brtn(abs(Brtn)>1e3) = nan;
    Br = Brtn(:,1); Bt = Brtn(:,2); Bn = Brtn(:,3); B_mod = sqrt(Br.^2 + Bt.^2 + Bn.^2);
    
    fld_plot_index = find(fld_Epoch >= time_plot_beg & fld_Epoch <= time_plot_end);
    fld_Epoch_plot = fld_Epoch(fld_plot_index);
    Br_plot = Br(fld_plot_index); Bt_plot = Bt(fld_plot_index); Bn_plot = Bn(fld_plot_index); B_mod_plot = B_mod(fld_plot_index);
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
    Np_spc = spdfcdfread(spc_dir,'Variables','np_moment');
    Np_spc(abs(Np_spc)>1e4) = nan;
    Tp_spc = spdfcdfread(spc_dir,'Variables','wp_moment');
    Tp_spc(abs(Tp_spc)>6e2) = nan;
    carr_lat = spdfcdfread(spc_dir,'Variables','carr_latitude');
    carr_lon = spdfcdfread(spc_dir,'Variables','carr_longitude');
    sc_pos_HCI = spdfcdfread(spc_dir,'Variables','sc_pos_HCI');
    sc_vel_HCI = spdfcdfread(spc_dir,'Variables','sc_vel_HCI');
    sc_pos_HCIx = sc_pos_HCI(:,1); sc_vel_HCIx = sc_vel_HCI(:,1);
    sc_pos_HCIy = sc_pos_HCI(:,2); sc_vel_HCIy = sc_vel_HCI(:,2);
    sc_pos_HCIz = sc_pos_HCI(:,3); sc_vel_HCIz = sc_vel_HCI(:,3);
    [sc_vel_RTNr,sc_vel_RTNt,sc_vel_RTNn] = calc_HCI2SCRTN(sc_vel_HCIx,sc_vel_HCIy,sc_vel_HCIz,sc_pos_HCIx,sc_pos_HCIy,sc_pos_HCIz);
    sc_vel_RTN = [sc_vel_RTNr,sc_vel_RTNt,sc_vel_RTNn];
    
    spc_plot_index = find(spc_Epoch >= time_plot_beg & spc_Epoch <= time_plot_end);
    spc_Epoch_plot = spc_Epoch(spc_plot_index);
    Vr_spc_plot = Vr_spc(spc_plot_index); Vt_spc_plot = Vt_spc(spc_plot_index); Vn_spc_plot = Vn_spc(spc_plot_index);
    Np_spc_plot = Np_spc(spc_plot_index); Tp_spc_plot = Tp_spc(spc_plot_index);
    carr_lat_plot = carr_lat(spc_plot_index); carr_lon_plot = carr_lon(spc_plot_index);
        %% get rid of general_flag ~= 0
        bad = find(general_flag~=0);
        badr = find(isnan(Vr_spc)); badt = find(isnan(Vt_spc)); badn = find(isnan(Vn_spc));
        badN = find(isnan(Np_spc)); badT = find(isnan(Tp_spc));
        
        spc_Epoch_good = spc_Epoch; spc_Epoch_good(cat(2,bad',badr)) = [];
        Vr_spc_good = Vr_spc; Vr_spc_good(cat(2,bad',badr)) = [];
        if isempty(Vr_spc_good) == 0
            Vr_spc = interp1(spc_Epoch_good,Vr_spc_good,spc_Epoch)';
        end
        
        spc_Epoch_good = spc_Epoch; spc_Epoch_good(cat(2,bad',badt)) = [];
        Vt_spc_good = Vt_spc; Vt_spc_good(cat(2,bad',badt)) = [];
        if isempty(Vt_spc_good) == 0
            Vt_spc = interp1(spc_Epoch_good,Vt_spc_good,spc_Epoch)';
        end
        
        spc_Epoch_good = spc_Epoch; spc_Epoch_good(cat(2,bad',badn)) = [];
        Vn_spc_good = Vn_spc; Vn_spc_good(cat(2,bad',badn)) = [];
        if isempty(Vn_spc_good) == 0
            Vn_spc = interp1(spc_Epoch_good,Vn_spc_good,spc_Epoch)';
        end
        
        spc_Epoch_good = spc_Epoch; spc_Epoch_good(cat(2,bad',badN')) = [];
        Np_spc_good = Np_spc; Np_spc_good(cat(2,bad',badN')) = [];
        if isempty(Np_spc_good) == 0
            Np_spc = interp1(spc_Epoch_good,Np_spc_good,spc_Epoch)';
        end
        
        spc_Epoch_good = spc_Epoch; spc_Epoch_good(cat(2,bad',badT')) = [];
        Tp_spc_good = Np_spc; Tp_spc_good(cat(2,bad',badT')) = [];
        if isempty(Tp_spc_good) == 0
            Tp_spc = interp1(spc_Epoch_good,Tp_spc_good,spc_Epoch)';
        end
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
    Np_spi = spdfcdfread(spi_dir,'Variables','DENS');
    Np_spi(abs(Np_spi)>1e4) = nan;
    Tp_spi = spdfcdfread(spi_dir,'Variables','TEMP');
    Tp_spi(abs(Tp_spi)>1e2) = nan;
    
    spi_plot_index = find(spi_Epoch >= time_plot_beg & spi_Epoch <= time_plot_end);
    spi_Epoch_plot = spi_Epoch(spi_plot_index);
    Vr_spi_plot = Vr_spi(spi_plot_index); Vt_spi_plot = Vt_spi(spi_plot_index); Vn_spi_plot = Vn_spi(spi_plot_index);
    Np_spi_plot = Np_spi(spi_plot_index); Tp_spi_plot = Tp_spi(spi_plot_index);
    %% spe_data: PA
    spe_file = ['psp_swp_spe_sf0_l3_pad_',plot_beg_lst(i_flap,1:4),plot_beg_lst(i_flap,6:7),plot_beg_lst(i_flap,9:10),'_v03.cdf'];
    spe_dir = [data_store,spe_file];
    spe_info = spdfcdfinfo(spe_dir);
    
    spe_Epoch = spdfcdfread(spe_dir,'Variables','Epoch');
    PA_bins = spdfcdfread(spe_dir,'Variables','PITCHANGLE');
    PA_bins = PA_bins(1,:);
    EGY_bins = spdfcdfread(spe_dir,'Variables','ENERGY_VALS');
    EGY_bins = EGY_bins(1,:);
    PA_eflux_all = spdfcdfread(spe_dir,'Variables','EFLUX_VS_PA_E');
    PA_eflux = squeeze(PA_eflux_all(:,9,:)); % row 9: 314 eV
    
    spe_plot_index = find(spe_Epoch >= time_plot_beg & spe_Epoch <= time_plot_end);
    spe_Epoch_plot = spe_Epoch(spe_plot_index);
    PA_eflux_plot = PA_eflux(:,spe_plot_index);
    %% determine LMN coordinate
    V_local = zeros(3,3,cross_num(i_flap));
    ratio21_arr = zeros(1,cross_num(i_flap));
    ratio32_arr = zeros(1,cross_num(i_flap));
    e_L_arr = zeros(3,cross_num(i_flap));
    e_M_arr = zeros(3,cross_num(i_flap));
    e_N_arr = zeros(3,cross_num(i_flap));
    e_L_sum = zeros(3,1);
    
    for i_cross = 1:cross_num(i_flap)
        time_calc_beg = datenum(calc_beg_lst(i_cross,:));
        time_calc_end = datenum(calc_end_lst(i_cross,:));
        [V_sub, D_sub,ratio21,ratio32] = get_LMN(Br,Bt,Bn,fld_Epoch,time_calc_beg,time_calc_end); % ratio21 >= 3 is needed
        ratio21_arr(i_cross) = ratio21;
        ratio32_arr(i_cross) = ratio32;
        e_L_sub = V_sub(:,3)/(sqrt(dot(V_sub(:,3), V_sub(:,3))));
        e_M_sub = V_sub(:,2)/(sqrt(dot(V_sub(:,2), V_sub(:,2))));
        e_N_sub = V_sub(:,1)/(sqrt(dot(V_sub(:,1), V_sub(:,1))));
        V_local(:,:,i_cross) = [e_N_sub,e_M_sub,e_L_sub];
        e_L_arr(:,i_cross) = e_L_sub;
        e_M_arr(:,i_cross) = e_M_sub;
        e_N_arr(:,i_cross) = e_N_sub;
        e_L_sum = e_L_sum + e_L_sub;
    end
    
    e_L = e_L_sum/sqrt(dot(e_L_sum,e_L_sum)); % size 3*1
    e_M = cross([0,0,1],e_L); % size 1*3
    e_N = cross(e_L,e_M); % size 1*3
    if mean_or_local == 0
        V = cat(2,e_N',e_M',e_L);
    elseif mean_or_local == 1
        V = squeeze(V_local(:,:,which_cross));
    end
    e_LMN_arr = cat(2,e_L_arr,e_M_arr,e_N_arr);
%     csvwrite([data_save,'e_LMN.csv'],e_LMN_arr); % e_L, e_M, e_N
    %% interpolation
    nTime = floor((time_plot_end - time_plot_beg)*86400/dTime);
    std_time = linspace(time_plot_beg, time_plot_end, nTime);
    
    std_Epoch = interp1(fld_Epoch_plot,fld_Epoch_plot,std_time,'pchip');
    Br_interp = interp1(fld_Epoch_plot,Br_plot,std_time,'linear'); Br_interp(isnan(Br_interp)) = mean(Br_interp,'omitnan');
    Bt_interp = interp1(fld_Epoch_plot,Bt_plot,std_time,'linear'); Bt_interp(isnan(Bt_interp)) = mean(Bt_interp,'omitnan');
    Bn_interp = interp1(fld_Epoch_plot,Bn_plot,std_time,'linear'); Bn_interp(isnan(Bn_interp)) = mean(Bn_interp,'omitnan');
    Vr_spi_interp = interp1(spi_Epoch_plot,Vr_spi_plot,std_time,'linear'); Vr_spi_interp(isnan(Vr_spi_interp)) = mean(Vr_spi_interp,'omitnan');
    Vt_spi_interp = interp1(spi_Epoch_plot,Vt_spi_plot,std_time,'linear'); Vt_spi_interp(isnan(Vt_spi_interp)) = mean(Vt_spi_interp,'omitnan');
    Vn_spi_interp = interp1(spi_Epoch_plot,Vn_spi_plot,std_time,'linear'); Vn_spi_interp(isnan(Vn_spi_interp)) = mean(Vn_spi_interp,'omitnan');
    Vr_spc_interp = interp1(spc_Epoch_plot,Vr_spc_plot,std_time,'linear'); Vr_spc_interp(isnan(Vr_spc_interp)) = mean(Vr_spc_interp,'omitnan');
    Vt_spc_interp = interp1(spc_Epoch_plot,Vt_spc_plot,std_time,'linear'); Vt_spc_interp(isnan(Vt_spc_interp)) = mean(Vt_spc_interp,'omitnan');
    Vn_spc_interp = interp1(spc_Epoch_plot,Vn_spc_plot,std_time,'linear'); Vn_spc_interp(isnan(Vn_spc_interp)) = mean(Vn_spc_interp,'omitnan');
    Np_spi_interp = interp1(spi_Epoch_plot,Np_spi_plot,std_time,'linear'); Np_spi_interp(isnan(Np_spi_interp)) = mean(Np_spi_interp,'omitnan');
    Tp_spi_interp = interp1(spi_Epoch_plot,Tp_spi_plot,std_time,'linear'); Tp_spi_interp(isnan(Tp_spi_interp)) = mean(Tp_spi_interp,'omitnan');
    Np_spc_interp = interp1(spc_Epoch_plot,Np_spc_plot,std_time,'linear'); Np_spc_interp(isnan(Np_spc_interp)) = mean(Np_spc_interp,'omitnan');
    Tp_spc_interp = interp1(spc_Epoch_plot,Tp_spc_plot,std_time,'linear'); Tp_spc_interp(isnan(Tp_spc_interp)) = mean(Tp_spc_interp,'omitnan');
    Brtn_interp = cat(2,Br_interp',Bt_interp',Bn_interp');
    B_mod_interp = sqrt(Br_interp.^2 + Bt_interp.^2 + Bn_interp.^2);
    Vrtn_spi_interp = cat(2,Vr_spi_interp',Vt_spi_interp',Vn_spi_interp');
    Vrtn_spc_interp = cat(2,Vr_spc_interp',Vt_spc_interp',Vn_spc_interp');
    %% calculate B_LMN,V_LMN,E_LMN,Pk,Pm,Pt
    % B_LMN and V_LMN to plot
    Bnml = Brtn * V; Vnml_spi = Vrtn_spi' * V; Vnml_spc = Vrtn_spc * V;
    BL = Bnml(:,3); BL_plot = BL(fld_plot_index);
    BM = Bnml(:,2); BM_plot = BM(fld_plot_index);
    BN = Bnml(:,1); BN_plot = BN(fld_plot_index);
    VL_spi = Vnml_spi(:,3); VL_spi_plot = VL_spi(spi_plot_index);
    VM_spi = Vnml_spi(:,2); VM_spi_plot = VM_spi(spi_plot_index);
    VN_spi = Vnml_spi(:,1); VN_spi_plot = VN_spi(spi_plot_index);
    VL_spc = Vnml_spc(:,3); VL_spc_plot = VL_spc(spc_plot_index);
    VM_spc = Vnml_spc(:,2); VM_spc_plot = VM_spc(spc_plot_index);
    VN_spc = Vnml_spc(:,1); VN_spc_plot = VN_spc(spc_plot_index);
    % B_LMN and V_LMN after interpolation
    Bnml_interp = Brtn_interp * V; Vnml_spi_interp = Vrtn_spi_interp * V; Vnml_spc_interp = Vrtn_spc_interp * V;
    BL_interp = Bnml_interp(:,3);
    BM_interp = Bnml_interp(:,2);
    BN_interp = Bnml_interp(:,1);
    VL_spi_interp = Vnml_spi_interp(:,3);
    VM_spi_interp = Vnml_spi_interp(:,2);
    VN_spi_interp = Vnml_spi_interp(:,1);
    VL_spc_interp = Vnml_spc_interp(:,3);
    VM_spc_interp = Vnml_spc_interp(:,2);
    VN_spc_interp = Vnml_spc_interp(:,1);
    % calculate Pk, Pm and Pt
    k_B = 1.38064852e-23;
    mu_0 = 1.2566370614e-6;
    Pk_spi_interp = Np_spi_interp * 1e6 .* Tp_spi_interp * 11605 * k_B * 1e9; % unit:nPa
    Pm_spi_interp = (Br_interp.^2 + Bt_interp.^2 + Bn_interp.^2) * 1e-18 / 2 / mu_0 * 1e9; % unit:nPa
    Pt_spi_interp = Pk_spi_interp + Pm_spi_interp; % unit:nPa
    Pk_spc_interp = Np_spc_interp * 1e6 .* Tp_spc_interp * 11605 * k_B * 1e9; % unit:nPa
    Pm_spc_interp = (Br_interp.^2 + Bt_interp.^2 + Bn_interp.^2) * 1e-18 / 2 / mu_0 * 1e9; % unit:nPa
    Pt_spc_interp = Pk_spc_interp + Pm_spc_interp; % unit:nPa
    % calculate E_LMN
    Enml_spi_interp = zeros(size(Bnml_interp));
    Enml_spc_interp = zeros(size(Bnml_interp));
    E_spi_mod_interp = zeros(1,nTime);
    E_spc_mod_interp = zeros(1,nTime);
    for i_time = 1 : nTime
        Enml_spi_interp(i_time,:) = cross(-Vnml_spi_interp(i_time,:),Bnml_interp(i_time,:)); % unit: uV/m
        E_spi_mod_interp(i_time) = sqrt(dot(Enml_spi_interp(i_time,:),Enml_spi_interp(i_time,:)));
        Enml_spc_interp(i_time,:) = cross(-Vnml_spc_interp(i_time,:),Bnml_interp(i_time,:)); % unit: uV/m
        E_spc_mod_interp(i_time) = sqrt(dot(Enml_spc_interp(i_time,:),Enml_spc_interp(i_time,:)));
    end
    EL_spi_interp = Enml_spi_interp(:,3);
    EM_spi_interp = Enml_spi_interp(:,2);
    EN_spi_interp = Enml_spi_interp(:,1);
    EL_spc_interp = Enml_spc_interp(:,3);
    EM_spc_interp = Enml_spc_interp(:,2);
    EN_spc_interp = Enml_spc_interp(:,1);
    % calculate Alfven speed
    mu_0 = 1.2566370614e-6;
    mp = 1.67262192e-27; % unit: kg
    Va_spi_interp = B_mod_interp ./ sqrt(mu_0 .* mp .* Np_spi_interp) * 1e-15; % unit: km/s
    Va_spc_interp = B_mod_interp ./ sqrt(mu_0 .* mp .* Np_spc_interp) * 1e-15; % unit: km/s
    Ea_spi_interp = Va_spi_interp .* B_mod_interp; % unit: uV/m
    Ea_spc_interp = Va_spc_interp .* B_mod_interp; % unit: uV/m
    %% calculate reconnection rate
%     time_lead_beg = datenum(lead_beg_lst(which_cross,:));
%     time_lead_end = datenum(lead_end_lst(which_cross,:));
%     time_tral_beg = datenum(tral_beg_lst(which_cross,:));
%     time_tral_end = datenum(tral_end_lst(which_cross,:));
%     lead_index = find(std_Epoch >= time_lead_beg & std_Epoch <= time_lead_end);
%     tral_index = find(std_Epoch >= time_tral_beg & std_Epoch <= time_tral_end);
%     lead_Epoch = std_Epoch(lead_index);
%     tral_Epoch = std_Epoch(tral_index);
%     
%     VN_spi_lead = VN_spi_interp(lead_index);
%     Va_spi_lead = Va_spi_interp(lead_index);
%     VN_spi_tral = VN_spi_interp(tral_index);
%     Va_spi_tral = Va_spi_interp(tral_index);
%     VN_inflow_spi = (nanmean(VN_spi_tral) - nanmean(VN_spi_lead)) / 2;
%     Va_spi = nanmean([Va_spi_lead, Va_spi_tral]);
%     Ma_spi = abs(VN_inflow_spi / Va_spi);
%     
%     VN_spc_lead = VN_spc_interp(lead_index);
%     Va_spc_lead = Va_spc_interp(lead_index);
%     VN_spc_tral = VN_spc_interp(tral_index);
%     Va_spc_tral = Va_spc_interp(tral_index);
%     VN_inflow_spc = (mean(VN_spc_tral) - mean(VN_spc_lead)) / 2;
%     Va_spc = mean([Va_spc_lead, Va_spc_tral]);
%     Ma_spc = abs(VN_inflow_spc / Va_spc);
    %% calculate shreshold of K-H instability
%     BL_lead = BL_interp(lead_index);
%     BL_tral = BL_interp(tral_index);
%     BM_lead = BM_interp(lead_index);
%     BM_tral = BM_interp(tral_index);
%     B_mod_lead = B_mod_interp(lead_index);
%     B_mod_tral = B_mod_interp(tral_index);
%     
%     VM_spi_lead = VM_spi_interp(lead_index);
%     VM_spi_tral = VM_spi_interp(tral_index);
%     Np_spi_lead = Np_spi_interp(lead_index);
%     Np_spi_tral = Np_spi_interp(tral_index);
%     rho_spi_lead = mp .* Np_spi_lead; % [kg/cm3]
%     rho_spi_tral = mp .* Np_spi_tral; % [kg/cm3]
%     
%     VM_spc_lead = VM_spc_interp(lead_index);
%     VM_spc_tral = VM_spc_interp(tral_index);
%     Np_spc_lead = Np_spc_interp(lead_index);
%     Np_spc_tral = Np_spc_interp(tral_index);
%     rho_spc_lead = mp .* Np_spc_lead; % [kg/cm3]
%     rho_spc_tral = mp .* Np_spc_tral; % [kg/cm3]
%     % get threshold by finding the min-value of function    
%     VM_spi_thres_fun = @(ang_V_shear_kt) 1e-15.*sqrt((1./mean(rho_spi_lead)+1./mean(rho_spi_lead))./mu_0./cos(ang_V_shear_kt).^2.* ...
%         ((mean(BM_lead).*cos(ang_V_shear_kt)+mean(BL_lead).*sin(ang_V_shear_kt)).^2+ ...
%          (mean(BM_tral).*cos(ang_V_shear_kt)+mean(BL_tral).*sin(ang_V_shear_kt)).^2)); % [km/s]
%     ang_spi_thres = fminbnd(VM_spi_thres_fun,0,2*pi);
%     VM_spi_thres = VM_spi_thres_fun(ang_spi_thres);
% 
%     VM_spc_thres_fun = @(ang_V_shear_kt) 1e-15.*sqrt((1./mean(rho_spc_lead)+1./mean(rho_spc_lead))./mu_0./cos(ang_V_shear_kt).^2).* ...
%         ((mean(BM_lead).*cos(ang_V_shear_kt)+mean(BL_lead).*sin(ang_V_shear_kt)).^2+ ...
%          (mean(BM_tral).*cos(ang_V_shear_kt)+mean(BL_tral).*sin(ang_V_shear_kt)).^2); % [km/s]
%     ang_spc_thres = fminbnd(VM_spc_thres_fun,0,2*pi);
%     VM_spc_thres = VM_spc_thres_fun(ang_spc_thres);
%     % shear velocity in fact
%     VM_shear_spi = abs(mean(VM_spi_tral) - mean(VM_spi_lead));
%     VM_shear_spc = abs(mean(VM_spc_tral) - mean(VM_spc_lead));
    %% plot figures: Brtn, Vrtn, Np&Tp, Pmkt, PA
    fig_num = 5;
    LineWidth = 1.5; FontSize = 10;
    w_base = 0.15; h_base  = 0.15;
    height = 0.12; width = 0.7; space = 0.002;
    subplot('position',[w_base,h_base+height*0+space*0,width,height])
    
%     subplot(fig_num,1,1)
    subplot('position',[w_base,h_base+height*4+space*4,width,height])
    plot(std_Epoch,Br_interp,'r','LineWidth',LineWidth); hold on
    plot(std_Epoch,Bt_interp,'g','LineWidth',LineWidth); hold on
    plot(std_Epoch,Bn_interp,'b','LineWidth',LineWidth); hold on
    plot(std_Epoch,B_mod_interp,'k','LineWidth',LineWidth); grid on
    ylabel('Brtn [nT]')
    for i_cross = 1 : cross_num(i_flap)
        mark_beg = datenum(mark_beg_lst(i_cross,:));
        mark_end = datenum(mark_end_lst(i_cross,:));
        shadow_marks(mark_beg,mark_end);
    end
    legend('Br','Bt','Bn','|B|','Location','eastoutside')
    xlim([time_plot_beg time_plot_end]);
    datetick('x','HH:MM');
    set(gca,'xticklabel',[],'tickdir','in','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
 
    if spc_or_spi(i_flap,:) == 'c'        
%         subplot(fig_num,1,2)
%         plot(std_Epoch,Vr_spc_interp-mean(Vr_spc_interp),'r','LineWidth',LineWidth); hold on
%         plot(std_Epoch,Vt_spc_interp-mean(Vt_spc_interp),'g','LineWidth',LineWidth); hold on
%         plot(std_Epoch,Vn_spc_interp-mean(Vn_spc_interp),'b','LineWidth',LineWidth); grid on
%         ylabel('Vrtn(spc) [km/s]')
%         legend('Vr','Vt','Vn','Location','eastoutside')
%         xlim([time_plot_beg time_plot_end]);
%         datetick('x','HH:MM');
%         set(gca,'tickdir','in','LineWidth',LineWidth,'FontSize',FontSize);
        
%         subplot(fig_num,1,3)
%         yyaxis left
%         plot(std_Epoch,Np_spc_interp,'r','LineWidth',LineWidth);
%         ax = gca; ax.YColor = 'r';
%         %     ylabel('N_p [cm^{-3}]')
%         yyaxis right
%         plot(std_Epoch,Tp_spc_interp,'b','LineWidth',LineWidth);
%         ax = gca; ax.YColor = 'b';
%         %     ylabel('T_p [eV]')
%         grid on
%         legend('N_p','T_p','Location','eastoutside')
%         xlim([time_plot_beg time_plot_end]);
%         datetick('x','HH:MM');
%         set(gca,'tickdir','in','LineWidth',LineWidth,'FontSize',FontSize);
      
        subplot(fig_num,1,3)
        plot(std_Epoch,VL_spc_interp-mean(VL_spc_interp),'r','LineWidth',LineWidth); hold on
        plot(std_Epoch,VM_spc_interp-mean(VM_spc_interp),'g','LineWidth',LineWidth); hold on
        plot(std_Epoch,VN_spc_interp-mean(VN_spc_interp),'b','LineWidth',LineWidth); grid on
        ylabel('V_{LMN} [km/s]')
        for i_cross = 1 : cross_num(i_flap)
            mark_beg = datenum(mark_beg_lst(i_cross,:));
            mark_end = datenum(mark_end_lst(i_cross,:));
            shadow_marks(mark_beg,mark_end);
        end
        legend('V_L','V_M','V_N','Location','eastoutside')
        xlim([time_plot_beg time_plot_end]);
        datetick('x','HH:MM');
        set(gca,'tickdir','in','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
        
        subplot(fig_num,1,4)
        plot(std_Epoch,EL_spc_interp,'r','LineWidth',LineWidth); hold on
        plot(std_Epoch,EM_spc_interp,'g','LineWidth',LineWidth); hold on
        plot(std_Epoch,EN_spc_interp,'b','LineWidth',LineWidth); hold on
        plot(std_Epoch,E_spc_mod_interp,'k','LineWidth',LineWidth); grid on
        ylabel('E_{LMN} [uV/m]')
        for i_cross = 1 : cross_num(i_flap)
            mark_beg = datenum(mark_beg_lst(i_cross,:));
            mark_end = datenum(mark_end_lst(i_cross,:));
            shadow_marks(mark_beg,mark_end);
        end
        legend('E_L','E_M','E_N','|E|','Location','eastoutside')
        xlim([time_plot_beg time_plot_end]);
        datetick('x','HH:MM');
        set(gca,'tickdir','in','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
        
        subplot(fig_num,1,5)
        plot(std_Epoch,Mr_spc_interp,'k','LineWidth',LineWidth); grid on
        ylabel('M_R')
        for i_cross = 1 : cross_num(i_flap)
            mark_beg = datenum(mark_beg_lst(i_cross,:));
            mark_end = datenum(mark_end_lst(i_cross,:));
            shadow_marks(mark_beg,mark_end);
        end
        legend('M_R','Location','eastoutside')
        xlim([time_plot_beg time_plot_end]);
        ylim([0 0.4]);
        datetick('x','HH:MM');
        set(gca,'tickdir','in','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
        
%         subplot(fig_num,1,4)
%         plot(std_Epoch,Pk_spc_interp,'r','LineWidth',LineWidth); hold on
%         plot(std_Epoch,Pm_spc_interp,'b','LineWidth',LineWidth); hold on
%         plot(std_Epoch,Pt_spc_interp,'k','LineWidth',LineWidth); grid on
%         ylabel('P [nPa]')
%         legend('P_k','P_m','P_t','Location','eastoutside')
%         xlim([time_plot_beg time_plot_end]);
%         datetick('x','HH:MM');
%         set(gca,'tickdir','in','LineWidth',LineWidth,'FontSize',FontSize);
    end
     
    if spc_or_spi(i_flap,:) == 'i'        
%         subplot(fig_num,1,2)
%         plot(std_Epoch,Vr_spi_interp-mean(Vr_spi_interp),'r','LineWidth',LineWidth); hold on
%         plot(std_Epoch,Vt_spi_interp-mean(Vt_spi_interp),'g','LineWidth',LineWidth); hold on
%         plot(std_Epoch,Vn_spi_interp-mean(Vn_spi_interp),'b','LineWidth',LineWidth); grid on
%         ylabel('Vrtn(spi) [km/s]')
%         legend('Vr','Vt','Vn','Location','eastoutside')
%         xlim([time_plot_beg time_plot_end]);
%         datetick('x','HH:MM');
%         set(gca,'tickdir','in','LineWidth',LineWidth,'FontSize',FontSize);
        
%         subplot(fig_num,1,4)
        subplot('position',[w_base,h_base+height*1+space*1,width,height])
        yyaxis left
        plot(std_Epoch,Np_spi_interp,'r','LineWidth',LineWidth);
        ax = gca; ax.YColor = 'r';
        ylabel('N_p [cm^{-3}]')
        for i_cross = 1 : cross_num(i_flap)
            mark_beg = datenum(mark_beg_lst(i_cross,:));
            mark_end = datenum(mark_end_lst(i_cross,:));
            shadow_marks(mark_beg,mark_end);
        end
        yyaxis right
        plot(std_Epoch,Tp_spi_interp,'b','LineWidth',LineWidth);
        ax = gca; ax.YColor = 'b';
        ylabel('T_p [eV]')
        grid on
%         legend('N_p','T_p','Location','eastoutside')
        xlim([time_plot_beg time_plot_end]);
        datetick('x','HH:MM');
        set(gca,'xticklabel',[],'tickdir','in','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
        
%         subplot(fig_num,1,3)
        subplot('position',[w_base,h_base+height*2+space*2,width,height])
        plot(std_Epoch,VL_spi_interp-mean(VL_spi_interp),'r','LineWidth',LineWidth); hold on
        plot(std_Epoch,VM_spi_interp-mean(VM_spi_interp),'g','LineWidth',LineWidth); hold on
        plot(std_Epoch,VN_spi_interp-mean(VN_spi_interp),'b','LineWidth',LineWidth); grid on
        ylabel('V_{LMN} [km/s]')
        for i_cross = 1 : cross_num(i_flap)
            mark_beg = datenum(mark_beg_lst(i_cross,:));
            mark_end = datenum(mark_end_lst(i_cross,:));
            shadow_marks(mark_beg,mark_end);
        end
        legend('V_L','V_M','V_N','Location','eastoutside')
        xlim([time_plot_beg time_plot_end]);
        datetick('x','HH:MM');
        set(gca,'xticklabel',[],'tickdir','in','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
        
%         subplot(fig_num,1,4)
%         plot(std_Epoch,EL_spi_interp,'r','LineWidth',LineWidth); hold on
%         plot(std_Epoch,EM_spi_interp,'g','LineWidth',LineWidth); hold on
%         plot(std_Epoch,EN_spi_interp,'b','LineWidth',LineWidth); hold on
%         plot(std_Epoch,E_spi_mod_interp,'k','LineWidth',LineWidth); grid on
%         ylabel('E_{LMN} [uV/m]')
%         legend('E_L','E_M','E_N','|E|','Location','eastoutside')
%         xlim([time_plot_beg time_plot_end]);
%         datetick('x','HH:MM');
%         set(gca,'tickdir','in','LineWidth',LineWidth,'FontSize',FontSize);
        
%         subplot(fig_num,1,5)
%         VN_spi_pos = VN_spi_interp - mean(VN_spi_interp);
%         VN_spi_pos(VN_spi_pos<0) = 0;
%         VN_spi_neg = VN_spi_interp - mean(VN_spi_interp);
%         VN_spi_neg(VN_spi_neg>0) = 0;
%         area(std_Epoch,VN_spi_pos,'FaceColor','r'); hold on;
%         area(std_Epoch,VN_spi_neg,'FaceColor','b'); grid on
%         ylabel('V_N')
%         legend('V_N^{>0}','V_N^{<0}','Location','eastoutside')
%         xlim([time_plot_beg time_plot_end]);
% %         ylim([-50 50]);
%         datetick('x','HH:MM');
%         set(gca,'tickdir','in','LineWidth',LineWidth,'FontSize',FontSize);
        
%         subplot(fig_num,1,4)
%         plot(std_Epoch,Pk_spi_interp,'r','LineWidth',LineWidth); hold on
%         plot(std_Epoch,Pm_spi_interp,'b','LineWidth',LineWidth); hold on
%         plot(std_Epoch,Pt_spi_interp,'k','LineWidth',LineWidth); grid on
%         ylabel('P [nPa]')
%         legend('P_k','P_m','P_t','Location','eastoutside')
%         xlim([time_plot_beg time_plot_end]);
%         datetick('x','HH:MM');
%         set(gca,'tickdir','in','LineWidth',LineWidth,'FontSize',FontSize);
    end
    
%     subplot(fig_num,1,2)
    subplot('position',[w_base,h_base+height*3+space*3,width,height])
    plot(std_Epoch,BL_interp,'r','LineWidth',LineWidth); hold on
    plot(std_Epoch,BM_interp,'g','LineWidth',LineWidth); hold on
    plot(std_Epoch,BN_interp,'b','LineWidth',LineWidth); grid on
    ylabel('B_{LMN} [km/s]')
    for i_cross = 1 : cross_num(i_flap)
        mark_beg = datenum(mark_beg_lst(i_cross,:));
        mark_end = datenum(mark_end_lst(i_cross,:));
        shadow_marks(mark_beg,mark_end);
    end
    legend('B_L','B_M','B_N','Location','eastoutside')
    xlim([time_plot_beg time_plot_end]);
    datetick('x','HH:MM');
    set(gca,'xticklabel',[],'tickdir','in','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
    
%     subplot(fig_num,1,5)
    subplot('position',[w_base,h_base+height*0+space*0,width,height])
    spe_x_plot = spe_Epoch_plot; spe_y_plot = PA_bins;
    [spe_X,spe_Y] = meshgrid(spe_x_plot,spe_y_plot);
    h = pcolor(spe_X,spe_Y,log10(PA_eflux_plot));
    set(h,'Linestyle','none');
    xlabel(['Time [HH:MM]']); ylabel({'314 eV'; 'e^- PA [deg.]'});
    for i_cross = 1 : cross_num(i_flap)
        mark_beg = datenum(mark_beg_lst(i_cross,:));
        mark_end = datenum(mark_end_lst(i_cross,:));
        vertical_marks(mark_beg,mark_end);
    end
%     legend('314','Location','eastoutside')
    xlim([time_plot_beg time_plot_end]);
    datetick('x','HH:MM');
    cb = colorbar; colormap jet;
    cb.Label.String = 'eV/cm2-s-ster-eV';
    colorbar('YTick',[8,9,10],'YTickLabel',{'10^8','10^9','10^{10}'});
    set(gca,'tickdir','out','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
        
    sgtitle(['PSP enounter ',num2str(encounter),' ',plot_beg_lst(i_flap,1:4),plot_beg_lst(i_flap,6:7),plot_beg_lst(i_flap,9:10),' ',plot_beg_lst(i_flap,12:19),'-',plot_end_lst(i_flap,12:19)]);
    if save_or_not == 1
        saveas(gcf,['C:\Users\22242\Documents\STUDY\Work\current_sheet_flapping\Encounter ',num2str(encounter),'\flapping_events\',plot_beg_lst(i_flap,1:4),plot_beg_lst(i_flap,6:7),plot_beg_lst(i_flap,9:10),'_',plot_beg_lst(i_flap,12:13),plot_beg_lst(i_flap,15:16),plot_beg_lst(i_flap,18:19),'-',plot_end_lst(i_flap,12:13),plot_end_lst(i_flap,15:16),plot_end_lst(i_flap,18:19),'.png']);
    end
end
%% functions
function [V,D,ratio21,ratio32] = get_LMN(Br,Bt,Bn,fld_Epoch,time_calc_beg,time_calc_end)
    Br_calc = Br(find(fld_Epoch >= time_calc_beg & fld_Epoch <= time_calc_end)); Br_calc(isnan(Br_calc)) = mean(Br_calc,'omitnan');
    Bt_calc = Bt(find(fld_Epoch >= time_calc_beg & fld_Epoch <= time_calc_end)); Bt_calc(isnan(Bt_calc)) = mean(Bt_calc,'omitnan');
    Bn_calc = Bn(find(fld_Epoch >= time_calc_beg & fld_Epoch <= time_calc_end)); Bn_calc(isnan(Bn_calc)) = mean(Bn_calc,'omitnan');
    M = [mean(Br_calc.*Br_calc) - mean(Br_calc)*mean(Br_calc), mean(Br_calc.*Bt_calc) - mean(Br_calc)*mean(Bt_calc), mean(Br_calc.*Bn_calc) - mean(Br_calc)*mean(Bn_calc);
         mean(Bt_calc.*Br_calc) - mean(Bt_calc)*mean(Br_calc), mean(Bt_calc.*Bt_calc) - mean(Bt_calc)*mean(Bt_calc), mean(Bt_calc.*Bn_calc) - mean(Bt_calc)*mean(Bn_calc);
         mean(Bn_calc.*Br_calc) - mean(Bn_calc)*mean(Br_calc), mean(Bn_calc.*Bt_calc) - mean(Bn_calc)*mean(Bt_calc), mean(Bn_calc.*Bn_calc) - mean(Bn_calc)*mean(Bn_calc)];
    [V, D] = eig(M);
    lambda1 = D(1,1);lambda2 = D(2,2); lambda3 = D(3,3); % lambda1 < lambda2 < lambda3; order: N, M, L
    ratio21 = lambda2 / lambda1; ratio32 = lambda3 / lambda2; % ratio21 usually should be limited
end
function [x_RTN,y_RTN,z_RTN] = calc_HCI2SCRTN(x_HCI,y_HCI,z_HCI,SC_HCIx,SC_HCIy,SC_HCIz)
%Change the coordiantes xyz in HCI frame to xyz in spacecraft RTN frame
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
function vertical_marks(epoch_beg,epoch_end)
% Plot vertical lines
    LineWidth = 2;
    hold on;
    xline(epoch_beg,'--k','LineWidth',LineWidth); hold on;
    xline(epoch_end,'--k','LineWidth',LineWidth); hold on;
end
function shadow_marks(epoch_beg,epoch_end)
% Plot shadow blocks
    hold on
    yl = ylim;
    vert = [epoch_beg,yl(1);epoch_end,yl(1);epoch_end,yl(2);epoch_beg,yl(2)];
    face = [1 2 3 4];
    patch('Faces',face,'Vertices',vert,'FaceColor','black','FaceAlpha',0.2,'EdgeColor','none');
    hold on
end