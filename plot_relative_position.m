clear; close all;
cfg = hcs_config();
plot_time = '2021-01-17 12:00:00';
Rs = 6.955e5; % [km]
%% Earth position in HCI coordinate, obtained from get_Earth_position.py
pos_earth_HCI = [41.31776271, -4.80444563, 0.98379865]; % (lon, lat, distance) in (deg, deg, AU)
earth_HCI_x = cosd(pos_earth_HCI(1));
earth_HCI_y = sind(pos_earth_HCI(1));
%% STEREO-A position in HCI coordinate, obtained from cdaweb
pos_stereo_HCI = [345, 1.8]; % (lon, lat) in (deg, deg)
stereo_HCI_x = cosd(pos_stereo_HCI(1));
stereo_HCI_y = sind(pos_stereo_HCI(1));
%% PSP position in HCI coordinate
% import data
spc_dir = [fullfile(cfg.psp_data_root, 'Encounter 7'), filesep];
spc_file = ['psp_swp_spc_l3i_',plot_time(1:4),plot_time(6:7),plot_time(9:10),'_v02.cdf'];
spc_path = [spc_dir,spc_file];
spc_info = spdfcdfinfo(spc_path);

% obtain PSP position
spc_Epoch = spdfcdfread(spc_path,'Variables','Epoch');
sc_pos_HCI = spdfcdfread(spc_path,'Variables','sc_pos_HCI');

num_epoch = length(spc_Epoch);
sc_pos_HCI = sc_pos_HCI / Rs;
pos_psp_HCI = sc_pos_HCI(num_epoch/2,:);
%% plot relative position
figure()
LineWidth = 2;
FontSize = 15;
sz = 50;
% Sun position
scatter(0,0,sz*10,'k','filled'); hold on
% Earth position
quiver(0,0,10*earth_HCI_x,10*earth_HCI_y,':r','LineWidth',LineWidth,'MaxHeadSize', 3, 'AutoScale', 'off'); hold on
plot([-20*earth_HCI_y,20*earth_HCI_y],[20*earth_HCI_x,-20*earth_HCI_x],'r','LineWidth',LineWidth); hold on
% STEREO position
quiver(0,0,10*stereo_HCI_x,10*stereo_HCI_y,':g','LineWidth',LineWidth,'MaxHeadSize', 3, 'AutoScale', 'off'); hold on
plot([-20*stereo_HCI_y,20*stereo_HCI_y],[20*stereo_HCI_x,-20*stereo_HCI_x],'g','LineWidth',LineWidth); hold on
% PSP position
scatter(pos_psp_HCI(1),pos_psp_HCI(2),sz,'b','filled'); hold on
plot(sc_pos_HCI(:,1),sc_pos_HCI(:,2),'b','LineWidth',LineWidth);
grid on

legend('SUN','Earth direction','L1 obs-plane','STEREO-A direction','STEREO-A obs-plane','PSP','Location','southwest');
axis equal
xlim([-25,25])
ylim([-25,25])
xlabel('HCI X [Rs]')
ylabel('HCI Y [Rs]')

annotation('textbox',[0.6,0,0.1,0.05],'String','plotted by plot\_relative\_position.m','FitBoxToText','on');

set(gca,'LineWidth',LineWidth,'FontSize',FontSize)
