clear; close all;
cfg = hcs_config();
figure();

% photosphere image by AIA
img = imread(fullfile(cfg.work_root, 'Encounter 7', 'flapping_events', ...
    '20210117', 'aia.png'));
% footpoint by backmapping
lon = 101.42851491767911;
lat = -27.155166512345676;

% range for plot 
min_x = 0;
max_x = 360;
min_y = -90;
max_y = 90;
% plot photosphere
imagesc([min_x max_x], [min_y max_y], flip(img,1)); hold on;
% plot footpoint
scatter(lon,lat,100,'filled','b');

% figure setting
axis equal
legend('footpoint')
xlim([min_x max_x]); xlabel('Longitude [deg.]');
ylim([min_y max_y]); ylabel('Latitude [deg.]');
set(gca,'YDir','normal','LineWidth',2,'FontSize',15,'TickDir','out');
