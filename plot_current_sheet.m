clear; close all;
cfg = hcs_config();
%% parameter
% unit conversion
deg2rad = pi/180; % from [deg.] to [rad.]
day2sec = 24*60*60; % from [day] to [sec]
km2RS = 1/696340; % from [km] to [Rs]
% scientific parameter
kB = 1.38064852e-23; % [m^2*kg*s^-2*K^-1]
G = 6.67259e-11; % [N*m^2*kg^-2]
Ms = 1.989e30; % [kg]
Rs = 6.9634e8; % [m]
mp = 1.6726219e-27; % [kg]
T = 1e6; % [K]
% critical value
ss = 2.5; % source surface [Rs]
vc = sqrt(2 * kB * T / mp); % [m/s]
rc = G * Ms * mp / 4 / kB / T; % [m]
%% import data: Br in the source surface
csv_store = [fullfile(cfg.legacy_work_root, '20210117'), filesep];
lon = importdata([csv_store,'longitude.csv']);
lon_arr = lon(:,1)'; % [deg.]
lat = importdata([csv_store,'latitude.csv']);
lat_arr = lat(1,:) - 90; % [deg.]
Br = importdata([csv_store,'Br.csv']);
Br = rot90(Br); % row: same longitude(phi), column: same latitude(theta)
%% standard grids
% angular grid
phi_num = 400;
phi_gap = 2 * pi / phi_num; % [rad.]
phi_std = 0 : phi_gap : 2*pi;
% radial grid
r_num = 20;
r_max = 22 * Rs; r_min = 2.5 * Rs;
r_gap = (r_max - r_min) / r_num; % unit: m
r_std = r_min : r_gap : r_max;
% modification range
modif_beg_ind = floor((phi_num + 1) / 3.6);
modif_end_ind = floor((phi_num + 1) / 3.5);
%% get original neutral line
% find reversing locations of Br
neut_arr = zeros(3,length(lon_arr));
count_arr = zeros(1,length(lon_arr));
for i_lon = 1 : length(lon_arr)
    for i_lat = 1 : length(lat_arr) - 1
        if Br(i_lat,i_lon) * Br(i_lat + 1,i_lon) < 0
            count_arr(i_lon) = count_arr(i_lon) + 1;
            neut_arr(count_arr(i_lon),i_lon) = lat_arr(i_lat);
        end
    end
end
% divide neutral line into pieces
piece_num = 5; piece_count = 0;
piece_beg = zeros(1,piece_num); piece_beg(1) = 1;
piece_end = zeros(1,piece_num); piece_end(end) = length(lon_arr);
for i_lon = 2 : length(lon_arr)
    if count_arr(i_lon) ~= count_arr(i_lon - 1)
        piece_count = piece_count + 1;
        piece_end(piece_count) = i_lon - 1;
        piece_beg(piece_count + 1) = i_lon;
    end
end
% divide into 5 major-pieces
neut_lon1 = lon_arr(piece_beg(1):piece_end(2)); % [deg.]
neut_lat1 = neut_arr(1,piece_beg(1):piece_end(2)); % [deg.]
neut_lon2 = lon_arr(piece_beg(2):piece_end(2));
neut_lat2 = neut_arr(2,piece_beg(2):piece_end(2));
neut_lon3 = lon_arr(piece_beg(2):piece_end(4));
neut_lat3 = [neut_arr(3,piece_beg(2):piece_end(2)),neut_arr(1,piece_beg(3):piece_end(4))];
neut_lon4 = lon_arr(piece_beg(4):piece_end(4));
neut_lat4 = neut_arr(2,piece_beg(4):piece_end(4));
neut_lon5 = lon_arr(piece_beg(4):piece_end(5));
neut_lat5 = [neut_arr(3,piece_beg(4):piece_end(4)),neut_arr(1,piece_beg(5):piece_end(5))];
%% smooth neutral line
neut_lon_arr = [neut_lon1,neut_lon2(end:-1:1),neut_lon3,neut_lon4(end:-1:1),neut_lon5];
neut_lat_arr = [neut_lat1,neut_lat2(end:-1:1),neut_lat3,neut_lat4(end:-1:1),neut_lat5];
neut_num_arr = 1 : length(neut_lon_arr);
neut_lon_sm = smooth(neut_lon_arr,15)';
neut_lat_sm = smooth(neut_lat_arr,15)';
% % check the neutral line in number
% subplot(2,2,1)
% plot(neut_num_arr,neut_lon_arr);
% grid on
% subtitle('longitude array to smooth')
% subplot(2,2,2)
% plot(neut_num_arr,neut_lat_arr);
% grid on
% subtitle('latitude array to smooth')
% subplot(2,2,3)
% plot(neut_num_arr,neut_lon_sm);
% grid on
% subtitle('longitude array smoothed')
% subplot(2,2,4)
% plot(neut_num_arr,neut_lat_sm);
% grid on
% subtitle('latitude array smoothed')

% smoothed pieces
neut_lon1_sm = neut_lon_sm(1:length(neut_lon1)); % [deg.]
neut_lat1_sm = neut_lat_sm(1:length(neut_lon1)); % [deg.]
neut_lon2_sm = neut_lon_sm(length(neut_lon1)+1:length(neut_lon1)+length(neut_lon2));
neut_lat2_sm = neut_lat_sm(length(neut_lon1)+1:length(neut_lon1)+length(neut_lon2));
neut_lon3_sm = neut_lon_sm(length(neut_lon1)+length(neut_lon2)+1:length(neut_lon1)+length(neut_lon2)+length(neut_lon3));
neut_lat3_sm = neut_lat_sm(length(neut_lon1)+length(neut_lon2)+1:length(neut_lon1)+length(neut_lon2)+length(neut_lon3));
neut_lon4_sm = neut_lon_sm(length(neut_lon1)+length(neut_lon2)+length(neut_lon3)+1:length(neut_lon1)+length(neut_lon2)+length(neut_lon3)+length(neut_lon4));
neut_lat4_sm = neut_lat_sm(length(neut_lon1)+length(neut_lon2)+length(neut_lon3)+1:length(neut_lon1)+length(neut_lon2)+length(neut_lon3)+length(neut_lon4));
neut_lon5_sm = neut_lon_sm(length(neut_lon1)+length(neut_lon2)+length(neut_lon3)+length(neut_lon4)+1:length(neut_lon1)+length(neut_lon2)+length(neut_lon3)+length(neut_lon4)+length(neut_lon5));
neut_lat5_sm = neut_lat_sm(length(neut_lon1)+length(neut_lon2)+length(neut_lon3)+length(neut_lon4)+1:length(neut_lon1)+length(neut_lon2)+length(neut_lon3)+length(neut_lon4)+length(neut_lon5));
%% interpolation into standard phi grid
% piece 1
phiS1 = interp1(neut_lon1_sm * deg2rad,neut_lon1_sm * deg2rad,phi_std,'linear'); % [rad.]
phiS1(isnan(phiS1)) = [];
thetaS1 = interp1(neut_lon1_sm * deg2rad,neut_lat1_sm * deg2rad,phi_std,'linear'); % [rad.]
thetaS1(isnan(thetaS1)) = [];
% piece 1, sub-piece 1: former
phiS1_form = phiS1(1:modif_beg_ind);
thetaS1_form = thetaS1(1:modif_beg_ind);
% piece 1, sub-piece 1: modify
phiS0 = phiS1(modif_beg_ind:modif_end_ind);
thetaS0 = thetaS1(modif_beg_ind:modif_end_ind);
% piece 1, sub-piece 1: later
phiS1_late = phiS1(modif_end_ind:end);
thetaS1_late = thetaS1(modif_end_ind:end);
% piece 2
phiS2 = interp1(neut_lon2_sm * deg2rad,neut_lon2_sm * deg2rad,phi_std,'linear');
phiS2(isnan(phiS2)) = [];
thetaS2 = interp1(neut_lon2_sm * deg2rad,neut_lat2_sm * deg2rad,phi_std,'linear');
thetaS2(isnan(thetaS2)) = [];
% piece 3
phiS3 = interp1(neut_lon3_sm * deg2rad,neut_lon3_sm * deg2rad,phi_std,'linear');
phiS3(isnan(phiS3)) = [];
thetaS3 = interp1(neut_lon3_sm * deg2rad,neut_lat3_sm * deg2rad,phi_std,'linear');
thetaS3(isnan(thetaS3)) = [];
% piece 4
phiS4 = interp1(neut_lon4_sm * deg2rad,neut_lon4_sm * deg2rad,phi_std,'linear');
phiS4(isnan(phiS4)) = [];
thetaS4 = interp1(neut_lon4_sm * deg2rad,neut_lat4_sm * deg2rad,phi_std,'linear');
thetaS4(isnan(thetaS4)) = [];
% piece 5
phiS5 = interp1(neut_lon5_sm * deg2rad,neut_lon5_sm * deg2rad,phi_std,'linear');
phiS5(isnan(phiS5)) = [];
thetaS5 = interp1(neut_lon5_sm * deg2rad,neut_lat5_sm * deg2rad,phi_std,'linear');
thetaS5(isnan(thetaS5)) = [];
%% modification of neutral line
modif_num = 30;
modif_amp = 0.08;
phiS_modif = linspace(phiS0(1),phiS0(end),modif_num);
thetaS_modif = interp1(phiS0,thetaS0,phiS_modif,'linear');
for i_phi = 2 : modif_num
    thetaS_modif(i_phi) = thetaS_modif(i_phi) - modif_amp * sin(i_phi / modif_num * 3 * pi);
end
%% calculate Parker Spiral
% combine neutral lines together to make corners smooth
phiS_tot = [phiS1_form,phiS_modif,phiS1_late,flip(phiS2),phiS3,phiS4,phiS5];
thetaS_tot = [thetaS1_form,thetaS_modif,thetaS1_late,flip(thetaS2),thetaS3,thetaS4,thetaS5];
% calculat numerical points of Parker Spiral
[x_tot,y_tot,z_tot,t] = CalcParker(r_num,r_gap,phiS_tot,thetaS_tot,r_std,vc,rc);
%% import data: orbit of PSP
data_store = [fullfile(cfg.psp_data_root, 'Encounter 7'), filesep];
spc_dir = [data_store,'psp_swp_spc_l3i_20210117_v02.cdf'];
spc_info = spdfcdfinfo(spc_dir);
spc_Epoch = spdfcdfread(spc_dir,'Variables','Epoch');
sc_pos_HCI = spdfcdfread(spc_dir,'Variables','sc_pos_HCI'); % unit: km
carr_lat = spdfcdfread(spc_dir,'Variables','carr_latitude'); % unit: deg
carr_lon_jet = spdfcdfread(spc_dir,'Variables','carr_longitude'); % unit: deg
carr_lon = carr_lon_jet - GetOmega(0) * t;
% extract time period of PSP trajectory
index_beg = 1; index_end = length(spc_Epoch);
pos_beg = sc_pos_HCI(index_beg,:); pos_end = sc_pos_HCI(index_end,:);
lat_beg = carr_lat(index_beg); lat_end = carr_lat(index_end);
lon_beg = carr_lon(index_beg); lon_end = carr_lon(index_end);
r_beg = sqrt(pos_beg(1).^2 + pos_beg(2).^2 + pos_beg(3).^2) * km2RS; % unit: Rs
r_end = sqrt(pos_end(1).^2 + pos_end(2).^2 + pos_end(3).^2) * km2RS; % unit: Rs
% turn to 3-D Cartesian coordinate system
x_beg = r_beg * cosd(lat_beg) * cosd(lon_beg);
y_beg = r_beg * cosd(lat_beg) * sind(lon_beg);
z_beg = r_beg * sind(lat_beg);
x_end = r_end * cosd(lat_end) * cosd(lon_end);
y_end = r_end * cosd(lat_end) * sind(lon_end);
z_end = r_end * sind(lat_end);
%% import data: Br, Bt, Bp in 3-D spherical coordinate
hdf_store = [fullfile(cfg.psi_data_root, '20210117'), filesep];
Br_file = [hdf_store,'br002.h5'];
Bt_file = [hdf_store,'bt002.h5'];
Bp_file = [hdf_store,'bp002.h5'];
Br_info = h5info(Br_file);
% original grid
rr = h5read(Br_file,'/scales1'); % distance from SUN [Rs.]
tt = h5read(Br_file,'/scales2'); % carrington latitude  [rad.]
pp = h5read(Br_file,'/scales3'); % carrington longitude [rad.]
rr(end) = []; % to match the size of Brtp
% import and neaten Brtp
Br_3d = h5read(Br_file,'/datas');
Bt_3d = h5read(Bt_file,'/datas');
Bp_3d = h5read(Bp_file,'/datas');
Br_3d(end,:,:) = []; Bt_3d(:,end,:) = []; Bp_3d(:,:,end) = []; % into the same size
% select needed range
rr_index = find((rr>1.2)&(rr<2.5)); rr_sub = rr(rr_index);
tt_index = find((tt>0.3)&(tt<2.7)); tt_sub = tt(tt_index);
pp_index = 1:10:length(pp); pp_sub = pp(pp_index);
Br_sub = Br_3d(rr_index,tt_index,pp_index); 
Bt_sub = Bt_3d(rr_index,tt_index,pp_index);
Bp_sub = Bp_3d(rr_index,tt_index,pp_index);
% change grid and data into 3-D Cartesian coordinate
[t_sub,r_sub,p_sub] = meshgrid(tt_sub,rr_sub,pp_sub);
x_sub = r_sub .* sin(t_sub) .* cos(p_sub);
y_sub = r_sub .* sin(t_sub) .* sin(p_sub);
z_sub = r_sub .* cos(t_sub);
Bx_sub = Br_sub .* sin(t_sub) .* cos(p_sub) + Bt_sub .* cos(t_sub) .* cos(p_sub) - Bp_sub .* sin(p_sub);
By_sub = Br_sub .* sin(t_sub) .* sin(p_sub) + Bt_sub .* cos(t_sub) .* sin(p_sub) + Bp_sub .* cos(p_sub);
Bz_sub = Br_sub .* cos(t_sub) - Bt_sub .* sin(t_sub);
% standard 3-D grids
xx_interp = -2.5:0.5:2.5;
yy_interp = -2.5:0.5:2.5;
zz_interp = -2.5:0.5:2.5;
[x_interp,y_interp,z_interp] = meshgrid(xx_interp,yy_interp,zz_interp);
Bx_interp = griddata(x_sub,y_sub,z_sub,Bx_sub,x_interp,y_interp,z_interp,'linear');
By_interp = griddata(x_sub,y_sub,z_sub,By_sub,x_interp,y_interp,z_interp,'linear');
Bz_interp = griddata(x_sub,y_sub,z_sub,Bz_sub,x_interp,y_interp,z_interp,'linear');
% start locations of streamline
start_r = 1.8;
start_gap = 10;
start_t1 = phiS1(1:start_gap:end); start_p1 = thetaS1(1:start_gap:end);
[start_x1,start_y1,start_z1] = Sphere2Carte(start_r,start_p1,start_t1);
start_t2 = phiS2(1:5:end); start_p2 = thetaS2(1:5:end);
[start_x2,start_y2,start_z2] = Sphere2Carte(start_r,start_p2,start_t2);
start_t3 = phiS3(1:5:end); start_p3 = thetaS3(1:5:end);
[start_x3,start_y3,start_z3] = Sphere2Carte(start_r,start_p3,start_t3);
start_t4 = phiS4(1:5:end); start_p4 = thetaS4(1:5:end);
[start_x4,start_y4,start_z4] = Sphere2Carte(start_r,start_p4,start_t4);
start_t5 = phiS5(1:5:end); start_p5 = thetaS5(1:5:end);
[start_x5,start_y5,start_z5] = Sphere2Carte(start_r,start_p5,start_t5);
% combine start locations together
start_x_tot = [start_x1,start_x2,start_x3,start_x4,start_x5];
start_y_tot = [start_y1,start_y2,start_y3,start_y4,start_y5];
start_z_tot = [start_z1,start_z2,start_z3,start_z4,start_z5];
%% plot figures
LineWidth = 3;
FontSize = 12;
%% figure 1: magnetogram in source surface
figure(1)
[X,Y] = meshgrid(lon_arr,lat_arr);

% plot Br in 2-D graph
h = pcolor(X,Y,Br);
hold on
set(h,'Linestyle','none');

% plot divided neutral line
plot(phiS_tot/deg2rad,thetaS_tot/deg2rad,'k','LineWidth',LineWidth)

% set properties of figure 1
cbr = colorbar;
colormap jet;
cbr.Label.String = 'Br [nT]';
xlabel('Carrington Longitude'); ylabel('Carrington Latitude');
% title('Neutral line in source surface')
axis equal
xlim([0 360]); ylim([-90 90]);
set(gca,'LineWidth',LineWidth,'FontSize',FontSize,'tickdir','out');
%% figure 2: heliospheric current sheet
figure(2)

% plot current sheet
surf_tot = surf(x_tot/Rs,y_tot/Rs,z_tot/Rs);
% set topography
topomap_hcs = turbo(64);
shading interp
% set material
material shiny
hold on

% plot solar surface
[sphx,sphy,sphz] = sphere(150);
% import aia data to cover solar surface
aiaimg = imread(fullfile(cfg.legacy_work_root, '20210117', 'aia.png'));
aiaimg = flip(aiaimg,1);
% set sphere surface
surf_aia = surf(sphx.*0.5,sphy.*0.5,sphz.*0.5);
surf_aia.FaceColor = 'texturemap';
surf_aia.CData = aiaimg;
surf_aia.EdgeColor = 'none';
% set topography
topo_aia = double(aiaimg);
topo_aia = (topo_aia-min(topo_aia))./(max(topo_aia)-min(topo_aia))*255+256;
% set material
material shiny
hold on

% add light
lig = light;
lig.Color = [1 1 1];
lighting phong

% set total colorbar
cmp = colormap(topomap_hcs);
caxis([-10 10]);
hold on

% plot PSP orbit (approximate to a straight line)
plot3([x_beg,x_end],[y_beg,y_end],[z_beg,z_end],'r','LineWidth',LineWidth);
hold on

% plot a circle (as a reference substance)
ref = 15;
cir_alpha = -pi : pi/30 : pi;
cirx = ref * sin(cir_alpha);
ciry = ref * cos(cir_alpha);
cirz = zeros(1,length(cirx));
plot3(cirx,ciry,cirz,'k','LineWidth',LineWidth/2);
hold on

% plot closed magnetic field lines
line1 = streamline(x_interp,y_interp,z_interp,Bx_interp,By_interp,Bz_interp,start_x_tot,start_y_tot,start_z_tot);
hold on
line2 = streamline(x_interp,y_interp,z_interp,-Bx_interp,-By_interp,-Bz_interp,start_x_tot,start_y_tot,start_z_tot);
set(line1,'LineWidth',LineWidth/2,'Color','b');
set(line2,'LineWidth',LineWidth/2,'Color','b');

% set properties
box off;
axis off;
axis equal;
legend('Heliospheric current sheet','Photosphere','PSP orbit','Reference Circle of 15 R_s','Close magnetic field lines','Location','best');
set(gca,'LineWidth',LineWidth/2,'FontSize',FontSize);
%% function
function omega = GetOmega(carr_lat)
% Get rotation rate of the Sun
%   input: Carrington latitude in rad.
%   output: Omega in rad./s
    deg2rad = pi/180; % from degree to radian
    day2sec = 24*60*60; % from day to second

    omega1 = 14.713; omega2 = -2.396; omega3 = -1.787; % unit deg./day
    omega = (omega1 + omega2 * sin(carr_lat).^2 + omega3 * sin(carr_lat).^4) * deg2rad / day2sec;
end
function v = SolveEqu(r,vc,rc)
% Solve the Parker Solution of solar wind
%   input: r: distance from calculate point to Sun in m
%          vc: critical velocity in m/s
%          rc: critical distance in m;
%   output: v: velovity of Parker Solution
    func = @(v) (v/vc)^2-2*log(v/vc)-4*log(r/rc)-4*(rc/r)+3;
    options = optimoptions('fsolve','Display','none');
    if r<rc
%         syms vv
%         v = double(vpasolve((vv/vc)^2-2*log(vv/vc)-4*log(r/rc)-4*(rc/r)+3,vv,[0,0.999*vc])); % unit: m/s
        v = fsolve(func,0.5*vc,options); % unit: m/s
    else
%         syms vv
%         v = double(vpasolve((vv/vc)^2-2*log(vv/vc)-4*log(r/rc)-4*(rc/r)+3,vv,[1.001,3*vc])); % unit: m/s
        v = fsolve(func,1.5*vc,options); % unit: m/s
    end
end
function [x,y,z,t] = CalcParker(r_num,r_gap,phi_std,thetaS,r_arr,vc,rc)
% Calculate numerical points array of Parker Spiral
%   input: r_num: number of radial distance array to calculate
%          r_gap: gap of radial distance array to calculate in m
%          phi_std: Carrington latitude array in the source in rad.
%          thetaS: Carrington longitude array in the source surface in rad.
%          r_arr: radial distance array to calculate in m
%          vc: critical velocity in m/s
%          rc: critical distance in m;
%   output: [x,y,z]: Magnetic field lines origined from (phi_arr,theta_arr) in Cartesian coordinate system in m
%                 t: time needed to travel along Parker Solution in s
    phi_num = length(phi_std);
    x = zeros(phi_num,r_num);
    y = zeros(phi_num,r_num);
    z = zeros(phi_num,r_num);
    vel = zeros(1,r_num);
    dr = r_gap;
    for i_phi = 1 : phi_num
        Omega = GetOmega(thetaS(i_phi));
        phi = phi_std(i_phi);
        theta = thetaS(i_phi);
        t = 0;
        for i_r = 1 : r_num
            V = SolveEqu(r_arr(i_r),vc,rc);
            dphi = Omega * dr / V;
            phi = phi - dphi;
            dt = dr / V;
            t = t + dt;
            x(i_phi,i_r) = r_arr(i_r) * cos(theta) * cos(phi);
            y(i_phi,i_r) = r_arr(i_r) * cos(theta) * sin(phi);
            z(i_phi,i_r) = r_arr(i_r) * sin(theta);
            vel(i_r) = V;
        end
    end
end
function [x,y,z] = Sphere2Carte(r,lat,lon)
% Change vectors from sphere coordinate into Cartesian coordinate
% input:  (r,lat,lon) in ([UNIT],[rad.],[rad.])
% output:     (x,y,z) in ([UNIT],[UNIT],[UNIT])
    x = r .* cos(lat) .* cos(lon);
    y = r .* cos(lat) .* sin(lon);
    z = r .* sin(lat);
end
