clear; close all;
cfg = hcs_config();
img_dir = [fullfile(cfg.helioviewer_data_root, '20210117', '171images'), filesep];
img_num = 97;
img_height = 928;
img_width = 1536;
cnl_num = 3;
min_tot = 300;
min_beg = 143;
min_end = 157;
img_beg = fix((min_beg / min_tot) * img_num);
img_end = fix((min_end / min_tot) * img_num);
color_lst = ['r','g','b'];
%% region in the shape of rectangle
slit_width = 5;
slit_height = 100;
slit_num = 3;
% northwest vertex of each rectangle
vertx = [780,800,820];
verty = [460,475,490];
% index range of each rectangle
indx = zeros(slit_num,slit_width + 1);
indy = zeros(slit_num,slit_height + 1);
for i_slit = 1 : slit_num
    indx(i_slit,:) = vertx(i_slit) : vertx(i_slit) + slit_width;
    indy(i_slit,:) = verty(i_slit) : verty(i_slit) + slit_height;
end
%% extract selected region of each image
img_plot = zeros(slit_height+1,img_num,cnl_num,slit_num);
for i_img = 1 : img_num
    image_file = [num2str(i_img), '.png'];
    img_mat = imread([img_dir, image_file]);
    for i_slit = 1 : slit_num
        slitRGB = img_mat(indy(i_slit,:),indx(i_slit,:),:);
        img_plot(:,i_img,:,i_slit) = mean(slitRGB(:,:,:),2);
    end
end
%% saturation adjustment
sat = 50;
sat_plot = zeros(size(img_plot));
for i_slit = 1 : slit_num
    sat_plot(:,:,:,i_slit) = SaturationAdjustment(uint8(squeeze(img_plot(:,:,:,i_slit))),sat);
end
%% image enhancement
enh_imadjust = zeros(size(img_plot));
enh_histeq = zeros(size(img_plot));
enh_adapthisteq = zeros(size(img_plot));
for i_slit = 1 : slit_num
    [enh_imadjust(:,:,:,i_slit),enh_histeq(:,:,:,i_slit),enh_adapthisteq(:,:,:,i_slit)] = ImgEnhance(uint8(squeeze(img_plot(:,:,:,i_slit))));
end
%% saturation adjustment
sat_imadjust = zeros(size(img_plot));
sat_histeq = zeros(size(img_plot));
sat_adapthisteq = zeros(size(img_plot));
for i_slit = 1 : slit_num
    [sat_imadjust(:,:,:,i_slit),sat_histeq(:,:,:,i_slit),sat_adapthisteq(:,:,:,i_slit)] = ImgEnhance(uint8(squeeze(sat_plot(:,:,:,i_slit))));
end
%% plot figures
LineWidth = 1;
FontSize = 10;
xtick = [1;(i_img-1)/5;(i_img-1)*2/5;(i_img-1)*3/5;(i_img-1)*4/5;i_img];
xticklabel = ['01:00';'02:00';'03:00';'04:00';'05:00';'06:00'];
%% plot original image of each region
% figure()
% for i_slit = 1 : slit_num
%     subplot(3,2,i_slit);
%     image(uint8(squeeze(img_plot(:,:,:,i_slit))));
%     hold on
%     PlotPeriod(img_beg,img_end,slit_height);
%     set(gca,'xtick',xtick,'xticklabel',xticklabel,'tickdir','out','LineWidth',LineWidth,'FontSize',FontSize);
%     subtitle(['slit ',num2str(i_slit)],'FontSize',FontSize*1.6);
% end
% sgtitle('Original image without enhancement','FontSize',FontSize*2);
%% plot image enhance of each region
% % imadjust
% figure()
% for i_slit = 1 : slit_num
%     subplot(3,2,i_slit);
%     image(enh_imadjust(:,:,:,i_slit));
%     hold on
%     PlotPeriod(img_beg,img_end,slit_height);
%     set(gca,'xtick',xtick,'xticklabel',xticklabel,'tickdir','out','LineWidth',LineWidth,'FontSize',FontSize);
%     subtitle(['slit ',num2str(i_slit)],'FontSize',FontSize*1.6);
% end
% sgtitle('Enhancement method: imadjust','FontSize',FontSize*2);
% % histeq
% figure()
% for i_slit = 1 : slit_num
%     subplot(3,2,i_slit);
%     image(enh_histeq(:,:,:,i_slit));
%     hold on
%     PlotPeriod(img_beg,img_end,slit_height);
%     set(gca,'xtick',xtick,'xticklabel',xticklabel,'tickdir','out','LineWidth',LineWidth,'FontSize',FontSize);
%     subtitle(['slit ',num2str(i_slit)],'FontSize',FontSize*1.6);
% end
% sgtitle('Enhancement method: histeq','FontSize',FontSize*2);
% % adapthisteq
% figure()
% for i_slit = 1 : slit_num
%     subplot(3,2,i_slit);
%     image(enh_adapthisteq(:,:,:,i_slit));
%     hold on
%     PlotPeriod(img_beg,img_end,slit_height);
%     set(gca,'xtick',xtick,'xticklabel',xticklabel,'tickdir','out','LineWidth',LineWidth,'FontSize',FontSize);
%     subtitle(['slit ',num2str(i_slit)],'FontSize',FontSize*1.6);
% end
% sgtitle('Enhancement method: adapthisteq','FontSize',FontSize*2);
%% plot saturation adjustment for each region
% % origin + saturation
% figure()
% for i_slit = 1 : slit_num
%     subplot(3,2,i_slit);
%     image(uint8(squeeze(sat_plot(:,:,:,i_slit))));
%     hold on
%     PlotPeriod(img_beg,img_end,slit_height);
%     set(gca,'xtick',xtick,'xticklabel',xticklabel,'tickdir','out','LineWidth',LineWidth,'FontSize',FontSize);
%     subtitle(['slit ',num2str(i_slit)],'FontSize',FontSize*1.6);
% end
% sgtitle('origin + saturation adjust','FontSize',FontSize*2);
% imadjust + saturation
figure()
for i_slit = 1 : slit_num
    subplot(3,1,i_slit);
    image(sat_imadjust(:,:,:,i_slit));
    hold on
    PlotPeriod(img_beg,img_end,slit_height);
    set(gca,'xtick',xtick,'xticklabel',xticklabel,'tickdir','out','XColor',color_lst(i_slit),'YColor',color_lst(i_slit),'XMinorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
%     subtitle(['Slit ',num2str(i_slit)],'FontSize',FontSize*1.6);
end
% sgtitle('imadjust + saturation adjust','FontSize',FontSize*2);
% % histeq + saturation
% figure()
% for i_slit = 1 : slit_num
%     subplot(3,2,i_slit);
%     image(sat_histeq(:,:,:,i_slit));
%     hold on
%     PlotPeriod(img_beg,img_end,slit_height);
%     set(gca,'xtick',xtick,'xticklabel',xticklabel,'tickdir','out','LineWidth',LineWidth,'FontSize',FontSize);
%     subtitle(['slit ',num2str(i_slit)],'FontSize',FontSize*1.6);
% end
% sgtitle('histeq + saturation adjust','FontSize',FontSize*2);
% % adapthisteq + saturation
% figure()
% for i_slit = 1 : slit_num
%     subplot(3,2,i_slit);
%     image(sat_adapthisteq(:,:,:,i_slit));
%     hold on
%     PlotPeriod(img_beg,img_end,slit_height);
%     set(gca,'xtick',xtick,'xticklabel',xticklabel,'tickdir','out','LineWidth',LineWidth,'FontSize',FontSize);
%     subtitle(['slit ',num2str(i_slit)],'FontSize',FontSize*1.6);
% end
% sgtitle('adapthisteq + saturation adjust','FontSize',FontSize*2);
%% function
function PlotPeriod(time_beg,time_end,slit_height)
    plot([time_beg,time_beg],[1,slit_height+1],'w');
    hold on
    plot([time_end,time_end],[1,slit_height+1],'w');
end
function [shadow_imadjust,shadow_histeq,shadow_adapthisteq] = ImgEnhance(shadow)
    shadow_lab = rgb2lab(shadow);
    max_luminosity = 100;
    L = shadow_lab(:,:,1)/max_luminosity;
    % imadjust
    shadow_imadjust = shadow_lab;
    shadow_imadjust(:,:,1) = imadjust(L)*max_luminosity;
    shadow_imadjust = lab2rgb(shadow_imadjust);
    % histeq
    shadow_histeq = shadow_lab;
    shadow_histeq(:,:,1) = histeq(L)*max_luminosity;
    shadow_histeq = lab2rgb(shadow_histeq);
    % adapthisteq
    shadow_adapthisteq = shadow_lab;
    shadow_adapthisteq(:,:,1) = adapthisteq(L)*max_luminosity;
    shadow_adapthisteq = lab2rgb(shadow_adapthisteq);
end

function Image_new = SaturationAdjustment(src,saturation)
% input: src: RGB matrix of the original image
%        saturation: saturation adjustment increment range from -100 to 100
% output: Imgae_new: RGB matrix after saturation adjust
Image = src;
Image = double(Image);
R = Image(:,:,1);
G = Image(:,:,2);
B = Image(:,:,3);
[row, col] = size(R);
R_new = R;
G_new = G;
B_new = B;
% change to percentage for adjusting
Increment = saturation;
Increment = Increment / 100;
for i = 1 : row
    for j = 1 : col
        rgbMax = max(R(i,j),max(G(i,j),B(i,j)));
        rgbMin = min(R(i,j),min(G(i,j),B(i,j)));
        Delta = (rgbMax-rgbMin) / 255;
        if Delta == 0 % unable to adjust
            continue;
        end
        value = (rgbMax + rgbMin)/255;
        % calculate saturation from lightness
        L = value/2; 
        if L < 0.5
            S=Delta / value;
        else
            S =Delta / (2 - value);
        end
        if Increment >= 0
            if Increment + S >= 1
                alpha = S;
            else
                alpha = 1 - Increment;
            end
          alpha = 1 / alpha - 1;
          R_new(i,j) = R(i,j) + (R(i,j) - L * 255) * alpha;
          G_new(i,j) = G(i,j) + (G(i,j) - L * 255) * alpha;
          B_new(i,j) = B(i,j) + (B(i,j) - L * 255) * alpha;
        else
          alpha = Increment;
          R_new(i,j) = L*255 + (R(i,j) - L * 255) * (1 + alpha);
          G_new(i,j) = L*255 + (G(i,j) - L * 255) * (1 + alpha);
          B_new(i,j) = L*255 + (B(i,j) - L * 255) * (1 + alpha); 
        end
    end
end     
Image_new(:,:,1) = R_new;
Image_new(:,:,2) = G_new;
Image_new(:,:,3) = B_new;
end
