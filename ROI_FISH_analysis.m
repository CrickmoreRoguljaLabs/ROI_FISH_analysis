%% Parameters
% Intensitiy threshold for channel 1
int_thresh_Ch1 = 800;

% Area threshold
area_thresh_Ch1 = 15;

% Intensitiy threshold for channel 2
int_thresh_Ch2 = 300;

% Area threshold for channel 2
area_thresh_Ch2 = 15;

% Background size Ch1
Ch1_bg_size = 15;

% Background size Ch2
Ch2_bg_size = 5;

%% Loading and drawing ROIs
%
% Read image
fn = fullfile('D:\Dropbox\Crickmore_research\Images\in situ\batch 3 creb\batch 3 creb\MAX_PB.tif');
im = stkread2(fn);

% Get ROI
ROI = getpoly(im(:,:,1), 'Draw ROI');
%}
%% Thresholding and watershedding
% Apply polygon
im_P = im(:,:,1);
im_P(~ROI) = 0;

% remove background for Ch1
CH1_b = imopen(im_P,strel('disk', Ch1_bg_size));
im_P2 = im_P - CH1_b;

% Apply intensity threshold
im_t = im_P2 >= int_thresh_Ch1;

% Standard watershedding steps
im_d = -bwdist(~im_t);
im_d(~im_t) = Inf;
im_L = watershed(im_d);
im_L(~im_t) = 0;

% Apply area threshold
im_L2 = areathresh(im_L,0,area_thresh_Ch1,3,8) > 0;

%% Labeling and checking channel 2
% Label Ch1
im_labeled_Ch1 = bwlabel(im_L2);

% Number of punctas in Ch1
N_pun_Ch1 = max(im_labeled_Ch1(:));

% Checking if the channels are positive
Ch_pos_1 = nan(N_pun_Ch1, 2);
Ch_val_1 = nan(N_pun_Ch1, 2);

% Start processing Ch2
im2 = im(:,:,2);
im2(~ROI) = 0;

% Removing Ch2 background
im2_b = imopen(im2,strel('disk', Ch2_bg_size));
im2_nobg = im2 - im2_b;

% Apply gaussian filter
h = fspecial('gaussian',6,3);
im2_nobg_gau = imfilter(im2_nobg,h);

% Threshold and label Ch2
im2_t = areathresh(im2_nobg_gau, int_thresh_Ch2,area_thresh_Ch2,3,8) > 0;
im_labeled_Ch2 = bwlabel(im2_t);

% Number of punctas in Ch2
N_pun_Ch2 = max(im_labeled_Ch2(:));

% Checking if the channels are positive
Ch_pos_2 = nan(N_pun_Ch2, 2);
Ch_val_2 = nan(N_pun_Ch2, 2);

% Looping through Ch1 indices
for i = 1 : N_pun_Ch1
    Ch_pos_1(i,1) = mean(im_L2(im_labeled_Ch1 == i));
    Ch_pos_1(i,2) = mean(im2_t(im_labeled_Ch1 == i));
    
    Ch_val_1(i,1) = mean(im_P2(im_labeled_Ch1 == i));
    Ch_val_1(i,2) = mean(im2_nobg(im_labeled_Ch1 == i));
end

% Looping through Ch2 indices
for i = 1 : N_pun_Ch2
    Ch_pos_2(i,1) = mean(im_L2(im_labeled_Ch2 == i));
    Ch_pos_2(i,2) = mean(im2_t(im_labeled_Ch2 == i));
    
    Ch_val_2(i,1) = mean(im_P2(im_labeled_Ch2 == i));
    Ch_val_2(i,2) = mean(im2_nobg(im_labeled_Ch2 == i));
end

%% Outputs
% Ch2 puncta that are Ch1 pos
Ch1p_Ch2_ind = Ch_pos_2(:,1) > 0

% Ch2+ Ch1+ Ch1 values
Ch1p_Ch2p_val = Ch_val_2(Ch1p_Ch2_ind,1)

% Ch2+ Ch1 values
Ch2p_val = Ch_val_2(:,1)

% Ch2- Ch1 values
Ch2n_val = Ch_val_1(Ch_pos_1(:,2) == 0,1)

