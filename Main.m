%%
% written by Boris Ferdman
% this functions performs phase retrieval via a direct numerical gradient
% associated paper: VIPR: Vectorial Implementation of Phase Retrieval for
% fast and accurate microscopic pixel-wise pupil estimation 2020

% License to use this code is granted to all interested, as long as the original author is
% referenced as such. The original author maintains the right to be solely associated with this work.
% Copyright by Boris Ferdman: borisferd@gmail.com
%
clear all;close all;clc
%% data set
data_set = 2; %1 - your data, 3-EPFL DH 1.5[um], 2 - 4[um] Tetrapod mask 
%% all user input defined here
VIPR_user_input;

%% load data
if data_set == 1
    VIPR_load_data;
elseif data_set == 2
    % load data and z positions
    load_data_stack;
    % z is always emitter radius
    z_pos = zeros(length(z_stack_pos),1)+IS.z_emit;
    % xy are 0
    xy = zeros(length(z_stack_pos),2);
elseif data_set == 3
    %load data
    EPFL_data_load;
    % how much to sample
    dI = 1;
    % choose a single stack
    IMG_T = DH_PSF(:,:,1:dI:end,5);
    % take minus z to match to NFP
    z_stack_pos = -z(1:dI:end,1);
    % image centering
    xy = xy(1:dI:end,:,1);
    z_pos = zeros(length(z_stack_pos),1)+IS.z_emit;
    IS.FOV_size = size(IMG_T,1);
    % do scalar model because we don't know really what the system is
    vec_model_flag = 0;
end

% define positions per image - [x,y,z,NFP]
q_cord = [xy';z_pos';z_stack_pos']';

%% Gaussian noise estimation per Z
if noisy_flag
    dx=IS.corner_size;
    for j = 1:size(IMG_T,3)
        tmp = [IMG_T(1:dx,1:dx,j);IMG_T(end-dx+1:end-1+1,end-dx+1:end-1+1,j);IMG_T(end-dx+1:end-1+1,1:dx,j);IMG_T(1:dx,end-dx+1:end-1+1,j)];
        %
        mean_stack(j) = mean(tmp(tmp>0));
        % Estimate the std
        std_stack(j) = std(tmp(tmp>0));
        
    end
    IMG_bu = IMG_T;
else
    IMG_bu = IMG_T;
    std_stack = zeros(1,size(q_cord,1));
    mean_stack = zeros(1,size(q_cord,1));
end
%% thr and reduce offset from the images
IMG_T = IMG_bu;
%
for j = 1:size(IMG_T,3)
    tmp = IMG_T(:,:,j)-mean_stack(j);
    % threshold on %
    if IS.I_thr_flag==1
        maskT = tmp>IS.I_thr.*max(tmp(:));
    else
        maskT = tmp>IS.I_thr.*std_stack(j);
    end
    % erode&dilate mask
    se = strel('disk',1,6);
    erodedI = imerode(double(maskT),se);
    se = strel('disk',1,6);
    mask(:,:,j) = imdilate(erodedI,se);
    
    IMG_T(:,:,j) = tmp.*mask(:,:,j);
    
    ff=figure(11)
    subplot(1,3,1)
    imagesc(tmp);
    axis  image
    title('input image');
    subplot(1,3,2)
    imagesc(IMG_T(:,:,j));
    axis  image
    title('thresholded image');
    subplot(1,3,3)
    imagesc(IMG_T(:,:,j)-tmp);
    axis  image
    title('diff');
end
close(11)

%% upsample if needed (usually if object space pixels are large compared to the wavelength)
tmp=[];
if IS.upsample_fact>1
    IS.upsample_fact = round(IS.upsample_fact); % only integers
    % change pixel size
    IS.Cam_psize = IS.Cam_psize./IS.upsample_fact;
    %upsample the data
    [Xg,Yg] = meshgrid([1:size(IMG_T,1)]-size(IMG_T,1)/2-0.5);
    %upsampled grid
    [Xup,Yup] = meshgrid([1:1/IS.upsample_fact:size(IMG_T,1)]-size(IMG_T,1)/2-0.5);
    %upsample it
    for j=1:size(IMG_T,3)
        tmp(:,:,j) = interp2(Xg,Yg,IMG_T(:,:,j),Xup,Yup,'spline');
    end
    %     tmp=imresize(IMG_T,IS.upsample_fact,'box');
    % make an odd grid
    if mod(size(tmp,1),2)==0
        IMG_T=tmp(1:end-1,1:end-1,:);
    else
        IMG_T = tmp;
    end
    %increase blur
    IS.gBlur = IS.gBlur*IS.upsample_fact;
    
end

%% Fresnel matching of the grid
% create an odd grid  size
if mod(floor(IS.lambda(1)*IS.f_4f./(IS.Cam_psize.*IS.SLM_psize)),2)==1
    IS.maskDiam_m = IS.lambda*IS.f_4f./(IS.Cam_psize*IS.SLM_psize)*IS.SLM_psize;
else
    IS.maskDiam_m = IS.lambda*IS.f_4f./(IS.Cam_psize*IS.SLM_psize)*IS.SLM_psize+IS.SLM_psize;
end
%% add zeros if BFP is too small
if IS.maskDiam_m./IS.SLM_psize < IS.FOV_size
    % fix this with box interpolation
    IS.maskDiam_m = IS.FOV*IS.SLM_psize;
end

%
mask_size = floor(IS.maskDiam_m/IS.SLM_psize);
%% pad data to match simulation size
tmp_stack = [];
for j = 1:size(IMG_T,3)
    tmp = padarray(IMG_T(:,:,j), floor([mask_size - (size(IMG_T,1)), mask_size - (size(IMG_T,2))]/2),'both');
    if mod(mask_size,2)==0
        tmp = padarray(tmp, floor([1,1]),'pre');
    end
    if mod(size(tmp,1),2)==0
        tmp = padarray(tmp, floor([1,1]),'pre');
    end
    % complete 0 to mean noise
    tmp_stack(:,:,j) = tmp;
    tmp = [];
end
% IMG_T = tmp_stack;

%% PR the mask
[maskRec,gB,Nph,I_mod] = PR_coverslip_data(IS,tmp_stack,q_cord,std_stack,gpu_flag,vec_model_flag,cost_function_flag,plot_flag,Alg_flag,est_gBlur_flag,noisy_flag,vec_model_pol);
%% outputs
% maskRec - the PR phase mask
% gB - estimation of the blur (if was enabled)
% Nph - estimation of image sum
% I_mod - simulated stack
