% clear all;close all
N_crop = 39;

%% load data images
Cpath = [pwd,'\EPFL DH data\'];
fname = 'sequence-as-stack-Beads-DH-Exp.tif';
info = imfinfo([Cpath,fname]);
imageStack = [];
numberOfImages = length(info);
for k = 1:numberOfImages
    currentImage = imread([Cpath,fname], k, 'Info', info);
    imageStack(:,:,k) = currentImage;
    %     imagesc(imageStack(:,:,k));
    %     drawnow
end

%% load xls data
[num,txt,raw] = xlsread([Cpath,'activations.csv']);
XY_data = num(:,3:4)/1000;
Z_data = num(:,5)/1000;
%% create grid
psize = 0.1;
line_x = (0:size(imageStack,1))*psize;
[X,Y] = meshgrid(line_x);

%% cut stack
xy=[];
for i = 1:6
    for j = 1:numberOfImages
        x_pos_nm = XY_data((i-1)*151+j,1);
        y_pos_nm = XY_data((i-1)*151+j,2);
        
        x_pos_px = round(x_pos_nm./psize);
        y_pos_px = round(y_pos_nm./psize);
        
        %residuals
        xy(j,1,i) = x_pos_nm-round(x_pos_nm./psize)*psize;
        xy(j,2,i) = y_pos_nm-round(y_pos_nm./psize)*psize;
        z(j,i) = Z_data((i-1)*151+j);
        
        
        tmp = imageStack(y_pos_px-floor(N_crop/2):y_pos_px+floor(N_crop/2),x_pos_px-floor(N_crop/2):x_pos_px+floor(N_crop/2),j);
        % rescale to match camera size
        %         DH_PSF(:,:,j,i) = imresize(tmp,10/16,'bilinear');
        tmp = tmp-min(tmp(:));
        DH_PSF(:,:,j,i) = tmp.*(tmp>0.03*max(tmp(:)));
        %         imagesc( DH_PSF(:,:,j,i));
        %         drawnow;
    end
end











