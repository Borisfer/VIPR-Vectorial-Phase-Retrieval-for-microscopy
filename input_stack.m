function [FOV_size,IMG_T] = input_stack(crop_flag,FOV_size)

% choose tiff_stack
[filename, PATHNAME] = uigetfile([pwd,'../*.tif'],'title','Choose tiff stack');
if filename == 0
    msgbox('Bad path name' );
    return
end

metadata = imfinfo([PATHNAME,filename]);
f = waitbar(0,'Wait for data to load');
for j = 1:size(metadata,1)
    waitbar(j*100/size(metadata,1),f,'Wait for data to load');
    IMG_T(:,:,j) = double(imread([PATHNAME,filename],j));
end
close(f)

if size(IMG_T,1) ~= size(IMG_T,2)
    msgbox('Chosen FOV is not square');
    return
end
%% crop the data
%% manual selection
if crop_flag == 1
    %% select ROI
    h = figure;
    imagesc(IMG_T(:,:,round(end/2)));daspect([1,1,1]);
    % h_rect = imrect(gca, [1 1 FOV_size FOV_size]); % create rectangle on the image
    % setResizable(h_rect,false);
    message = sprintf('Drag and click on the center pixel choose ROI');
    % uiwait(msgbox(message));
    % % R = [     10     6   114   114];
    [x,y] =  getpts(h);
%     R = [ceil(floor(x)-FOV_size/2+0.5)     ceil(floor(y)-FOV_size/2+0.5)   FOV_size   FOV_size];
    % R = round(wait(h_rect)); % get position 
    W = min([FOV_size,min(round(floor(x)+FOV_size/2),size(IMG_T,2))-max(round(floor(x)-FOV_size/2),1)...
        ,min(round(floor(y)+FOV_size/2),size(IMG_T,1))-max(round(floor(y)-FOV_size/2),1)]);
    %
    CL_x = max(round(floor(x)-W/2),1);
    CL_y = max(round(floor(y)-W/2),1);
    CR_x = min(round(floor(x)+W/2),size(IMG_T,2));
    CR_y = min(round(floor(y)+W/2),size(IMG_T,1));
    %
    if CR_x-CL_x ~= CR_y-CL_y
            msgbox('Chosen FOV size is too large, image not cropped');
            return
    end
    %
    IMG_T = IMG_T(CL_y:CR_y,CL_x:CR_x,:);
    %
    FOV_size = size(IMG_T,1);
    %
    close(h)
elseif crop_flag ==2
    % create max-projection
    IMG_maxP = max(IMG_T,[],3);
    % calculate CoG
    [X,Y] = meshgrid([1:size(IMG_T,1)]);
    %
    [x] = round(sum(IMG_maxP(:).*X(:))./sum(IMG_maxP(:)));
    [y] = round(sum(IMG_maxP(:).*Y(:))./sum(IMG_maxP(:)));
    %
    R = [floor(x)-FOV_size/2     floor(y)-FOV_size/2   FOV_size   FOV_size];
    % R = round(wait(h_rect)); % get position
    W = max(R(3),R(4));
    %
    %
    IMG_T = IMG_T(max(1,round(R(2))):min(round(R(2))+W,size(IMG_T,1)),...
        max(1,round(R(1))):min(round(R(2))+W,size(IMG_T,2)),:);
    %
    
    %
    FOV_size = size(IMG_T,1);
else
    FOV_size = size(IMG_T,1);
end

%% select ROI

