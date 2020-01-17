function [FOV_size,IMG_T] = input_stack

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

%% select ROI
FOV_size = size(IMG_T,1);

