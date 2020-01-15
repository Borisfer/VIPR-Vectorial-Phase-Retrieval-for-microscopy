
%% load your data here
[FOV_size,IMG_T] = input_stack;
IS.FOV_size = FOV_size; % size of ROI used

%% load image coordinates here
% NFP positions input
load(['TP images','\NFP_pos_TP_images.mat']);
z_stack_pos = z_stack_pos';
%     z is always emitter radius
z_pos = zeros(length(z_stack_pos),1)+IS.z_emit;
% xy are 0
xy = zeros(length(z_stack_pos),2);

%
if size(xy,1) ~= size(IMG_T,3) || size(z_pos,1) ~= size(IMG_T,3) || size(z_stack_pos,1) ~= size(IMG_T,3)
    msgbox('Chosen FOV is not square');
    return
end