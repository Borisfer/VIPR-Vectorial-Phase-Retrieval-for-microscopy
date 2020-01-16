
%% load your data here (automatically when Main is ran)
[FOV_size,IMG_T] = input_stack;
IS.FOV_size = FOV_size; % size of ROI used

%% load image coordinates here (change to match your data)
% NFP positions input
load(['TP images','\NFP_pos_TP_images.mat']);
z_stack_pos = z_stack_pos';
% z position of molecule (default to emitter radius)
z_pos = zeros(length(z_stack_pos),1)+IS.z_emit;
% xy default 0
xy = zeros(length(z_stack_pos),2);

%
if size(xy,1) ~= size(IMG_T,3) || size(z_pos,1) ~= size(IMG_T,3) || size(z_stack_pos,1) ~= size(IMG_T,3)
    msgbox('coordinate matrix does not match in size with amount of loaded images');
    return
end