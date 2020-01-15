%% load data TP 4um PR SLM aligned 
Cpath = [pwd,'\TP images\'];

filelist_n = dir(Cpath);
filelist_n = filelist_n(~cellfun(@(x) x==0, {filelist_n.bytes}));
filelist_n = filelist_n(cellfun(@(x) x(1)=='T', {filelist_n.name}));

metadata = imfinfo([Cpath,filelist_n.name]);
for qq = 1:size(metadata,1)
    stack = imread([Cpath,filelist_n.name],qq);
    IMG(:,:,qq) = (double(stack));
    % read z position from file (stage read)
end
% correct z-positions offset 
load([Cpath,'NFP_pos_TP_images.mat']);
z_stack_pos = z_stack_pos';

%% select ROI
IS.FOV_size = 80;
%pre - defined coordinates 
R = [1 1 80 80];
W = 79;
% cut FOV 
IMG_T = IMG(round(R(2)):round(R(2))+W,round(R(1)):round(R(1))+W,:);



