%% load data TP 4um PR SLM aligned 
Cpath = [pwd,'\TP images\'];

filelist_n = dir(Cpath);
filelist_n = filelist_n(~cellfun(@(x) x==0, {filelist_n.bytes}));
filelist_n = filelist_n(cellfun(@(x) x(1)=='T', {filelist_n.name}));

for qq = 1:size(filelist_n,1)
    stack = bfopen([Cpath,filelist_n(qq).name]);
    metadata = imfinfo([Cpath,filelist_n(qq).name]);
    stack = stack{:,1};
    %(camera alignment position compared to the model)
    IMG(:,:,qq) = flipud(rot90(double(stack{1,1})));
    % read z position from file (stage read)
    z_stack_pos(qq) = metadata.UnknownTags(5).Value(3);
end
% correct z-positions offset 
z_stack_pos = (z_stack_pos-z_stack_pos(1)+ 3);


%% select ROI
IS.FOV_size = 80;
%pre - defined coordinates 
R = [37 36 80 80];
W = 80;
% cut FOV 
IMG_T = IMG(round(R(2)):round(R(2))+W,round(R(1)):round(R(1))+W,:);



