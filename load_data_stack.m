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
h = figure;
imagesc(IMG(:,:,round(end/2)));daspect([1,1,1]);
FOV_size = IS.FOV_size;
h_rect = imrect(gca, [1 1 FOV_size FOV_size]); % create rectangle on the image
setResizable(h_rect,false);
message = sprintf('Drag and double click on the rectangle box to choose ROI');
uiwait(msgbox(message));
R = round(wait(h_rect)); % get position
W = max(R(3),R(4));

IMG_T = IMG(round(R(2)):round(R(2))+W,round(R(1)):round(R(1))+W,:);
FOV_size = size(IMG_T,1);

close(h)


