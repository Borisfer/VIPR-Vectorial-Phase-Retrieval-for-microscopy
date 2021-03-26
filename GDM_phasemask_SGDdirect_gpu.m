function [x,gB_est,Nph] = GDM_phasemask_SGDdirect_gpu(data,fun, x,GenFunc,IS,opt_flag,plot_flag,est_gBlur_flag,q,std_stack0,Cost_fun_init,step_size,gpu_flag,noisy_flag)

% function calculates the estimated phase mask
%% inputs:
% data :  full z_stack (cropped in a square odd grid)
% fun: function handle to calculate the gradient
% x : initial Phase mask
% GenFunc : function handle to plot PSFs
% IS : structure with the optical parameters and hyperparameter(defined in main function)
% opt_flag : gradient scheme % 1 - ADAM, 2 - Adamax , 3- Nadam, 4 - Nesterov, 5- SVRG, 6- Vanilla SGD
% plot_flag: 1-plot every 30 iterations, 0- don't plot anything (faster)
% est_gBlur_flag : 1- estimate gBlur after 1/3 of the iterations
% q : (x,y,z,NFP) position [um] per image (matrix size NX4)
% std_stack0: estimated std per z-stack image
% Cost_fun_init: function handle to calculate initial point
% step_size: step size for optimization
% noisy_flag : 0 for design PSFs, 1 for PR;
% gpu_flag : 1 - use GPU, 0 - use CPU

%% output
% x - estimated phase mask
% gB_est - estimated blur kernel
% Nph - estimated signal per image

setappdata(0,'UseNativeSystemDialogs',false);
iter = 1;
f1=zeros(size(x));f2=f1;V=f1;f11=[];f22=[];grad=zeros(size(x,1));O=[];
load('periodiccolormap2');

%% opt params
maxiter = IS.SGDiter+IS.last_iter;
beta1 = 0.9; % adam g
beta2 = 0.999; % adam g^2
count_noise = 0;
epsilon = 1e-8;
gamma = 0.9; % momentum decay
mu_hat = f1;
Count_sample = 20;
Nph_opt_flag = 0;
count_Nph = 0;
gB = IS.gBlur;
std_stack = std_stack0;
%
point_num = min(IS.point_num,size(data,3));
%% estimate intensities 
% image sums without mean bg corrected to matching
for z_ind = 1:size(q,1)
    tmp = data(:,:,z_ind);
    Nph(z_ind) = sum(sum((tmp).*(tmp~=0)));
end
% 
if IS.Photobleach_cost_mult ==1
    p0 = polyfit(1:length(Nph),Nph./max(Nph),1);
    cost_altr = fliplr(polyval(p0,1:length(Nph)));
else
    cost_altr = ones(1,length(Nph));
end
%% run non repeating SGD per Epoch
ind_opt_perm = randperm(size(q,1));
Epoch_length = floor(size(q,1)/point_num);
gB_est = IS.gBlur;

gBlur_gpu = ( fspecial('gaussian',[5 5],gB_est));
if gpu_flag==1
    gBlur_gpu = gpuArray(gBlur_gpu);
end

sample_flag = 1;
%% do init if no init mask was given
if sum(x(:)) == 0
    [guess1] = Cost_fun_init(q,data,std_stack0,Nph,gBlur_gpu);
    x = x - step_size*guess1/10^6;
end
%% define normalization inline function for plotting 
imn = @(im) (im - min(im(:)))./(max(im(:)) - min(im(:)));
%% run iterations
while iter<=maxiter
    %   randominze noise and blur
    % randomize the noise in the slice
    if noisy_flag
        std_stack = max(std_stack0 + normrnd(0,mean(std_stack0)/8,1),mean(std_stack0)/2);
        % randomize Blur kernel
        gB = min(max(normrnd(gB_est,gB_est/8),1/2*gB_est),1.5);
        gBlur_gpu = ( fspecial('gaussian',[5 5],gB));
        if gpu_flag==1
            gBlur_gpu = gpuArray(gBlur_gpu);
        end
    end
    
    %% Pick z positions
    pool_z = [];
    
    if sample_flag==1 % SGD
        ind_opt = ind_opt_perm(1:point_num);
        ind_opt_perm(1:point_num) = [];
    elseif sample_flag ==2 % global gradient
        ind_opt = 1:size(q,1);
    elseif sample_flag ==3 % at worst corr2
        % directed to correct at bad corr2 positions
        grid_sampling_space = 1; %3
        numz = size(q,1);
        grid_sampling = 1:grid_sampling_space:numz;
        corr_coarse = (zeros(1,length(grid_sampling)));
        
        for z_ind = 1:length(grid_sampling)
            %             zz = q(grid_sampling(z_ind));
            tmp = Nph(grid_sampling(z_ind))*real(GenFunc(x,q(grid_sampling(z_ind),:)));
            tmp = imfilter(tmp,gBlur_gpu,'replicate');
            %crop
            tmp_data = gather(data(round(end/2-IS.FOV_size/2)+1:round(end/2+IS.FOV_size/2),round(end/2-IS.FOV_size/2)+1:round(end/2+IS.FOV_size/2),z_ind));
            % thr mask
            
            mask = tmp_data>max(tmp_data(:)).* IS.thr_corr;
            tmp = gather(tmp.*mask);
            tmp_data = gather(tmp_data.*mask);  
            % corr2
            corr_coarse(z_ind) = corr2((tmp)./norm(tmp(:)),(tmp_data)./norm(tmp_data(:)));
        end
        [~,ind_opt] = sort(corr_coarse);
        ind_opt = [ind_opt(1:point_num)];
    elseif sample_flag ==4 % adaptive sampling
        ind_opt=[];
        Q_idx=[];
        Q_opt=[];
        % seperate into quarters
        Q_idx = round(linspace(1,size(data,3),7));
        % directed to correct at bad corr2 positions
        grid_sampling_space = 1; %3
        numz = size(q,1);
        grid_sampling = 1:grid_sampling_space:numz;
        
        if Count_sample==20
            corr_coarse = (zeros(1,length(grid_sampling)));
            
            for z_ind = 1:length(grid_sampling)
                %             zz = q(grid_sampling(z_ind));
                tmp = Nph(grid_sampling(z_ind))*real(GenFunc(x,q(grid_sampling(z_ind),:)));
                tmp = imfilter(tmp,gBlur_gpu,'replicate');
                % crop
                tmp_data = (data(round(end/2-IS.FOV_size/2)+1:round(end/2+IS.FOV_size/2),round(end/2-IS.FOV_size/2)+1:round(end/2+IS.FOV_size/2),z_ind));                % thr mask
                % thr
                mask = tmp_data>max(tmp_data(:)).* IS.thr_corr;
                tmp = gather(tmp.*mask);
                tmp_data = gather(tmp_data.*mask);
                % corr2
                corr_coarse(z_ind) = (1-corr2((tmp)./norm(tmp(:)),(tmp_data)./norm(tmp_data(:))));
            end
        else
            Count_sample=0;
        end
        %
        for j=1:length(Q_idx)-1
            norm_accu_Q(j) = (mean(corr_coarse(Q_idx(j):Q_idx(j+1))));
        end
        % normalization
        p_accu_grad = norm_accu_Q./sum(norm_accu_Q);
        CDF_sampling = [0,cumsum(p_accu_grad)];
        % generate uniform random points
        R_uni = rand(point_num,1);
        CDF_interps = interp1(CDF_sampling,[0:1/(length(CDF_sampling)-1):1],R_uni'); %spans zero to one
        %new uniform dest;
        Q_opt = ceil(CDF_interps*(length(CDF_sampling)-1));
        %sample uniformly from the selected Qs
        for j=1:point_num
            Q = Q_opt(j);
            ind_opt(j) = randi(Q_idx(Q+1)-Q_idx(Q)+1)+Q_idx(Q)-1;
        end
        ind_opt = unique(ind_opt);
        Count_sample = Count_sample+1;
    end
    %     ind_opt = 1:point_num;
    data_pool = data(:,:,(ind_opt));
    q_pool = q((ind_opt),:);
    std_pool = std_stack((ind_opt));
    Nph_pool = Nph(ind_opt);
    cost_altr_pool = cost_altr(ind_opt);
    
    %% comp grad
    % grad
    [out,grad,Nph_max_fact] = fun(x,q_pool,data_pool,std_pool,Nph_pool,gBlur_gpu,Nph_opt_flag,cost_altr_pool);
    
    % correct intensity estimate 
    Nph(ind_opt) = Nph(ind_opt).*Nph_max_fact';
    % update cost
    O(end+1) = gather(out);%
    
    % use gradient scheme
    if opt_flag ==1 %ADAM
        f1 = beta1.*f1 + (1-beta1).*grad;
        f1_ub=f1./(1-beta1.^(iter));
        f2 = beta2.*f2 + (1-beta2).*(grad.^2);
        f2_ub=f2./(1-beta2.^(iter));
        %    calc the final V step
        V = step_size*(f1_ub)./(sqrt(f2_ub)+epsilon);
    elseif opt_flag == 2 %Adamax
        f1 = beta1.*f1 + (1-beta1).*grad;
        f2 = bsxfun(@max,beta2.*f2,abs(grad)).*IS.circmask_large;
        V = step_size*(f1)./(f2+eps).*IS.circmask_large;
    elseif opt_flag ==3 % Nadam
        f1 = beta1.*f1 + (1-beta1).*grad;
        %         f1_hat = f1./(1-beta1.^iter);
        f2 = beta2.*f2 + (1-beta2).*grad.^2;
        %         f2_hat  = f2./(1-beta2.^iter);
        V = step_size./(sqrt(f2)+epsilon).*(beta1.*f1+(1-beta1).*grad).*IS.circmask_large;
    elseif opt_flag == 4 % Nesterov
        V = gamma.*V+step_size.*grad;
    elseif opt_flag == 5 % vanilla SGD
        V = step_size.*(grad);
    end
        
    % update mask
    x_new = (x - V);
    %update
    x = x_new;
    %     masks(:,:,iter) = x;
    %     fprintf(['cost = ',num2str(out),'\n']);
    count_noise = count_noise+1;
    
    %% plot every 30 iterations
    % phase mask
    if  mod(iter,30) == 1 && plot_flag
        
        fig= findobj('name','Cost');
        if ~isempty(fig)
            figure(fig);
        else
            f=figure('units','inches');
            pos = get(gcf,'pos');
            set(gcf,'pos',[pos(1)-4 pos(2)-3 4 2]);
            set(f,'NumberTitle','off','name','Cost');
        end
        
        plot(O)
        xlabel('iter');
        ylabel('cost');
        
        % sample PSF's
        countS=1;
        fig= findobj('name','Sample PSFs');
        if ~isempty(fig)
            figure(fig);
        else
            f=figure('units','inches');
            pos = get(gcf,'pos');
            set(gcf,'pos',[pos(1)+4 pos(2)-4 8 6]);
            set(f,'NumberTitle','off','name','Sample PSFs');
        end
        
        for z_ind = 1:length(ind_opt)
            zz = (q_pool(z_ind,4));
            % calc the PSF
            tmp = Nph_pool(z_ind)*GenFunc(x,q_pool(z_ind,:));
            tmp = imfilter(tmp,gBlur_gpu,'replicate');
            tmp = tmp(round(end/2-IS.plotsize/2)+1:round(end/2+IS.plotsize/2),round(end/2-IS.plotsize/2)+1:round(end/2+IS.plotsize/2));
            % plots
            subplot(length(ind_opt),2,countS);
            countS=countS+1;
            imagesc([imn(real(tmp))]);
            colormap(hot)
            title(['model z = ',num2str(zz)]);
            daspect([1 1 1]);colorbar;
            subplot(length(ind_opt),2,countS);
            countS=countS+1;
            imagesc([imn(data(round(end/2-IS.plotsize/2)+1:round(end/2+IS.plotsize/2),round(end/2-IS.plotsize/2)+1:round(end/2+IS.plotsize/2),(ind_opt(z_ind))))]);
            title(['data z = ',num2str(zz)]);
            daspect([1 1 1]);colorbar;
            colormap(hot)
        end
        % plot the current phase mask
        fig= findobj('name','Phase mask');
        if ~isempty(fig)
            figure(fig);
        else
            f=figure('units','inches');
            pos = get(gcf,'pos');
            set(gcf,'pos',[pos(1)-4 pos(2) 4 4]);
            set(f,'NumberTitle','off','name','Phase mask');
        end
        
        % plot phase mask
        x_plot = gather(real(x));
        x_plot = mod(x_plot,2*pi)-pi;
        imagesc(x_plot,'alphadata',IS.circmask_large);
        axis image;
        colormap(periodiccolormap);
        colorbar;
        title('current phase mask');
        drawnow;
    end

    %% start correcting image sum
    if (iter == floor(IS.SGDiter/4)) && IS.update_Signal==1
            Nph_opt_flag = 1;
    end
    %% estimate gblur
    if (est_gBlur_flag>0 && (iter == floor(IS.SGDiter/3)))
        % generate PSF stack
        gBlur_gpu = ( fspecial('gaussian',[5 5],gB_est));
        if gpu_flag==1
            gBlur_gpu = gpuArray(gBlur_gpu);
        end
        
        %% create the current full data set
        for j = 1:size(data,3)
            % create model and data(with no bg) stacks on CPU
            model_stack(:,:,j) = gather(Nph(j)*(GenFunc(x,q(j,:))));
            % crop data
            data_stack(:,:,j) = gather(data(round(end/2-IS.FOV_size/2)+1:round(end/2+IS.FOV_size/2),round(end/2-IS.FOV_size/2)+1:round(end/2+IS.FOV_size/2),j));
            %
            thr_dat(j) = gather(max(max(data_stack(:,:,j))).*0.1);% create thr vec
        end
        
        %% estimate gBlur
        options = optimoptions(@fmincon,'Display','iter','gradObj','off','algorithm','sqp');
        
        cost = @(gB) CostgBlur(gB,data_stack,model_stack,IS.gBlur_cost,gather(std_stack0),(thr_dat));
        gB_est = fmincon(cost,IS.gBlur,[],[],[],[],0.3,10,[],options);
        % create new blur kernel
        gBlur_gpu = ( fspecial('gaussian',[5 5],gB_est));
        if gpu_flag==1
            gBlur_gpu = gpuArray(gBlur_gpu);
        end
        %
        clear model_stack
        clear data_stack
    end
    
    %% re-permute each epoch
    if count_noise == Epoch_length
        ind_opt_perm = randperm(size(q,1));
        count_noise = 0;
        %dont overfit the  noise
        if Nph_opt_flag==1
            count_Nph = count_Nph+1;
        end
        if count_Nph == 2
            Nph_opt_flag = 0;
        end
    end
    
    %% switch sampling proc
    if sample_flag == 1 && iter == IS.SGDiter
%         noisy_flag = 0;
        sample_flag = IS.last_iter_flag;
        %create constant blur and noise for the rest of iterations 
        gB = gB_est;
        std_stack = std_stack0;
        gBlur_gpu = ( fspecial('gaussian',[5 5],gB));
        if gpu_flag==1
            gBlur_gpu = gpuArray(gBlur_gpu);
        end
    end
    iter = iter+1;
end
