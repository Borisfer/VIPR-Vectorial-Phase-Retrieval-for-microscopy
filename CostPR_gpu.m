function [out,grad,Nph_max_fact] = CostPR_gpu(phase_mask,IS,q,data,std_stack,Nph_vec,cost_function_flag...
    ,gBlur ,k,opt_phase_mat,g_bfp_init,vec_model_flag,int_cos,circ,circ_sc,Nph_opt_flag,cost_altr)

% function calculates the gradient at each point

%% inputs:

% phase_mask : the phase mask at the current iteration
% IS : structure with the optical parameters and hyperparameter(defined in main function)
% q : (x,y,z,NFP) position [um] per image (matrix size NX4)
% data :  full z_stack (cropped in a square odd grid)
% std_stack: estimated std per z-stack image
% Nph_vec: estimation of signal per z-stack image
% cost_function_flag: flag to choose cost function
% gBlur : gaussian blur 
% k : wavenumber
% opt_phase_mat : phase matrices whihc multiply q
% g_bfp_init : BFP amplitude components
% vec_model_flag :   1 - vectorial model, 0 - scalar model
% int_cos : apodization component
% circ : aperture matrix limited by NA
% circ_sc : aperture matrix limited by RI of medium


%% output
% out : cost at current point
% grad - gradient update to the phase mask 

%% preallocates
grad = phase_mask*0;
Iimg = phase_mask*0;

%% cropping of the fft to image size
N = size(circ,1);
Nph_max_fact = q(:,1)*0+1;

for z_ind = 1:size(q,1)
    % zero the image
    Iimg = Iimg*0;
    
    % position parameters
    x0 = q(z_ind,1);
    y0 = q(z_ind,2);
    z0 = q(z_ind,3);
    NFP = q(z_ind,4);
    % data at the z-stack
    data_z0 = (data(:,:,z_ind));
    % SNR 
    Nph = Nph_vec(z_ind);
    std_Gbg = std_stack(z_ind);
    
    % add noise to image 
%     Noise_mat = normrnd(0,std_Gbg,size(data_z0));
%     data_z0 = data_z0+Noise_mat;
    
    %% calc PSF
    % BFP phase pattern (position + phase mask)  
    BFP_Phase = k*(x0*opt_phase_mat(:,:,1)+y0*opt_phase_mat(:,:,2)+z0.*opt_phase_mat(:,:,3)-NFP*opt_phase_mat(:,:,4));
    BFP_Phase = exp(1i.*(BFP_Phase+phase_mask));
    
    if vec_model_flag
        g_bfp = g_bfp_init*0;
        g_img = g_bfp_init*0;
        
        if sum(abs(IS.p_vec) == 0) == 3 % freely rotating - superposition solution
            % BFP corrected intensity
            for g_id = 1:size(IS.g_bfp,3)
                g_bfp(:,:,g_id) = g_bfp_init(:,:,g_id) .*circ .*BFP_Phase;
            end
            
            for g_id = 1:size(IS.g_bfp,3)
                g_img(:,:,g_id) = 1/N.*(fft2(g_bfp(:,:,g_id)));
                Iimg = Iimg+g_img(:,:,g_id).*conj(g_img(:,:,g_id));
            end
        else % coherent dipole 
            % normalization
            for div_pol = 1:floor(size(IS.g_bfp,3)/3)
                for g_id = (1:3)+3*(div_pol-1)
                    g_bfp(:,:,g_id) = g_bfp_init(:,:,g_id) .*BFP_Phase.*circ;
                end
            end
            for div_pol = 1:floor(size(IS.g_bfp,3)/3)
                for g_id = (1:3)+3*(div_pol-1)
                   
                    g_img(:,:,g_id) = 1/N.*(fft2((g_bfp(:,:,g_id))).*IS.p_vec(g_id-3*(div_pol-1)));
                end
                Iimg = Iimg+abs(sum(g_img,3)).^2;
            end
        end
        
    else % scalar model 
        
        % scalar bfp
        Ebfp = BFP_Phase.*circ_sc.*int_cos;
        
        % generate image from pupil function
        Eimg = 1/N*fft2(Ebfp); % image plane
        Iimg = Eimg.*conj(Eimg); % intensity at image plane
    end
    
    %% shift and blur the PSF
    tmp  = (fftshift(real(Iimg)));
    %% crop to correct intensity 
    % crop image
    if mod(size(tmp,1),2)==1
        if mod((IS.FOV_size_crop)/2,1)
            crop_tmp = real(tmp(end/2-floor(IS.FOV_size_crop)/2+1:end/2+floor(IS.FOV_size_crop)/2,end/2-floor(IS.FOV_size_crop)/2+1:end/2+floor(IS.FOV_size_crop)/2));
        else
            crop_tmp = real(tmp(end/2-floor(IS.FOV_size_crop)/2+0.5+1:end/2+floor(IS.FOV_size_crop)/2+1,end/2-floor(IS.FOV_size_crop)/2+0.5+1:end/2+floor(IS.FOV_size_crop)/2+1));
        end
    else
        crop_tmp = real(tmp(round(size(tmp,1)/2-round(IS.FOV_size_crop-1)/2+1/2):round(size(tmp,1)/2+floor(IS.FOV_size_crop-1)/2+1/2),round(size(tmp,2)/2-round(IS.FOV_size_crop-1)/2+1/2):round(size(tmp,2)/2+floor(IS.FOV_size_crop-1)/2+1/2)));
    end
    %
    normfact = Nph./sum(crop_tmp(:));
    tmp = tmp.*normfact;
    % blur
    tmp = imfilter(tmp,gBlur,'replicate');
    %% correct intensity after crop
    if Nph_opt_flag
        %      
        sort_dat = sort(data_z0(:));
        sort_tmp = sort(tmp(:));
        Nph_max_fact(z_ind) = mean((sort_dat(end-50:end))./sort_tmp(end-50:end));
                
        % limit to a reasonable 0.85 to 1.2
        if Nph_max_fact(z_ind) < 0.85
            Nph_max_fact(z_ind) = 0.85;
        elseif Nph_max_fact(z_ind) > 1.2
            Nph_max_fact(z_ind) = 1.2;
        end
        
        %
        tmp = tmp.*Nph_max_fact(z_ind);
    end
    %% calc  cost and derivative  
    if cost_function_flag==1
        %% L1 loss
        cost(z_ind) = sum(abs(tmp(:) - data_z0(:)));
        % derivative
        dcost_dtmp = sign(tmp - data_z0).*(data_z0~=0);
    elseif cost_function_flag==2
        %% L2 Loss
        cost(z_ind) = sum((tmp(:) - data_z0(:)).^2);
        % derivative
        dcost_dtmp = 2*(tmp - data_z0).*(data_z0~=0);
    elseif cost_function_flag==3
        %% MLE
        cost(z_ind) = -sum(data_z0(:).*log(tmp(:))-tmp(:));
        % derivative
        dcost_dtmp = (1-data_z0./tmp).*(data_z0~=0);
    elseif cost_function_flag==4
        % double gaussian MLE
        cost(z_ind) = sum((log(sqrt(2*pi*(tmp(:)+std_Gbg.^2+eps)))+0.5*((data_z0(:)-tmp(:)).^2)./(tmp(:)+std_Gbg.^2)).*(data_z0(:)~=0));
        % derivative
        denom = 2*(tmp+std_Gbg.^2+eps).^2;
%         nom = ((tmp+std_Gbg.^2)-2*(tmp+std_Gbg.^2).*(data_z0-tmp)-(data_z0-tmp).^2).*(data_z0~=0);
        nom = ((tmp+std_Gbg.^2)-2*(tmp+std_Gbg.^2).*(data_z0-tmp)-(data_z0-tmp).^2);
        dcost_dtmp = nom./denom;

    end
    cost(z_ind) = cost(z_ind)*cost_altr(z_ind);
    % filter the derivative
    dcost_dIimg = cost_altr(z_ind)*imfilter(dcost_dtmp,gBlur,'replicate');
  
    %% calc gradient
    if vec_model_flag
        if sum(abs(IS.p_vec) == 0) == 3 % freely rotating - superposition solution
            for g_id = 1:size(g_img,3)
                grad_tmp(:,:,g_id) = 2*1/N*real((fft2(ifftshift(dcost_dIimg).*1i.*conj(g_img(:,:,g_id))).*g_bfp(:,:,g_id))).*Nph_max_fact(z_ind).*normfact; %%
            end
        else % coherent 
            for div_pol = 1:floor(size(IS.g_bfp,3)/3)
                Conj_fact = conj(sum(g_img(:,:,(1:3)+3*(div_pol-1)),3));
                for g_id = (1:3)+3*(div_pol-1)
                    grad_tmp(:,:,g_id) = 2*1/N*real(fft2(ifftshift(dcost_dIimg).*1i.*Conj_fact).*g_bfp(:,:,g_id).*IS.p_vec(g_id-3*(div_pol-1))).*Nph_max_fact(z_ind).*normfact; %%
                end
            end
        end
        grad(:,:,z_ind) = sum((grad_tmp),3);
    else % scalar model 
        grad(:,:,z_ind) = 2*1/N*real(fft2(ifftshift(dcost_dIimg).*1i.*conj(Eimg)).*Ebfp).*Nph_max_fact(z_ind).*normfact;
    end
    
end

out = mean(cost);
% remove SAF artifacts if the model is scalar
if IS.SAF_flag 
    grad = sum(grad,3);
else
    grad = sum(grad,3).*circ_sc;
end





















