function [grad] = CostPR_gpu_init(IS,q,data,std_stack,Nph_init...
     ,k,opt_phase_mat,g_bfp_init,vec_model_flag,int_cos,circ,circ_sc)

% function calculates the initial point 
%% inputs:
% q : (x,y,z,NFP) position [um] per image (matrix size NX4)
% IS : structure with the optical parameters and hyperparameter(defined in main function)
% data :  full z_stack (cropped in a square odd grid)
% std_stack: estimated std per z-stack image
% Nph_init: estimation of signal per z-stack image
% opt_phase_mat : phase matrices whihc multiply q
% g_bfp_init : BFP amplitude components
% vec_model_flag :   1 - vectorial model, 0 - scalar model
% int_cos : apodization component
% circ : aperture matrix limited by NA
% circ_sc : aperture matrix limited by RI of medium


%% output
% grad : initialization of phase 

%% preallocates
grad = circ*0;
Iimg = circ*0;

% cropping of the fft to image size
N_crop = min(IS.FOV_size,size(circ,1));
N = size(circ,1);


for z_ind = 1:size(q,1)
    % zero the image
    Iimg = Iimg*0;
    
    % slice  parameters
    x0 = q(z_ind,1);
    y0 = q(z_ind,2);
    z0 = q(z_ind,3);
    NFP = q(z_ind,4);
    
    %sum and noise of image
    data_z0 = (data(:,:,z_ind)).*(data(:,:,z_ind)~=0);
    Nph = Nph_init(z_ind);
    std_Gbg = std_stack(z_ind);
    
    %% calc PSF
    
    % back focal plane phase 
    BFP_Phase = k*(x0*opt_phase_mat(:,:,1)+y0*opt_phase_mat(:,:,2)+z0.*opt_phase_mat(:,:,3)-NFP*opt_phase_mat(:,:,4));
    BFP_Phase = exp(1i.*(BFP_Phase));
    
    if vec_model_flag 
        g_bfp = g_bfp_init*0;
        g_img = g_bfp_init*0;
        
        if sum(IS.p_vec == 0) == 3 % freely rotating - superposition solution
            
            % BFP corrected intensity
            for g_id = 1:size(IS.g_bfp,3)
                g_bfp(:,:,g_id) = g_bfp_init(:,:,g_id) .*circ .*BFP_Phase;
            end
            normfact = 1;
            
            %calc back focal plane Green's tensor  
            for g_id = 1:size(IS.g_bfp,3)
                g_img(:,:,g_id) = 1/N.*(fft2(g_bfp(:,:,g_id)*normfact));
                Iimg = Iimg+g_img(:,:,g_id).*conj(g_img(:,:,g_id));
            end
        else % vectorial coherent dipole (with known orientation)
            % normalization
            for div_pol = 1:floor(size(IS.g_bfp,3)/3)
                for g_id = 1:3
                    g_bfp(:,:,g_id) = g_bfp_init(:,:,g_id).*BFP_Phase.*circ;
                    P_mat(:,:,g_id) = g_bfp(:,:,g_id).*IS.p_vec(g_id);
                end
%                 P_fact(div_pol) = sum(sum(abs(sum(P_mat,3)).^2));
            end
%             normfact = sqrt(Nph)./sqrt(sum(P_fact));
            normfact =1;
            for div_pol = 1:floor(size(IS.g_bfp,3)/3)
                for g_id = (1:3)+3*(div_pol-1)
                    g_img(:,:,g_id) = 1/N.*(fft2((g_bfp(:,:,g_id).*normfact)).*IS.p_vec(g_id-3*(div_pol-1)));
                end
                Iimg = Iimg+abs(sum(g_img,3)).^2;
            end
        end
        
    else
        
        % scalar bfp
        Ebfp = BFP_Phase.*circ_sc.*int_cos;

        % generate E field in image plane from pupil function
        Eimg = 1/N*fft2(Ebfp); % image plane
        Iimg = Eimg.*conj(Eimg); % intensity at image plane
    end
    %% shift and blur the PSF
    tmp  = fftshift(real(Iimg));
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
    
    %% calc gradient
    if vec_model_flag 
        if sum(IS.p_vec == 0) == 3 % freely rotating - superposition solution
            for g_id = 1:size(g_img,3)
                grad_tmp(:,:,g_id) = 2*1/N*real((fft2(ifftshift(data_z0).*1i.*conj(g_img(:,:,g_id))).*g_bfp_init(:,:,g_id))).*normfact; %%
            end
        else %  vectorial coherent dipole (with known orientation)
            for div_pol = 1:floor(size(IS.g_bfp,3)/3)
                Conj_fact = conj(sum(g_img(:,:,(1:3)+3*(div_pol-1)),3));
                for g_id = (1:3)+3*(div_pol-1)
                    grad_tmp(:,:,g_id) = 2*1/N*real(fft2(ifftshift(data_z0).*1i.*Conj_fact).*g_bfp_init(:,:,g_id)).*normfact; %%
                end
            end
        end
        grad(:,:,z_ind) = sum((grad_tmp),3);
    else  % scalar solution
        grad(:,:,z_ind) = 2*1/N*real(fft2(ifftshift(data_z0).*1i.*conj(Eimg)).*Ebfp).*normfact;
    end
    
end

% create a weigting to the SAF light

if IS.SAF_flag
    %     grad = sum(grad,3).*circ_gpu ;
    grad = mean(grad,3);
else
    grad = mean(grad,3).*circ_sc;
end










