function [Iimg] = PSF_generator(phase_mask,xyz,IS,nomFocus,opt_phase_mat,g_bfp,circ,circ_sc,int_cos,vec_model_flag)
%% inputs
% phase_mask - phase mask design in BFP
% xyz - vector 1X3 of cartesian emitter position [um]
% IS - input structure defined in the main function by init_input
% nomFocus - nominal focal plane position
% opt_phase_mat,normfact_vec_gpu,g_bfp_gpu,circ_gpu,circ_gpu_sc,int_cos - optical matrices defined by create_gpu_mats
% vec_model_flag - use vectorial (1) or scalar models (0)

%% outputs
% Iimg - image intensity pattern

%% ESTABLISH COORDINATE SYSTEMS - bfp and CCD plane

Iimg = phase_mask*0;
N = size(IS.circmask,1);

% get cropping coordinates
N_crop = IS.FOV_size;

%wave number
k = 2 * pi/IS.lambda;

% BFP phase component 
BFP_Phase = k*(xyz(3).*opt_phase_mat(:,:,3)+opt_phase_mat(:,:,1)*xyz(1) + ...
    opt_phase_mat(:,:,2).*xyz(2)+opt_phase_mat(:,:,4).*(-nomFocus));
BFP_Phase(isnan(BFP_Phase))=0;
%add phase mask
BFP_Phase = exp(1i.*(BFP_Phase+phase_mask));

if vec_model_flag == 1
    g_img = opt_phase_mat*0;
    
    % calc image pattern
    if sum(abs(IS.p_vec) == 0) == 3 % freely rotating - superposition solution

        for g_id = 1:size(IS.g_bfp,3)
            g_bfp(:,:,g_id) = g_bfp(:,:,g_id).*BFP_Phase.*circ;
            g_img(:,:,g_id) = 1/N.*fftshift(fft2(ifftshift(g_bfp(:,:,g_id))));
            Iimg = Iimg+g_img(:,:,g_id).*conj(g_img(:,:,g_id));
        end
        
    else
        
        % calc diffraction pattern
        for div_pol = 1:floor(size(IS.g_bfp,3)/3)
            for g_id = (1:3)+3*(div_pol-1)
                g_bfp(:,:,g_id) = g_bfp(:,:,g_id).*BFP_Phase.*circ;
                g_img(:,:,g_id) = 1/N.*fftshift(fft2(ifftshift(g_bfp(:,:,g_id)))).*IS.p_vec(g_id-3*(div_pol-1));
            end
            Iimg = Iimg+abs(sum(g_img,3)).^2;
        end
        
    end
else
    circ_sc = circ_sc.*int_cos;
    %
    Pupil = BFP_Phase.*circ_sc;
    
    % normalization
    normfact = 1;
    Pupil = Pupil*normfact;
    
    % generate image from pupil function
    Eimg_OS = 1/N.*fftshift(fft2(ifftshift(Pupil))); % oversampled image plane guess
    Iimg = Eimg_OS.*conj(Eimg_OS); %
    
end
% crop image
if mod(size(Iimg,1),2)==1
    if mod((N_crop)/2,1)
        Iimg = real(Iimg(end/2-floor(N_crop)/2+1:end/2+floor(N_crop)/2,end/2-floor(N_crop)/2+1:end/2+floor(N_crop)/2));
    else
        Iimg = real(Iimg(end/2-floor(N_crop)/2+0.5+1:end/2+floor(N_crop)/2+1,end/2-floor(N_crop)/2+0.5+1:end/2+floor(N_crop)/2+1));
    end
else
    Iimg = real(Iimg(round(size(Iimg,1)/2-round(N_crop-1)/2+1/2):round(size(Iimg,1)/2+floor(N_crop-1)/2+1/2),round(size(Iimg,2)/2-round(N_crop-1)/2+1/2):round(size(Iimg,2)/2+floor(N_crop-1)/2+1/2)));
end

% rescale in image space
Iimg = Iimg.*IS.Signal./sum(Iimg(:));

%  blur
if IS.gBlur > 0
    h = (fspecial('gaussian',[5 5],IS.gBlur));
    if isa(Iimg,'gpuArray')
        h = gpuArray(h);
    end
    Iimg = imfilter(Iimg,h,'replicate');
end

% add bg
Iimg = (Iimg)+(IS.bg);








