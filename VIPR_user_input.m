%% user defined flags
gpu_flag = 1; % 1 - use GPU, 0 - on CPU
vec_model_flag = 1; % 1 - vectorial model, 0 - scalar model
cost_function_flag = 4; % optimization cost 1 - L1, 2 - L2, 3 - Poiss MLE, 4 - Sum of gaussians MLE
plot_flag = 1; % plot while SGD runs, slows down ~ X4
Alg_flag = 1  ; % gradient method : 1 - ADAM, 2 - Adamax , 3- Nadam, 4 - Nesterov ,5- Vanilla SGD
vec_model_pol = 'y' ; %'x' or 'y' for having a  polarizer, 'b' for full vectorial
noisy_flag = 1; % 0- for design PSFs, 1 - for PR;
est_gBlur_flag = 1; % 1- estimate gBlur after 1/3 of the iterations

%% define optical parameters
IS.Signal=1; % for generator - keep 1
IS.bg=0; % for generator - keep 0

% input optical parameters
if data_set==1
     %optical parameters
    %
    IS.M=100; % objective magnification
    IS.NA=1.45; % objective NA
    IS.lambda=[605]*10^-3; % wavelength [um]
    IS.SLM_psize = 20; % pixel size of SLM [um]
    IS.Cam_psize = 16; % pixel size of camera [um]
    IS.gBlur=0.75; %initial guess of the blurring
    IS.n_glass=1.518; % RI of immersion oil
    IS.n_med=1.33; % RI of sample medium
    IS.f_4f = 15e4; % focal length of 4-f  system (if 2D imaging- use tube f)
    IS.FOV_size = 80; % size of ROI used
    IS.SAF_flag = 1; % include SAF or not (1=include - recommended)
    IS.Iris_pNA = 1; % optional iris to limit BFP, in range [0,1] where 1 is the NA
    % emitter size
    IS.z_emit = 0.023; % emitter radius [um]
    
    % polarization of dipole (0,0,0) is for freely rotating
    IS.p_vec = [0,0,0]; % incoherent
    
    % for coherent - normalize  the vector
%     IS.p_vec = [1,-1,1];
    % IS.p_vec = [1,0,0];
elseif data_set == 2
    IS = init_input;
    
elseif data_set == 3 % 
    IS = init_input_EPFL;
end

if sum(IS.p_vec == 0) ~= 3
    IS.p_vec = IS.p_vec./norm(IS.p_vec);
end

%% optimization parameters
% the parameters in this section define the optimization proccess

% pre-proc parameters
IS.I_thr_flag = 2; % 1- thr above IS.thr*max(I) per image, else - thr above IS.thr*background_std
IS.I_thr = 1; %  threshold parameter
IS.corner_size = 10; % size of corners to  estimate noise [pixels]

% hyper-params
IS.SGDiter = 250; %  how  many iterations to SGD
IS.step_size = 5e-1; % step  size (try 3e-1 for ADAM and 3e-8 for SGD)
IS.point_num = 3; % size of mini-batch per SGD iteration

% additional options
IS.gBlur_cost = 2; % cost to estimate gBlur if est_gBlur_flag=1, (1-4 same as cost_function_flag , 5 - by corr2)
% option not to use SGD
IS.last_iter = 50; % how many iterations  to run not with SGD (at end of optimization)
IS.last_iter_flag = 3; % 1 - contuine SGD, 2 - global gradient, 3- batch the lowest correlation points, 4- adaptive sampling with side info on corr
IS.thr_corr = 0.01; % threshold for correlation calc (used if last_iter_flag = 3)
IS.upsample_fact = 1; % how much to upsample the data (usually leave at 1)
IS.update_Signal = 0; % 1 - update signal at second half of iterations (needs more iterations, but is more accurate), 0 - keep the image sum
IS.Photobleach_cost_mult = 0; % add to cost the SNR consideration
% plot sizes for PSF
IS.plotsize = 99 ; % size of psf plots [pixels]
