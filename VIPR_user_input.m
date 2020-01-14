%% user defined flags
Prior_mask_flag = 0; % 1 loading a predefined phase mask for initialization
gpu_flag = 1; % 1 - use GPU, 0 - on CPU
vec_model_flag = 1; % 1 - vectorial model, 0 - scalar model
cost_function_flag = 4; % optimization cost 1 - L1, 2 - L2, 3 - Poiss MLE, 4 - Sum of gaussians MLE
plot_flag = 1; % plot while SGD runs, slows down ~ X4
Alg_flag = 1  ; % gradient method : 1 - ADAM, 2 - Adamax , 3- Nadam, 4 - Nesterov ,5- Vanilla SGD
vec_model_pol = 'y' ; %'x' or 'y' for having a  polarizer, 'b' for full vectorial
noisy_flag = 1; % 0- for design PSFs, 1 - for PR;

%% load data stack and z positions and define optical parameters
% input optical parameters


if data_set==1
    %optical parameters ( some parameters are a guess as they are not provided by the challange)
    IS = init_input_EPFL;
else
    %optical parameters 
    IS = init_input;
end 
%% optimization parameters
% the parameters in this section define the optimization proccess

% pre-proc parameters
IS.I_thr_flag = 2; % 1- thr above IS.thr*max(I) per image, else - thr above IS.thr*background_std
IS.I_thr = 2; %  threshold parameter
IS.corner_size = 10; % size of corners to  estimate noise [pixels]

% hyper-params
IS.SGDiter = 300; %  how  many iterations to SGD
IS.step_size = 3e-1; % step  size (try 3e-1 for ADAM and 3e-8 for SGD)
IS.point_num = 3; % size of mini-batch per SGD iteration

% additional options
est_gBlur_flag = 1; % 1- estimate gBlur after 1/3 of the iterations
IS.gBlur_cost = 2; % cost to estimate gBlur if est_gBlur_flag=1, (1-4 same as cost_function_flag , 5 - by corr2)
% option not to use SGD
IS.last_iter = 100; % how many iterations  to run not with SGD (at end of optimization)
IS.last_iter_flag = 1; % 1 - contuine SGD, 2 - global gradient, 3- batch the lowest correlation points, 4- adaptive sampling with side info on corr
IS.thr_corr = 0.01; % threshold for correlation calc (used if last_iter_flag = 3)
IS.upsample_fact = 1; %
IS.update_Signal = 1; % 1 - update signal at second half of iterations (needs more iterations, but is more accurate), 0 - keep the image sum 
% plot sizes for PSF
IS.plotsize = 99 ; % size of psf plots [pixels]

%% initialize xyz positions to emitter radius
if data_set == 1
    %load data
    EPFL_data_load;
    % how much to sample 
    dI = 1;
    % choose a single stack
    IMG_T = DH_PSF(:,:,1:dI:end,5);
    % take minus z to match to NFP 
    z_stack_pos = -z(1:dI:end,1)';
    % image centering 
    xy = xy(1:dI:end,:,1);
    IS.FOV_size = size(IMG_T,1);
    % do scalar model because we don't know really what the system is 
    vec_model_flag = 0;
else
    % load data and z positions
    load_data_stack;
    % xy positions set to zero 
    xy = zeros(length(z_stack_pos),2);
end
z_pos = zeros(length(z_stack_pos),1)+IS.z_emit;

%% full coordinate definition
q_cord = [xy';z_pos';z_stack_pos]'; %x,y,z,NFP

%% add prior mask if needed
if Prior_mask_flag == 1
    % create phase mask here which fits the aperture (= aperture)
%     maskInit = ;
end

