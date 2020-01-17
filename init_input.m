function IS = init_input

% defined signal modifiers ( don't change)
IS.Signal=1;
IS.bg=0;
% 
IS.M=100; % objective magnification
IS.NA=1.45; % objective NA 
IS.lambda=[605]*10^-3; % wavelength [um]
IS.SLM_psize = 20; % pixel size of SLM [um] 
IS.Cam_psize = 16; % pixel size of camera [um]
IS.gBlur=0.6; %initial guess of the blurring 
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
% IS.p_vec = [1,0,0];
% IS.p_vec = IS.p_vec./norm(IS.p_vec);




