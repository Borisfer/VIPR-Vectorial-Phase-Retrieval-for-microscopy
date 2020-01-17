function IS = init_input_EPFL

% defined signal modifiers ( don't change)
IS.Signal=1;
IS.bg=0;
% 
IS.M=100; % objective magnification
IS.NA=1.49; % objective NA 
IS.lambda=[670]*10^-3; % wavelength [um]
IS.SLM_psize = 30; % pixel size of SLM [um] 
IS.Cam_psize = 10; %pixel size of camera [um]
IS.gBlur=1; %initial guess of the blurring 
IS.n_glass=1.518; % RI of immersion oil 
IS.n_med=1.33; % RI of sample medium 
IS.f_4f = 10e4; % focal length of 4-f  system (if 2D imaging- use tube f)
IS.FOV_size = 80;
IS.SAF_flag = 1; % include SAF or not (1=include - recommended)
IS.Iris_pNA = 1; % optional iris to limit BFP, in range [0,1] where 1 is the NA 
% emitter size
IS.z_emit = 0.05; % emitter radius [um]
% IS.z_emit = 0.02; % emitter radius [um]

% polarization of dipole (0,0,0) is for freely rotating
IS.p_vec = [0,0,0]; % incoherent

% for coherent - normalize  the vector 
% IS.p_vec = [1,-1,1];
% IS.p_vec = [1,0,0];
% IS.p_vec = IS.p_vec./norm(IS.p_vec);




