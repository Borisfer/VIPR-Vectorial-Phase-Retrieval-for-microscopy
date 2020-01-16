function IS = pre_comp_opt(IS,vec_model_flag,pol)

% function calculates all the optical matrices that can be pre-computed 
%
% Inputs
% IS          -   Input structure including all the optical parameters
% vec_model_flag -   1 to use the vectorial model, 0- use the scalar model
% pol         -   polarization direction for the vectorial model
%               
% Outputs
% IS          -   Updated input structure. 
%

% physical size of back focal plane objects: mask and E field
maskdiam_px = IS.maskDiam_m/IS.SLM_psize;
BFPdiam_m = 2*IS.f_4f*IS.NA/sqrt(IS.M^2 - IS.NA^2); % diameter of E field in bfp (region of E field support) [um]
xPhys = ((1:maskdiam_px) - ceil(maskdiam_px/2))*IS.SLM_psize;
[IS.XI,IS.ETA] = meshgrid(xPhys,xPhys); % physical coordinates in SLM space [um]

% establish angular coordinates in back focal plane
xAng = linspace(-1,1,maskdiam_px)*maskdiam_px/(BFPdiam_m/IS.SLM_psize);
[XX,YY] = meshgrid(xAng,xAng); % each pixel is (BFPdiam_px/2) wide
r = sqrt(XX.^2+YY.^2); % radial coordinate s.t. r = 1 at edge of E field support
mask_NA = r<=IS.Iris_pNA;
% calculate angles in bfp for both refractive indices
IS.XX = XX;
IS.YY = YY;

% wave propagation in glass
IS.sin_theta = IS.NA/IS.n_glass*r;
IS.cos_theta = sqrt(1-IS.sin_theta.^2);
% circmask_glass = ~(abs(imag(IS.cos_theta))>0);

% wave propagation in medium
sin_theta2 = IS.n_glass/IS.n_med*IS.sin_theta; %Snell's law
IS.cos_theta2 = sqrt(1-sin_theta2.^2);
circmask_med = ~(abs(imag(IS.cos_theta2))>0);

% circmask - only UAF, circmask_large - UAF and SAF 
IS.circmask = (mask_NA).* circmask_med;
IS.circmask_large =  (mask_NA);

% amplitude correction
IS.int_cos = 1./sqrt(IS.cos_theta);
% IS.int_cos = 1;

if IS.SAF_flag == 0
    tmp = IS.circmask;
else
    tmp = IS.circmask_large;
end
% phase matrices 
IS.psi(:,:,1) =  IS.XI*IS.M/IS.f_4f; % x shift phase
IS.psi(:,:,2) =  IS.ETA*IS.M/IS.f_4f; % y shift phase
IS.psi(:,:,3) =  IS.n_med*IS.cos_theta2; % depth phase
IS.psi(:,:,4) =  IS.n_glass*IS.cos_theta; % defocus phase
% wave-number 
k = 2*pi/IS.lambda;

if vec_model_flag
    % fresnel coeffs
    IS.tp = (2 * IS.n_med .* IS.cos_theta2) ./ (IS.n_glass.*IS.cos_theta2 +IS.n_med .* IS.cos_theta);
    IS.ts = (2 * IS.n_med .* IS.cos_theta2) ./ (IS.n_med.*IS.cos_theta2 + IS.n_glass .* IS.cos_theta);

    % transfer coeficients
    IS.c1  = (IS.n_glass/IS.n_med)^2 .* (IS.cos_theta./IS.cos_theta2) .* IS.tp ;
    IS.c2 =  (IS.n_glass/IS.n_med) .* IS.tp ;
    IS.c3 =  (IS.n_glass/IS.n_med) .* (IS.cos_theta./IS.cos_theta2) .* IS.ts ;
    IS.c1(isnan(IS.c1)) = 0;
    IS.c2(isnan(IS.c2)) = 0;
    IS.c3(isnan(IS.c3)) = 0;
    
    % BFP angular coords
    sin_phi = (IS.YY ./sqrt(IS.XX.^2+IS.YY.^2));
    cos_phi = (IS.XX ./sqrt(IS.XX.^2+IS.YY.^2));
    
    if mod(size(sin_phi,1),2)==1
        sin_phi(ceil(end/2),ceil(end/2))=1;
        cos_phi(ceil(end/2),ceil(end/2))=1;
    end
    
    
    %% green matrix BFP
    IS.g_bfp_xx = (IS.c3 .* (sin_phi.^2) + IS.c2.* (cos_phi.^2).* IS.cos_theta ) ;
    IS.g_bfp_xy = (- sin_phi .* cos_phi .*(IS.c3-IS.c2.*IS.cos_theta));
    IS.g_bfp_xz = (-IS.c1 .* cos_phi .* IS.sin_theta );
    IS.g_bfp_yx = (-  sin_phi .* cos_phi.*(IS.c3-IS.c2.*IS.cos_theta) );
    IS.g_bfp_yy = (IS.c3.* (cos_phi.^2)+IS.c2.* (sin_phi.^2).*IS.cos_theta);
    IS.g_bfp_yz = (-IS.c1 .* sin_phi .* IS.sin_theta );
    % SAF decay calc 
    IS.BFP_decay = abs(exp(1i.*k.*IS.z_emit.*IS.psi(:,:,3))).*tmp;
    %     IS.BFP_decay = 1;
    
    if strcmp(pol,'x')
        IS.g_bfp(:,:,1) =  IS.int_cos.*IS.g_bfp_xx.*tmp;
        IS.g_bfp(:,:,2) =  IS.int_cos.*IS.g_bfp_xy.*tmp;
        IS.g_bfp(:,:,3) =  IS.int_cos.*IS.g_bfp_xz.*tmp;
        
    elseif  strcmp(pol,'y')
        IS.g_bfp(:,:,1) =  IS.int_cos.*IS.g_bfp_yx.*tmp;
        IS.g_bfp(:,:,2) =  IS.int_cos.*IS.g_bfp_yy.*tmp;
        IS.g_bfp(:,:,3) =  IS.int_cos.*IS.g_bfp_yz.*tmp;
        
    else
        IS.g_bfp(:,:,1) =  IS.int_cos.*IS.g_bfp_xx.*tmp;
        IS.g_bfp(:,:,2) =  IS.int_cos.*IS.g_bfp_xy.*tmp;
        IS.g_bfp(:,:,3) =  IS.int_cos.*IS.g_bfp_xz.*tmp;
        IS.g_bfp(:,:,4) =  IS.int_cos.*IS.g_bfp_yx.*tmp;
        IS.g_bfp(:,:,5) =  IS.int_cos.*IS.g_bfp_yy.*tmp;
        IS.g_bfp(:,:,6) =  IS.int_cos.*IS.g_bfp_yz.*tmp;
    end
        
    %% vector model normalization factor in fourier space (Parseval)
    if sum(IS.p_vec)==0 % free rotating
        for g_id = 1:size(IS.g_bfp,3)
            P_fact(g_id) = sum(sum(abs(IS.g_bfp(:,:,g_id).*IS.BFP_decay).^2));
        end
        normfact = sqrt((IS.Signal))./sqrt(sum(P_fact));
    else
        for div_pol = 1:floor(size(IS.g_bfp,3)/3)
            for g_id = 1:3
                P_mat(:,:,g_id) = IS.g_bfp(:,:,g_id).*IS.BFP_decay.*IS.p_vec(g_id);
            end
            P_fact(div_pol) = sum(sum(abs(sum(P_mat,3)).^2));
        end
        normfact = sqrt(gpuArray(IS.Signal))./sqrt(sum(P_fact));
    end
    % circmask
    IS.circmask_opt = tmp;
    IS.normfact = normfact;
end

%% scalar model normalization factor in fourier space (Parseval)
N = size(IS.circmask,1);
tmp = IS.circmask;
IS.circmask_opt_sc = tmp ;
IS.circmask_opt = IS.circmask_large ;

Q = tmp.*IS.int_cos;
IS.Parseval_fact = 1/sqrt((sum(Q(:).^2)*(N^2))); % for scalar model
%
