function [y_est,y_hat] = BM3D_denoising(z, sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BM3D is algorithm for attenuation of additive white Gaussian noise from 
%  grayscale images. This algorithm reproduces the results from the article:
%
%  [1] K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image denoising 
%   by sparse 3D transform-domain collaborative filtering," 
%   IEEE Transactions on Image Processing, December, 2006, in review,
%   preprint at http://www.cs.tut.fi/~foi/GCF-BM3D.
%
%  FUNCTION INTERFACE:
%
%  [PSNR, y_est] = BM3D(y, z, sigma, profile, print_to_screen)
%
%  ! The function can work without any of the input arguments, 
%   in which case, the internal default ones are used !
% 
%  BASIC USAGE EXAMPLES:
%
%     Case 1) Using the default parameters (i.e., image name, sigma, etc.)
% 
%      [PSNR, y_est] = BM3D;
% 
%     Case 2) Using an external noisy image:
%
%      % Read a grayscale image and scale its intensities in range [0,1]
%      y = im2double(imread('Cameraman256.png')); 
%      % Generate the same seed used in the experimental results of [1]
%      randn('seed', 0);
%      % Standard deviation of the noise --- corresponding to intensity 
%      %  range [0,255], despite that the input was scaled in [0,1]
%      sigma = 25;
%      % Add the AWGN with zero mean and standard deviation 'sigma'
%      z = y + (sigma/255)*randn(size(y));
%      % Denoise 'z'. The denoised image is 'y_est', and 'NA = 1' because 
%      %  the true image was not provided
%      [NA, y_est] = BM3D(1, z, sigma); 
%      % Compute the putput PSNR
%      PSNR = 10*log10(1/mean((y(:)-y_est(:)).^2))
%      % show the noisy image 'z' and the denoised 'y_est'
%      figure; imshow(z);   
%      figure; imshow(y_est);
% 
%     Case 3) If the original image y is provided as the first input 
%      argument, then some additional information is printed (PSNRs, 
%      figures, etc.). That is, "[NA, y_est] = BM3D(1, z, sigma);" in the
%      above code should be replaced with:
% 
%      [PSNR, y_est] = BM3D(y, z, sigma);
% 
% 
%  INPUT ARGUMENTS (OPTIONAL):
%
%     1) y (matrix M x N): Noise-free image (needed for computing PSNR),
%                           replace with the scalar 1 if not available.
%     2) z (matrix M x N): Noisy image (intensities in range [0,1] or [0,255])
%     3) sigma (double)  : Std. dev. of the noise (corresponding to intensities
%                          in range [0,255] even if the range of z is [0,1])
%     4) profile (char)  : 'np' --> Normal Profile 
%                          'lc' --> Fast Profile
%     5) print_to_screen : 0 --> do not print output information (and do 
%                                not plot figures)
%                          1 --> print information and plot figures
%
%  OUTPUTS:
%     1) PSNR (double)          : Output PSNR (dB), only if the original 
%                                 image is available, otherwise PSNR = 0                                               
%     2) y_est (matrix M x N): Final estimate (in the range [0,1])
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright © 2006 Tampere University of Technology. All rights reserved.
% This work should only be used for nonprofit purposes.
%
% AUTHORS:
%     Kostadin Dabov, email: dabov _at_ cs.tut.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% In case, a noisy image z is not provided, then use the filename 
%%%%  below to read an original image (might contain path also). Later, 
%%%%  artificial AWGN noise is added and this noisy image is processed 
%%%%  by the BM3D.
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Quality/complexity trade-off profile selection
%%%%
%%%%  'np' --> Normal Profile (balanced quality)
%%%%  'lc' --> Low Complexity Profile (fast, lower quality)
%%%%
%%%%  'high' --> High Profile (high quality, not documented in [1])
%%%%
if (exist('profile') ~= 1)
    profile         = 'high'; %% default profile
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Specify the std. dev. of the corrupting noise
%%%%
if (exist('sigma') ~= 1),
    sigma               = 25; %% default standard deviation of the AWGN
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Following are the parameters for the Normal Profile.
%%%%

%%%% Select transforms ('dct', 'dst', 'hadamard', or anything that is listed by 'help wfilters'):
transform_2D_HT_name     = 'bior1.5'; %% transform used for the HT filt. of size N1 x N1
transform_2D_Wiener_name = 'dct';     %% transform used for the Wiener filt. of size N1_wiener x N1_wiener
transform_3rd_dim_name   = 'haar';    %% transform used in the 3-rd dim, the same for HT and Wiener filt.

%%%% Hard-thresholding (HT) parameters:
N1                  = 8;   %% N1 x N1 is the block size used for the hard-thresholding (HT) filtering
Nstep               = 3;   %% sliding step to process every next reference block
N2                  = 16;  %% maximum number of similar blocks (maximum size of the 3rd dimension of a 3D array)
Ns                  = 39;  %% length of the side of the search neighborhood for full-search block-matching (BM), must be odd
tau_match           = 3000;%% threshold for the block-distance (d-distance)
lambda_thr2D        = 0;   %% threshold parameter for the coarse initial denoising used in the d-distance measure
lambda_thr3D        = 2.7; %% threshold parameter for the hard-thresholding in 3D transform domain
beta                = 2.0; %% parameter of the 2D Kaiser window used in the reconstruction

%%%% Wiener filtering parameters:
N1_wiener           = 8;
Nstep_wiener        = 3;
N2_wiener           = 32;
Ns_wiener           = 39;
tau_match_wiener    = 400;
beta_wiener         = 2.0;

%%%% Block-matching parameters:
stepFS              = 1;  %% step that forces to switch to full-search BM, "1" implies always full-search
smallLN             = 'not used in np'; %% if stepFS > 1, then this specifies the size of the small local search neighb.
stepFSW             = 1;
smallLNW            = 'not used in np';
thrToIncStep        = 8;  %% used in the HT filtering to increase the sliding step in uniform regions

if strcmp(profile, 'lc') == 1,

    Nstep               = 6;
    Ns                  = 25;
    Nstep_wiener        = 5;
    N2_wiener           = 16;
    Ns_wiener           = 25;

    thrToIncStep        = 3;
    smallLN             = 3;
    stepFS              = 6*Nstep;
    smallLNW            = 2;
    stepFSW             = 5*Nstep_wiener;

end

if (strcmp(profile, 'vn') == 1) | (sigma > 40),

    transform_2D_HT_name = 'dct'; 
    
    N1                  = 12;
    Nstep               = 4;
 
    N1_wiener           = 11;
    Nstep_wiener        = 6;

    lambda_thr3D        = 2.8;
    lambda_thr2D        = 2.0;
    thrToIncStep        = 3;
    tau_match_wiener    = 3500;
    tau_match           = 5000;
    
    Ns_wiener           = 39;
    
end

decLevel = 0;        %% dec. levels of the dyadic wavelet 2D transform for blocks (0 means full decomposition, higher values decrease the dec. number)
thr_mask = ones(N1); %% N1xN1 mask of threshold scaling coeff. --- by default there is no scaling, however the use of different thresholds for different wavelet decompoistion subbands can be done with this matrix

if strcmp(profile, 'high') == 1, %% this profile is not documented in [1]
    
    decLevel     = 1; 
    Nstep        = 2;
    Nstep_wiener = 2;
    lambda_thr3D = 2.5;
    vMask = ones(N1,1); vMask((end/4+1):end/2)= 1.01; vMask((end/2+1):end) = 1.07; %% this allows to have different threhsolds for the finest and next-to-the-finest subbands
    thr_mask = vMask * vMask'; 
    beta         = 2.5;
    beta_wiener  = 1.5;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Note: touch below this point only if you know what you are doing!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Check whether to dump information to the screen or remain silent
dump_output_information = 1;
if (exist('print_to_screen') == 1) & (print_to_screen == 0),
    dump_output_information = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create transform matrices, etc.
%%%%
[Tfor, Tinv]   = getTransfMatrix(N1, transform_2D_HT_name, decLevel);     %% get (normalized) forward and inverse transform matrices
[TforW, TinvW] = getTransfMatrix(N1_wiener, transform_2D_Wiener_name, 0); %% get (normalized) forward and inverse transform matrices

if (strcmp(transform_3rd_dim_name, 'haar') == 1) | (strcmp(transform_3rd_dim_name(end-2:end), '1.1') == 1),
    %%% If Haar is used in the 3-rd dimension, then a fast internal transform is used, thus no need to generate transform
    %%% matrices.
    hadper_trans_single_den         = {};
    inverse_hadper_trans_single_den = {};
else
    %%% Create transform matrices. The transforms are later applied by
    %%% matrix-vector multiplication for the 1D case.
    for hpow = 0:ceil(log2(max(N2,N2_wiener))),
        h = 2^hpow;
        [Tfor3rd, Tinv3rd]   = getTransfMatrix(h, transform_3rd_dim_name, 0);
        hadper_trans_single_den{h}         = single(Tfor3rd);
        inverse_hadper_trans_single_den{h} = single(Tinv3rd');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2D Kaiser windows used in the aggregation of block-wise estimates
%%%%
Wwin2D           = kaiser(N1, beta) * kaiser(N1, beta)'; % Kaiser window used in the aggregation of the HT part
Wwin2D_wiener    = kaiser(N1_wiener, beta_wiener) * kaiser(N1_wiener, beta_wiener)'; % Kaiser window used in the aggregation of the Wiener filt. part
%Wwin2D=ones(N1);Wwin2D_wiener=ones(N1_wiener);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% If needed, read images, generate noise, or scale the images to the 
%%%% [0,1] interval
%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 1. Produce the basic estimate by HT filtering
%%%%
tic;
y_hat = bm3d_thr(z, hadper_trans_single_den, Nstep, N1, N2, lambda_thr2D,...
	lambda_thr3D, tau_match*N1*N1/(255*255), (Ns-1)/2, (sigma/255), thrToIncStep, single(Tfor), single(Tinv)', inverse_hadper_trans_single_den, single(thr_mask), Wwin2D, smallLN, stepFS );
estimate_elapsed_time = toc;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 2. Produce the final estimate by Wiener filtering (using the 
%%%%  hard-thresholding initial estimate)
%%%%
tic;
y_est = bm3d_wiener(z, y_hat, hadper_trans_single_den, Nstep_wiener, N1_wiener, N2_wiener, ...
    'unused arg', tau_match_wiener*N1_wiener*N1_wiener/(255*255), (Ns_wiener-1)/2, (sigma/255), 'unused arg', single(TforW), single(TinvW)', inverse_hadper_trans_single_den, Wwin2D_wiener, smallLNW, stepFSW, single(ones(N1_wiener)) );
wiener_elapsed_time = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate the final estimate's PSNR, print it, and show the
%%%% denoised image next to the noisy one
%%%%
y_est = double(y_est);


return;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some auxiliary functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [Transf_Matrix, invTransf_Matrix] = getTransfMatrix (N, transform_type, Nden)
%
% Create forward and inverse transform matrices, which allow for perfect
% reconstruction. The forward transform matrix is normalized so that the 
% l2-norm of each basis element is 1.
%
% [Tforward, Tinverse] = getTransfMatrix (N, transform_type, dec_levels)
%
%  INPUTS:
%
%   N           --> Size of the transform (for wavelets, must be 2^K)
%
%   transform_type  --> 'dct', 'dst', 'hadamard', or anything that is 
%                       listed by 'help wfilters' (bi-orthogonal wavelets)
%                       'DCrand' -- an orthonormal transform with a DC and all
%                       the other basis elements of random nature
%
%   dec_levels      --> If a wavelet transform is generated, this is the
%                       desired decomposition level. Must be in the
%                       range [0, log2(N)-1], where "0" implies
%                       full decomposition.
%
%  OUTPUTS:
%
%   Tforward        --> (N x N) Forward transform matrix
%
%   Tforward        --> (N x N) Inverse transform matrix
%

if exist('Nden') ~= 1,
    Nden = 0;
end

if N == 1,
    Transf_Matrix = 1;
elseif strcmp(transform_type, 'dct') == 1,
    Transf_Matrix    = dct(eye(N));
elseif strcmp(transform_type, 'dst') == 1,
    Transf_Matrix    = dst(eye(N));
elseif strcmp(transform_type, 'DCrand') == 1,
    x = randn(N); x(1:end,1) = 1; [Q,R] = qr(x); 
    if (Q(1) < 0), 
        Q = -Q; 
    end;
    Transf_Matrix = Q';
elseif strcmp(transform_type, 'hadamard') == 1,
    Transf_Matrix    = hadamard(N);
else %% wavelet transform
   
    %%% Set periodic boundary conditions, to preserve bi-orthogonality
    dwtmode('per','nodisp');  

    [LO_D, HI_D, LO_R, HI_R] = wfilters(transform_type);
    for i = 1:N
        Transf_Matrix(i,:)=waverec(circshift([1 zeros(1,N-1)],[Nden i-1]), 2.^[Nden Nden:log2(N)], LO_D, -HI_D);  %% construct transform matrix
    end
end

%%% Normalize the basis elements
Transf_Matrix = (Transf_Matrix' * diag(sqrt(1./sum(Transf_Matrix.^2,2))))'; 

%%% Compute the inverse transform matrix
invTransf_Matrix = inv(Transf_Matrix);

return;

