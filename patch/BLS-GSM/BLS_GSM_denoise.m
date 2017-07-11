function im_d=BLS_GSM_denoise(im,sig)
[Ny,Nx] = size(im);


% Noise Parameters (assumed additive, Gaussian and independent of the
% original image)
%sig = 20;		        % standard deviation
PS = ones(size(im));	% power spectral density (in this case, flat, i.e., white noise)
seed = 0;               % random seed

% Pyramidal representation parameters
Nsc = ceil(log2(min(Ny,Nx)) - 4);  % Number of scales (adapted to the image size)
% if opt==1
% Nor = 3;				            % Number of orientations (for X-Y separable wavelets it can only be 3)
% repres1 = 'uw';                     % Type of pyramid (shift-invariant version of an orthogonal wavelet, in this case)
% repres2 = 'daub1';                  % Type of wavelet (daubechies wavelet, order 2N, for 'daubN'; in this case, 'Haar')
% else
Nor = 8;                           % 8 orientations
repres1 = 'fs';                    % Full Steerable Pyramid, 5 scales for 512x512
repres2 = '';                      % Dummy parameter when using repres1 = 'fs'   
%end
% Model parameters (optimized: do not change them unless you are an advanced user with a deep understanding of the theory)
blSize = [5 5];	    % n x n coefficient neighborhood of spatial neighbors within the same subband
                    % (n must be odd): 
parent = 0;			% including or not (1/0) in the neighborhood a coefficient from the same spatial location
                    % and orientation as the central coefficient of the n x n neighborhood, but
                    % next coarser scale. Many times helps, but not always.
boundary = 1;		% Boundary mirror extension, to avoid boundary artifacts 
covariance = 1;     % Full covariance matrix (1) or only diagonal elements (0).
optim = 1;          % Bayes Least Squares solution (1), or MAP-Wiener solution in two steps (0)

% Uncomment the following 4 code lines for reproducing the results of our IEEE Trans. on Im. Proc., Nov. 2003 paper
% This configuration is slower than the previous one, but it gives slightly better results (SNR)
% on average for the test images "lena", "barbara", and "boats" used in the cited article.


parent = 1;                        % Include a parent in the neighborhood


% Call the denoising function
tic; im_d = denoi_BLS_GSM(im, sig, PS, blSize, parent, boundary, Nsc, Nor, covariance, optim, repres1, repres2, seed); toc

