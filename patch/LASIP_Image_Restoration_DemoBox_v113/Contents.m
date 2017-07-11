% Anisotropic Nonparametric Image Restoration DemoBox (for Matlab 6.5 or later)
% Version 1.13, 3 July 2012 - Tampere University of Technology
%
%
% Image Denoising.
%   demo_DenoisingGaussian              - Anisotropic LPA-ICI Denoising demo (AWGN).
%   demo_RecursiveDenoisingGaussian     - Recursive Anisotropic LPA-ICI Denoising demo (AWGN).
%   demo_DenoisingSignDepNoise          - Anisotropic LPA-ICI Denoising demo (Signal-Dependent Noise).
%   demo_MR_FilteringGaussian           - Anisotropic Multiresolution LPA denoising demo (AWGN).
%
% Inverse Filtering.
%   demo_DeblurringGaussian             - Anisotropic LPA-ICI Regularized Deconvolution (AWGN).
%   demo_DeblurringPoisson              - Anisotropic LPA-ICI Regularized Deconvolution (Poissonian Noise).
%   demo_InverseHalftoning              - Anisotropic LPA-ICI Inverse Halftoning (Error-diffusion).
%
% Differentiation.
%   demo_AnisotropicGradient            - Anisotropic Gradient demo (Riemann Surface).
%
% Adaptive Scale Selection.
%   function_ICI                        - Intersection of Confidence Intervals (ICI) algorithm.
%
% Kernel design.
%   demo_CreateLPAKernels               - Creates some LPA kernels.
%   demo_CreateMRLPAKernels             - Creates and displays some Multiresolution LPA kernels.
%   function_CreateLPAKernels           - Creates LPA kernels.
%
% Visualization.
%   function_AnisSect_explorer          - Displays currently used anisotropic estimation neighborhoods
%   utility_DrawLPAKernels              - Displays currently used LPA kernels.
%
% Other functions.
%   function_Errors                     - Computes Error-criteria (SNR, MSE, etc.).
%   function_stdEst2D                   - Estimates noise standard deviation (AWGN).
%   function_ErrorDiffusion             - Generates Error-diffusion Halftones (e.g. Floyd-Steinberg, Jarvis).

% Copyright 2003-2012 Tampere University of Technology, Tampere, Finland
% 