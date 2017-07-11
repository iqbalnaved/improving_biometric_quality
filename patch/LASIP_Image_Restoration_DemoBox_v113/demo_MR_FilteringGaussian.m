% Anisotropic multiresolution (MR) LPA Denoising Demo (demo_MR_FilteringGaussian)
% based on the MR analysis, thresholding and synthesis.
%
% Alessandro Foi, Vladimir Katkovnik - Tampere University of Technology - 2003-2005
% -----------------------------------------------------------------------------------
%
% Performs the MR anisotropic LPA denoising on observations which are
% contaminated by additive Gaussian White noise.
% Multiscale kernels are used for MR signal analysis and thresholding for noise
% removal.
%
%
% Observation model:
%
%    z = y + n
%
% z : noisy observation
% y : true image (assumed as unknown)
% n : Gaussian white noise
%
% Other key variables:
% y_hat   : anisotropic "fused" estimate
%
%
% Katkovnik, V., “Multiresolution nonparametric regression: a new approach to pointwise spatial adaptation,”
% Digital Signal Processing vol. 15, pp. 73–116, 2005.
%

clear all, close all, global z y lenh sigma h1





addnoise=1;                % add noise to observation
sigma_noise=0.1;           % standard deviation of the noise


%--------------------------------------------------------------------------
% LPA ORDER, KERNEL SIZES, THRESHOLDS
%--------------------------------------------------------------------------
m=[0,0];        % THE VECTOR ORDER OF LPA;

ndir=8;  % number of directions
h1=[1 3 5 11 17];  % LPA KERNELS SIZES
h2=ones(size(h1));
lenh=length(h1);
maxh1=max(h1);

thresholding_type=31;
threshold_parameter=1.3;

q=0.2;   % Larger q decreases the effective threshold and makes the estimate noisy

%%% The following threshold types are used with recommended threshold parameter values:
%%%  1 HARD THRESHOLDING, threshold_parameter=2.5;
%%%  2 SOFT THRESHOLDING, threshold_parameter=2;
%%%  3 STEIN rule, threshold_parameter=2;
%%% 31 - SMOOTHED STEIN rule, threshold_parameter=1.3;
%%% 21 - ABRAMOVICH rule. The proportion of the zero coefficients erroneously included
%%%                       in the model defined by the parameter q=0.2;

%--------------------------------------------------------------------------
% WINDOWS PARAMETERS
%--------------------------------------------------------------------------
sig_winds=[ones(size(h1))*0.6 ; ones(size(h2))*0.6];    % Gaussian parameter

window_type=1;  % window=1 for uniform, window=2 for Gaussian

TYPE=10;            % TYPE IS A SYMMETRY OF THE WINDOW
%                   % 00 SYMMETRIC
%                   % 10 NONSYMMETRIC ON X1 and SYMMETRIC ON X2
%                   % 11 NONSYMMETRIC ON X1,X2  (Quadrants)


%--------------------------------------------------------------------------
% MODELLING
%--------------------------------------------------------------------------
y=im2double(imread('image_Cameraman256.png'));
%y=im2double(imread('image_Lena512.png'));
%y=im2double(imread('image_Cheese128.png'));
%y=im2double(imread('image_Boats512.png'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Anisotropic multiresolution (MR) LPA Denoising Demo ')
disp('-----------------------------------------------------------------------------------')
disp(' ')

tic
%---------------------------------------------------------
% Images SIMULATION
%---------------------------------------------------------
[size_z_1,size_z_2]=size(y);
if addnoise==1
    init=0;%2055615866; 
    randn('seed', init);
    n=sigma_noise*randn(size_z_1,size_z_2);
    z = y + n;
else
    z = y;
end
sigma=function_stdEst2D(z);

%---------------------------------------------------------
% Kernels construction
%---------------------------------------------------------

% calling kernel creation function
[kernels, kernels_higher_order]=function_CreateLPAKernels(m,h1,h2,TYPE,window_type,ndir,sig_winds,1);
disp('kernels created');
toc
YICI_Final1=0; var_inv=0;
htbar=timebar('MR filtering running','Progress');
%%%%%%%%%%%%%%%%% Loops over directions   %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s1=1:ndir     % directional index

    for s2=1:lenh                % kernel size index
        gh=kernels_higher_order{s1,lenh-s2+1,1}(:,:,1);   % gets single kernel from the cell array (starting from the largest !!)
        if s2==1, clear GH
            [maxh11,maxh12]=size(gh); GH(1:maxh11,1:maxh11,1:lenh)=0;  end
        [N1gh, N2gh]=size(gh);
        LL= [(maxh11+1)/2-(N1gh-1)/2:(maxh11+1)/2+(N1gh-1)/2] ;
        GH(LL,LL,s2)=gh;
    end
    %%%% Filtering %%%%%%%%%%%%%%%%%%
    [Y_out,Var]=function_MR_Filtering(GH,thresholding_type,threshold_parameter);
    axes('position',[0.5 0.94 0.0001 0.0001]); axis off, title(['\theta=', num2str(360/ndir*(s1-1)),'^o']);
    %%%% Fusing of directional estimates
    YICI_Final1=YICI_Final1+Y_out./(Var+eps);
    var_inv=var_inv+1./(Var+eps);
    %%%%%%% Variance Map is the pointwise variance of the directional
    %%%%%%% estimates obtained after thresholding
    if s1==1
        figure
        figure_number=gcf;
        subplot(2,2,1), imshow(Var, []), title(['MR Variance Map (\theta=', num2str(360/ndir*(s1-1)),'^o)']), end
    if s1==round(1+1*ndir/4)
        if gcf~=figure_number,         figure(figure_number), end
        subplot(2,2,2), imshow(Var, []),  title(['MR Variance Map (\theta=', num2str(360/ndir*(s1-1)),'^o)']), end
    if s1==(1+2*ndir/4)
        if gcf~=figure_number,         figure(figure_number), end
        subplot(2,2,3), imshow(Var, []), title(['MR Variance Map (\theta=', num2str(360/ndir*(s1-1)),'^o)']), end
    if s1==(1+3*ndir/4)
        if gcf~=figure_number,         figure(figure_number), end
        subplot(2,2,4), imshow(Var, []),  title(['MR Variance Map (\theta=', num2str(360/ndir*(s1-1)),'^o)']), end
    timebar(htbar,s1/ndir);
end     %for s1 directions
toc
close(htbar);
y_hat=YICI_Final1./var_inv;
y_hat=max(0,min(1,y_hat));
figure
subplot(1,3,1), imshow(y), title('True image'),
subplot(1,3,2), imshow(z), title('Noisy image'),
subplot(1,3,3), imshow(y_hat), title('MR estimate'),
%%%%%%%%%%%%%   END OF ALGORITHM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[error_criteria,Error_labels]=function_Errors(y,y_hat,z,2);
