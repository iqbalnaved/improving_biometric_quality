% Anisotropic LPA-ICI Inverse Halftoning Demo (demo_InverseHalftoning)
%
% Alessandro Foi - Tampere University of Technology - 2003-2005
% -----------------------------------------------------------------------
%
% Reconstructs a continuous-tone image from a given error-diffusion halftone
% image. Inverse-halftoning is performed using the Anisotropic LPA-ICI
% deconvolution based on the RI (Regularized Inverse) and RWI (Regularized
% Wiener Inverse) estimates with ICI-driven adaptive scale selection.
%
% The two files in which the filtering is performed are needed:
% function_InverseHalftoning_RI  and  function_InverseHalftoning_RW
%
%
% Error-diffusion halftoning is modeled as a convolutional process of the
% form        z  =  p * y  +  q * n
% where  z  is the halftone,  p  and  q  are convolutional kernels (which
% depend on the error-diffusion kernel), and  n  is Gaussian white noise.
% Here  *  stands for the convolution operation.
%  (see  Kite et al., “Digital Image Halftoning as 2-D Delta-Sigma
%                          Modulation“, Proc. IEEE ICIP, 1997.       )
%
% The error-diffusion kernel is assumed as known.
%
% Other key variables:
% zRI      : "unfiltered" regularized inverse estimate
% y_hat_RI : filtered RI estimate
% zRW      : "unfiltered" regularized Wiener inverse estimate
% y_hat_RW : filtered RW estimate (final estimate)
%
%
% This code implements the algorithm and replicates the results presented in
% Foi, A., V. Katkovnik, K. Egiazarian, and J. Astola,
% “Inverse halftoning based on the anisotropic LPA-ICI deconvolution”,
% Proc. Int. TICSP Workshop Spectral Meth. Multirate Signal Process.,
% SMMSP 2004, Vienna, pp. 49-56, 2004.
%


clear all
close all

% ----------------------------------------------------------------------
% INPUT SIGNAL SETTING
%-----------------------------------------------------------------------
% Error-diffusion kernel
err_kernel=1;  % 1 is Floyd-Steinberg , 2 is Jarvis et al.
% Continuous-tone image
y=im2double(imread('image_Boats512.png'));
% y=im2double(imread('image_Lena512.png'));
% y=im2double(imread('image_Peppers512.png'));




% ----------------------------------------------------------------------
% Algorithm's main options
% ----------------------------------------------------------------------
estimate_sigma=0;       % estimate variance of the noise or use true value? (1 or 0)
do_wiener=1;            % do_wiener = 0 does not perform RW stage ( = 1 performs, dafault )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALGORITHM PARAMETERS (it is recommended not to modify the following parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------
% LPA WINDOWS PARAMETERS
%---------------------------------------------------------
% h1RI=[1 3 5 6 11];
% h1RW=[1 3 5 8 17];
h1RI=[2 3 4 5 6 8]; % set of scales used for LPA in RI
h2RI=[1 2 2 2 3 3];

%h1RW=[1 2 3 5  6 8];
%h2RW=[1 1 1 1  2 4 ];
h1RW=[1 2 3 5 4 4 6 5 8]; % set of scales used for LPA in RW
h2RW=[1 1 1 1 2 3 2 3 4];
% h2RI=ones(size(h1RI)); % To have line LPA kernels in RI
% h2RW=ones(size(h1RW)); % To have line LPA kernels in RW
% h2RW=max(1,ceil(h1RW*tan(0.5*pi/8))); %% sectorial kernels

alphaorderRI=-0;      % LPA order-mixture parameter (-1 zero order, 0 first order)
alphaorderRW=-1;      % LPA order-mixture parameter (-1 zero order, 0 first order)

% ICI threshold
% -------------
% GammaParameterRI  ICI Gamma for RI
% GammaParameterRW  ICI Gamma for RW
% (the larger the Gamma, the larger the adaptive kernel size and smoothing)
%
GammaParameterRI=[0.34];
GammaParameterRW=[0.33];

% Regularization Parameters:
% --------------------------
% Regularization_epsilon_RI  Regularization parameter of Regularized Inverse operator of RI stage
% Regularization_epsilon_RW  Regularization parameter of Regularized Wiener Inverse operator of RW stage
Regularization_epsilon_RI=2;
Regularization_epsilon_RW=34;

disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Anisotropic LPA-ICI Inverse Halftoning Demo ')
disp('-----------------------------------------------------------------------------------')
disp(' ')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  GENERATES BLURRED AND NOISY OBSERVATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yN,xN]=size(y);
size_z_1=yN; size_z_2=xN;

if err_kernel==1
    error_kernel=[0 0 7;3 5 1];  %% Floyd-Steinberg
    K = 2.03;   %% Gain constant
    disp(' ')
    disp(' Floyd-Steinberg error-diffusion')
    disp(' ')
elseif err_kernel==2
    error_kernel=[0 0 0 7 5;3 5 7 5 3;1 3 5 3 1];  %% Jarvis et al.
    K = 4.45;   %% Gain constant
    disp(' ')
    disp(' Jarvis et al. error-diffusion')
    disp(' ')
end
error_kernel=error_kernel/sum(error_kernel(:));  %% normalizes error-kernel
z=function_ErrorDiffusion(y,error_kernel);   %%% generates halftone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Kite's error-diffusion convolutional model  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Z = PY + QN
h = zeros(size_z_1,size_z_2);
h(1:size(error_kernel,1),1:size(error_kernel,2))=error_kernel;  %% zero-pads error-kernel
h=circshift(h,[0 -(size(error_kernel,2)-1)/2]);
H = fft2(h);
P=K./(1+(K-1)*H);
Q=(1-H)./(1+(K-1)*H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
imshow(z);
title('Error-Diffusion Halftone  z')
drawnow

version -release; % get matlab release
matlab_R=str2num(ans);

lenhRI=length(h2RI); % number of scales in RI
lenhRW=length(h1RW); % number of scales in RW
ndirRI=8; % number of directions
ndirRW=8; % number of directions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    FILTERING STARTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%  RI Estimation %%%%%%%%%%%%%%%%%%%%%
disp('  ')
disp('starting RI ...     ')
function_InverseHalftoning_RI ; % RI WORKING PROGRAM
RI_toc=toc;
disp(sprintf(repmat('\b',[1 2+20]))),
disp(['RI completed in ',num2str(RI_toc),' seconds.           '])
[Err_RI,Err_labels]=function_Errors(y,y_hat_RI,z);  %% computes error criteria
%%% PRINT RESULTS TO SCREEN
number_of_digits=7;
if do_wiener>=1
    disp('starting RW ...     ')
    Err_RW=repmat('wait ...',[size(Err_RI,1) 1]);
else
    Err_RW=[repmat(' ',[size(Err_RI,1) floor((number_of_digits-1)/2)]),repmat('--',[size(Err_RI,1) 1]),repmat(' ',[size(Err_RI,1) ceil((number_of_digits-1)/2)])];
end
disp(' ');
tab_title=[repmat(' ',[1 size(Err_labels,2)+1+floor(number_of_digits/2)]),'RI',repmat(' ',[1 number_of_digits+1]),'RW'];
tab_content=[Err_labels,repmat(' ',[size(Err_RI,1) 1]),num2str(Err_RI,number_of_digits),repmat(' ',[size(Err_RI,1) 2]),Err_RW];
disp(tab_title);
disp(tab_content);
figure,  imshow(y_hat_RI), title('LPA-ICI Regularized Inverse (RI) estimate'),
if do_wiener>=1    %%%%%%%%%%%%  RW Estimation %%%%%%%%%%%%%%%%%%%%%
    Wiener_Pilot=abs(fft2(y_hat_RI));   %%% WIENER PILOT ESTIMATE
    function_InverseHalftoning_RW % RW WORKING PROGRAM
    RW_toc=toc;
    Err_RW=function_Errors(y,y_hat_RW,z);  %% computes error criteria
    %%% PRINT RESULTS TO SCREEN
    disp(sprintf(repmat('\b',[1 numel(tab_content)+size(tab_content,1)+numel(tab_title)+5+20]))),
    disp(['RW completed in ',num2str(RW_toc-RI_toc),' seconds.   Total time: ',num2str(RW_toc),' seconds.'])
    disp(' ');
    disp(tab_title);
    disp([Err_labels,repmat(' ',[size(Err_RI,1) 1]),num2str(Err_RI,number_of_digits),repmat(' ',[size(Err_RI,1) 2]),num2str(Err_RW,number_of_digits)]);
    disp(' ')
    figure, imshow(y_hat_RW);title('LPA-ICI Regularized Wiener Inverse (RWI) estimate')
else
    disp(' ')
    disp('(!!!) RW skipped (!!!)                 ( do_wiener=0 )')
    disp(' ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%     FILTERING ENDS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SOME FIGURES WITH RANDOMLY SELECTED DETAILS OF IMAGES FROM THE ALGORITHM
if (do_wiener>=1)&(ndirRW==8)&(ndirRI==8)
    for aaa=1:3 % HOW MANY DIFFERENT FIGURES (DIFFERENT DETAILS) TO SHOW
        screensize = get(0,'screensize');       % User's screen size [1 1 width height]
        size_patch=screensize([4,3])/14;
        range1=(1+round((size_z_1-1-size_patch(1))*rand))+[1:size_patch(1)];
        range2=(1+round((size_z_2-1-size_patch(2))*rand))+[1:size_patch(2)];

        figure
        subplot(3,3,1)
        imshow(y(range1,range2)),
        if matlab_R<14
            title('original image  y')
        else
            title('original image  $y$','interpreter','latex')
        end

        subplot(3,3,2)
        imshow(z(range1,range2)),
        if matlab_R<14
            title('Error-Diffusion Halftone  z')
        else
            title('Error-Diffusion Halftone \ \  $z$','interpreter','latex')
        end

        MaxPQ=max(max(abs(P(:))),max(abs(Q(:))));
        subplot(3,6,5)
        mesh(fftshift(abs(P(1:ceil(size_z_1/30):end,1:ceil(size_z_2/30):end)))), axis square tight
        set(gca,'xtick',[],'xticklabel',[]);
        set(gca,'ytick',[],'yticklabel',[]);
        axis_get=axis;
        axis_get(5:6)=[0,MaxPQ];
        axis(axis_get);
        view([-25,10]);
        camproj('perspective'); box on
        axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.93 0.0001 0.0001]);  axis off,
        if matlab_R<14
            title('|P|')
        else
            title('$|P|$','interpreter','latex')
        end
        

        subplot(3,6,6)
        mesh(fftshift(abs(Q(1:ceil(size_z_1/30):end,1:ceil(size_z_2/30):end)))), axis square tight
        set(gca,'xtick',[],'xticklabel',[]);
        set(gca,'ytick',[],'yticklabel',[]);
        axis_get=axis;
        axis_get(5:6)=[0,MaxPQ];
        axis(axis_get);
        view([-25,10]);
        camproj('perspective'); box on
        axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.93 0.0001 0.0001]);  axis off,
        if matlab_R<14
            title('|Q|')
        else
            title('$|Q|$','interpreter','latex')
        end
        


        subplot(3,3,4)
        imshow(zRI(range1,range2)),
        if matlab_R<14
            title('Regularized Inverse  z^{RI}')
        else
            title('Regularized Inverse \ \  $z^{RI}$','interpreter','latex')
        end

        subplot(3,3,7)
        imshow(zRW(range1,range2)),
        if matlab_R<14
            title('Regularized Wiener Inverse  z^{RWI}')
        else
            title('Regularized Wiener Inverse \ \  $z^{RWI}$','interpreter','latex')
        end

        subplot(3,3,6)
        imshow(y_hat_RI(range1,range2)),
        if matlab_R<14
            title('Filtered Regularized Inverse estimate  y^{\^RI}')
        else
            title('Filtered Regularized Inverse estimate \ \  $\hat{y}^{RI}$','interpreter','latex')
        end

        subplot(3,3,9)
        imshow(y_hat_RW(range1,range2)),
        if matlab_R<14
            title('Filtered Regularized Wiener Inverse estimate  y^{\^RWI}')
        else
            title('Filtered Regularized Wiener Inverse estimate \ \  $\hat{y}^{RWI}$','interpreter','latex')
        end

        subplot(3,3,5)
        imshow(h_optRI(range1,range2),[]),
        if matlab_R<14
            title('Adaptive scales  h^{+}( \cdot ,7\pi/4)    (RI)')
        else
            title('Adaptive scales \ \  $h^{+}(\ \cdot \ ,7\pi/4)$ \ \  (RI)','interpreter','latex')
        end

        subplot(3,3,8)
        imshow(h_opt_Q(range1,range2,5),[]),
        if matlab_R<14
            title('Adaptive scales  h^{+}( \cdot ,\pi)    (RWI)')
        else
            title('Adaptive scales \ \  $h^{+}(\ \cdot \ ,\pi)$ \ \  (RWI)','interpreter','latex')
        end

    end %%% LOOP ON RANDOM DETAILS
end %%% FIGURES WITH DETAILS
%%% end of code

