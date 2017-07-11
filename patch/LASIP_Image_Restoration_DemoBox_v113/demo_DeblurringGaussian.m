% Anisotropic LPA-ICI Deconvolution Demo (demo_DeblurringGaussian)
%
% Alessandro Foi - Tampere University of Technology - 2003-2006
% -----------------------------------------------------------------------
%
% Performs deblurring (deconvolution) from observations which are
% blurred and noisy. The RI (Regularized Inverse) and RWI (Regularized
% Wiener Inverse) Deconvolution Algorithm with Anisotropic LPA-ICI
% adaptive estimate selection is used.
%
% The two files in which the filtering is performed are needed:
% function_DeblurringGaussian_RW  and  function_DeblurringGaussian_RI
%
%
%
% Observation model:
% z=y*v+n
% z : blurred noisy observation
% y : true image (assumed as unknown)
% v : point-spread function (assumed as known)
% n : gaussian white noise
% * : convolution
%
% Other key variables:
% zRI      : "unfiltered" regularized inverse estimate
% y_hat_RI : filtered RI estimate
% zRW      : "unfiltered" regularized Wiener inverse estimate
% y_hat_RW : filtered RW estimate (final estimate)
%
%
% This code implements the algorithm and replicates the results presented in
% Katkovnik, V., A. Foi, K. Egiazarian, and J. Astola,
% “Directional varying scale approximations for anisotropic signal processing”,
% Proc. of XII European Signal Process. Conf., EUSIPCO 2004, pp. 101-104, 2004.
%

clear all
close all

% ----------------------------------------------------------------------
% INPUT SIGNAL SETTING
%-----------------------------------------------------------------------
% Four experiments are prepared as described in the aforementioned paper.
%
% Each experiment has either different Point-Spread Function (PSF) or noise
% variance.
%
ExpNumber=1;   % set experiment number ( 1, 2, 3 or 4 )

% ----------------------------------------------------------------------
% Algorithm's main options
% ----------------------------------------------------------------------
estimate_sigma=1;       % estimate variance of the noise or use true value? (1 or 0)
do_wiener=1;            % do_wiener = 0 does not perform RW stage ( = 1 performs, dafault )
estimate_derivative=1;  % estimates directional derivates at RW stage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALGORITHM PARAMETERS (it is recommended not to modify the following parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------
% LPA WINDOWS PARAMETERS
%---------------------------------------------------------
h1RI=[1 3 5 6 11]; % set of scales used for LPA in RI
h1RW=[1 3 5 8 17]; % set of scales used for LPA in RW
h2RI=ones(size(h1RI)); % To have line LPA kernels in RI
%h2RW=ones(size(h1RW)); % To have line LPA kernels in RW
h2RW=max(1,ceil(h1RW*tan(0.5*pi/8))); %% sectorial kernels
lenhRI=length(h2RI); % number of scales in RI
lenhRW=length(h1RW); % number of scales in RW
ndirRI=8; % number of directions
ndirRW=8; % number of directions
% NOTE: some parts of this demo make explicit use of 8 directions

% ICI threshold
% -------------
% GammaParameterRI  ICI Gamma for RI  (default 1.65)
% GammaParameterRW  ICI Gamma for RW  (default 1.35)
% (the larger the Gamma, the larger the adaptive kernel size and smoothing)
%
GammaParameterRI=1.35;
GammaParameterRW=1.25;

alphaorderRI=-1;      % LPA order-mixture parameter (-1 zero order, 0 first order)
alphaorderRW=-1;    % LPA order-mixture parameter (-1 zero order, 0 first order)

Precooked=1;          %% Use optimized "precooked" kernels for RI stage (default)
%% [ number of directions (ndirRI) must be either 8 or 4 ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  EXPERIMENT SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Observation Parameters
% ----------------------
% BSNR is the desired Blurred-Signal-to-Noise-Ratio
% v is the blur Point-Spread Function (PSF)
%
% Regularization Parameters:
% --------------------------
% Regularization_epsilon_RI  Regularization parameter of Regularized Inverse operator of RI stage
% Regularization_epsilon_RW  Regularization parameter of Regularized Wiener Inverse operator of RW stage
%

if ExpNumber==1 % Experiment 1
    BSNR=40; org_sigma=-1;             % noise level;  org_sigma=-1  forces the algorithm to use specified BSNR value (BSNR=40)
    v=ones(9); v=v./sum(v(:));
    Regularization_epsilon_RI=0.014;
    Regularization_epsilon_RW=0.11;
    y=im2double(imread('image_Cameraman256.png'));
end
if ExpNumber==2 % Experiment 2
org_sigma=sqrt(2)/255;            % noise level   sigma^2=2
    s1=0; for a1=-7:7; s1=s1+1; s2=0; for a2=-7:7; s2=s2+1; v(s1,s2)=1/(a1^2+a2^2+1); end, end;     v=v./sum(v(:));
    Regularization_epsilon_RI=0.045;
    Regularization_epsilon_RW=0.14;
    y=im2double(imread('image_Cameraman256.png'));
end
if ExpNumber==3 % Experiment 3
    org_sigma=sqrt(8)/255;             % noise level  sigma^2=8
    s1=0; for a1=-7:7; s1=s1+1; s2=0; for a2=-7:7; s2=s2+1; v(s1,s2)=1/(a1^2+a2^2+1); end, end;     v=v./sum(v(:));
    Regularization_epsilon_RI=0.07;
    Regularization_epsilon_RW=0.10;
    y=im2double(imread('image_Cameraman256.png'));
end
if ExpNumber==4 % Experiment 4
    org_sigma=7/255;                   % noise level  sigma^2=49
    v1=[1 4 6 4 1]/16; v=v1'*v1; v=v./sum(v(:));
    Regularization_epsilon_RI=0.27;
    Regularization_epsilon_RW=0.22;
    y=im2double(imread('IMAGE_lena512.png'));
end
if min(abs(ExpNumber-[1 2 3 4]))>0
    disp(' ');disp('ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'),
    disp('Experiment Number (ExpNumber) has to be 1, 2, 3 or 4'),disp(' ');
    break
end
disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Anisotropic LPA-ICI Deconvolution Demo ')
disp('-----------------------------------------------------------------------------------')
disp(' ')
disp(['Performing Experiment #',num2str(ExpNumber),' ...'])
if estimate_derivative==1&do_wiener==1
    disp('  Also directional derivatives will be estimated during RW stage.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  GENERATES BLURRED AND NOISY OBSERVATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yN,xN]=size(y);
size_z_1=yN; size_z_2=xN;
[ghy,ghx]=size(v);
%%  BLURRING  %%%%%%%%%
big_v=zeros(yN,xN); big_v(1:ghy,1:ghx)=v; big_v=circshift(big_v,-round([(ghy-1)/2 (ghx-1)/2])); % pad PSF with zeros to whole image domain, and centers it.
V=fft2(big_v); % Frequency responde of the PSF
y_blur=real(ifft2(V.*fft2(y))); % performs blurring (convolution is obtained by product in frequency domain)
%%%  ADDING NOISE to BLURRED SIGNAL  %%%%%%%%%%%%
init=0;
randn('seed',init);  %%% FIX SEED FOR RANDOM PROCESSES  (OPTIONAL)
if org_sigma==-1;   %% use BSNR in order to define value of sigma
    org_sigma=sqrt(norm(y_blur(:)-mean(y_blur(:)),2)^2 /(size_z_2*size_z_1*10^(BSNR/10))); % sigma of noise in blurred image given the desired BSNR
end
n = org_sigma*randn(size(y_blur)); % white Gaussian noise with variance org_sigma^2
z=y_blur+n; % observation
clear y_blur n big_v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  ESTIMATION OF THE SIGMA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if estimate_sigma==0
    sigma=org_sigma;
else
    sigma_est=function_stdEst2D(z,2);
    sigma=sigma_est;
    disp(['Estimated noise sigma = ',num2str(sigma_est),'   (true = ',num2str(org_sigma),',  misestimation = ',num2str(100*(sigma_est-org_sigma)/org_sigma),'%)']);
end
figure
imshow(z),
title('blurred and noisy observation  z')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Algorithm's Structural Parameters (DO NOT TOUCH THESE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lenhRI=length(h2RI); % number of scales in RI
lenhRW=length(h1RW); % number of scales in RW

version -release; % get matlab release
matlab_R=str2num(ans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    FILTERING STARTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%  RI Estimation %%%%%%%%%%%%%%%%%%%%%
disp('  ')
disp('starting RI ...     ')
function_DeblurringGaussian_RI ; % RI WORKING PROGRAM
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
    function_DeblurringGaussian_RW % RW WORKING PROGRAM
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
        drawnow
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
            title('blurred and noisy observation  z')
        else
            title('blurred and noisy observation \ \  $z$','interpreter','latex')
        end

        subplot(3,3,3)
        med_v=zeros(2*ceil(max(size(v))/3)+max(size(v))); med_v(ceil(max(size(v))/3)+1+floor((max(size(v))-size(v,1))/2):ceil(max(size(v))/3)+floor((max(size(v))-size(v,1))/2)+size(v,1),ceil(max(size(v))/3)+1+floor((max(size(v))-size(v,2))/2):ceil(max(size(v))/3)+floor((max(size(v))-size(v,2))/2)+size(v,2))=v;
        bar3(med_v,1,'w'),

        set(gca,'xtick',[size(med_v,1)/2],'xticklabel',[]);
        set(gca,'ytick',[size(med_v,1)/2],'yticklabel',[]);
        set(gca,'ztick',[0 max(v(:))],'zticklabel',[]);
        axis square tight
        box on
        camproj('perspective');
        axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.93 0.0001 0.0001]);  axis off,
        if matlab_R<14
            title('Point-Spread Function (PSF)  v')
        else
            title('Point-Spread Function (\textit{PSF}) \ \  $v$','interpreter','latex')
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
        axes('position',[0.5 0.94 0.0001 0.0001]); axis off, title('Detailed view of a fragment of the image')
        if estimate_derivative==1   %%% FIGURES WITH DERIVATIVE DETAILS
            sum_of_abs_derivative=sum(abs(yder),3); %%% EDGE MAP
            drawnow
            figure
            subplot(4,4,1)
            imshow(yd_hat_Q(range1,range2,4),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{+^{3\pi/4}}')
            else
                title('$\hat{\partial}_{_{+}3\pi/4}$','interpreter','latex')
            end

            subplot(4,4,2)
            imshow(yd_hat_Q(range1,range2,3),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{+^{\pi/2}}')
            else
                title('$\hat{\partial}_{_{+}\pi/2}$','interpreter','latex')
            end

            subplot(4,4,3)
            imshow(yd_hat_Q(range1,range2,2),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{+^{\pi/4}}')
            else
                title('$\hat{\partial}_{_{+}\pi/4}$','interpreter','latex')
            end

            subplot(4,4,5)
            imshow(yd_hat_Q(range1,range2,5),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{+^{\pi}}')
            else
                title('$\hat{\partial}_{_{+}\pi}$','interpreter','latex')
            end

            subplot(4,4,7)
            imshow(yd_hat_Q(range1,range2,1),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{+^{0}}')
            else
                title('$\hat{\partial}_{_{+}0}$','interpreter','latex')
            end

            subplot(4,4,9)
            imshow(yd_hat_Q(range1,range2,6),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{+^{5\pi/4}}')
            else
                title('$\hat{\partial}_{_{+}5\pi/4}$','interpreter','latex')
            end

            subplot(4,4,10)
            imshow(yd_hat_Q(range1,range2,7),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{+^{3\pi/2}}')
            else
                title('$\hat{\partial}_{_{+}3\pi/2}$','interpreter','latex')
            end

            subplot(4,4,11)
            imshow(yd_hat_Q(range1,range2,8),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{+^{7\pi/4}}')
            else
                title('$\hat{\partial}_{_{+}7\pi/4}$','interpreter','latex')
            end

            subplot(4,4,13)
            imshow(sum_of_abs_derivative(range1,range2),[0,1.3*max1])
            if matlab_R<14
                title('\Sigma_i|\partial_{\theta_i}|   (edge map)')
            else
                title('$\Sigma_i|\hat{\partial}_{\theta_i}|$ \ \ (edge map)','interpreter','latex')
            end

            subplot(4,4,6)
            imshow(z(range1,range2))
            if matlab_R<14
                title('blurred and noisy observation  z')
            else
                title('blurred and noisy observation \ \  $z$','interpreter','latex')
            end

            subplot(4,4,4)
            imshow(yder(range1,range2,2),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{\pi/4}')
            else
                title('$\hat{\partial}_{\pi/4}$','interpreter','latex')
            end

            subplot(4,4,8)
            imshow(yder(range1,range2,1),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{0}')
            else
                title('$\hat{\partial}_{0}$','interpreter','latex')
            end

            subplot(4,4,12)
            imshow(-yder(range1,range2,4),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{7\pi/4}')
            else
                title('$\hat{\partial}_{7\pi/4}$','interpreter','latex')
            end

            subplot(4,4,16)
            imshow(-yder(range1,range2,3),[min1,max1]*0.75)
            if matlab_R<14
                title('\partial_{3\pi/2}')
            else
                title('$\hat{\partial}_{3\pi/2}$','interpreter','latex')
            end

            subplot(4,4,14)
            imshow(y(range1,range2)),
            if matlab_R<14
                title('original image  y')
            else
                title('original image  $y$','interpreter','latex')
            end

            subplot(4,4,15)
            imshow(y_hat_RW(range1,range2)),
            if matlab_R<14
                title('Filtered Regularized Wiener Inverse estimate  y^{\^RWI}')
            else
                title('Filtered Regularized Wiener Inverse estimate \ \  $\hat{y}^{RWI}$','interpreter','latex')
            end
            axes('position',[0.5 0.94 0.0001 0.0001]); axis off, title('Detailed view of a fragment of the image')
        end %%% DERIVATIVE FIGURE
    end %%% LOOP ON RANDOM DETAILS
end %%% FIGURES WITH DETAILS
%%% end of code

