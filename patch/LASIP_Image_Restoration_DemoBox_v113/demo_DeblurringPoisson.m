% Anisotropic LPA-ICI Poissonian Deconvolution Demo (demo_DeblurringPoissonian)
%
% Alessandro Foi - Tampere University of Technology - 2004-2006
% -----------------------------------------------------------------------
%
% Performs deblurring (deconvolution) from observations which are
% blurred and noisy. Noise is modeled as a Poisson process.
%
% The RI (Regularized Inverse) and RWI (Regularized
% Wiener Inverse) Deconvolution Algorithm with Anisotropic LPA-ICI
% adaptive estimate selection is used.
%
% The two files in which the filtering is performed are needed:
% function_DeblurringPoisson_RW  and  function_DeblurringPoisson_RI
%
%
%
% Observation model:
% z ~ P(lambda x y * v)
%
% z      : blurred noisy observation
% P      : Poissonian distribution
% y      : true image (assumed as unknown)
% lambda : scaling factor (assumed as known)  Larger value corresponds to higher BSNR.
% v      : point-spread function (assumed as known)
% *      : convolution
% x      : multiplication
%
%
% Other key variables:
% zRI      : "unfiltered" regularized inverse estimate
% y_hat_RI : filtered RI estimate
% zRW      : "unfiltered" regularized Wiener inverse estimate
% y_hat_RW : filtered RW estimate (final estimate)
%
%
% This code implements the algorithm and replicates the results presented in
% Foi, A., S. Alenius, M. Trimeche, V. Katkovnik, and K. Egiazarian,
% “A spatially adaptive Poissonian image deblurring”,
% Proc. of IEEE 2005 Int. Conf. Image Processing, ICIP 2005, September 2005.
%

clear all
close all


lambda=17600; %%% multiplicator for the image. Larger value corresponds higher BSNR.  lambda=17600  produce a randomness of the noise comparable to real pictures taken with a cameraphone.


% ICI & Reg Parameters Initiation
GammaParameterRI=[1.35];             % GammaParameterRI  ICI Gamma for RI
Regularization_epsilon_RI=0.027;
GammaParameterRW=[1.35];             % GammaParameterRW  ICI Gamma for RW
Regularization_epsilon_RW=0.077;


alphaorderRI=-0.8;                  %coefficients for order-mixture  -1 is zero order, 0 is first order (positive coefficients sharpen the image (not reccomended))
alphaorderRW=-1.0;                  %coefficients for order-mixture  -1 is zero order, 0 is first order (positive coefficients sharpen the image (not reccomended))


Precooked=1;          %% Use optimized "precooked" kernels for RI stage (default)
%                     %% [ number of directions (ndirRI) must be either 8 or 4 ]

do_wiener=1;          %% do_wiener = 0 does not perform RW stage ( = 1 performs, dafault )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------
% WINDOWS PARAMETERS
% h1--Corresponds TO X1
% h2--Corresponds to X2

h1RI=[1 3 5 6 11];                    % set of scales used for LPA in RI
h1RW=[1 3 5 8 13];                    % set of scales used for LPA in RW

h2RI=ones(size(h1RI));                %% To have line kernels for RI
h2RW=max(1,ceil(h1RW*tan(0.5*pi/8))); %% sectorial kernels
%h2RW=ones(size(h1RW));               %% To have line kernels also for RW

ndirRI=8; % number of directions
ndirRW=8; % number of directions

tic


disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Anisotropic LPA-ICI Poissonian Deconvolution Demo ')
disp('-----------------------------------------------------------------------------------')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD TRUE IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=im2double(imread('image_Cameraman256.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  GENERATES BLURRED AND NOISY OBSERVATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yN,xN]=size(y);
y=y*lambda;            %%% rescales image
size_z_1=yN; size_z_2=xN;
v=ones(9); v=v./sum(v(:));   % v is the blur Point-Spread Function (PSF)
[ghy,ghx]=size(v);
%%  BLURRING  %%%%%%%%%
big_v=zeros(yN,xN); big_v(1:ghy,1:ghx)=v; big_v=circshift(big_v,-round([(ghy-1)/2 (ghx-1)/2])); % pad PSF with zeros to whole image domain, and centers it.
V=fft2(big_v); % Frequency responde of the PSF
y_blur=real(ifft2(V.*fft2(y))); % performs blurring (convolution is obtained by product in frequency domain)
%%%  GENERATE NOISY BLURRED SIGNAL  %%%%%%%%%%%%
init=2043152662; randn('seed', init); rand('seed', init);  %%% FIX SEED FOR RANDOM PROCESSES  (OPTIONAL)
z=poissrnd(y_blur.*(y_blur>0));   %%  Poissonian process with mean y_blur

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Algorithm's Structural Parameters (DO NOT TOUCH THESE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lenhRI=length(h2RI); % number of scales in RI
lenhRW=length(h1RW); % number of scales in RW

version -release; % get matlab release
matlab_R=str2num(ans);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURE OF OBSERVATION SHOWING SIGNAL DEPENDANT NOISE CHARACTERISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state',sum(100*clock));
line1=(1+round((size_z_1-1)*rand));
line2=(1+round((size_z_1-1)*rand));
line3=(1+round((size_z_1-1)*rand));
line4=(1+round((size_z_1-1)*rand));
figure
subplot (1,3,1)
imshow(z/lambda)
line([1 size(z,2)],[line1 line1],'Color',[1 0 0])
line([1 size(z,2)],[line2 line2],'Color',[0 1 0])
line([1 size(z,2)],[line3 line3],'Color',[0 0 1])
line([1 size(z,2)],[line4 line4],'Color',[0.8 0.8 0])
title('blurred and noisy observation  z')
subplot (2,3,2)
plot(z(line1,:),'Color',[1 0 0]), axis([[1 size(z,2)],[0 lambda]])
title('cross-section')
subplot (2,3,3)
plot(z(line2,:),'Color',[0 0.9 0]),axis([[1 size(z,2)],[0 lambda]])
title('cross-section')
subplot (2,3,5)
plot(z(line3,:),'Color',[0 0 1]),axis([[1 size(z,2)],[0 lambda]])
title('cross-section')
subplot (2,3,6)
plot(z(line4,:),'Color',[0.7 0.7 0]),axis([[1 size(z,2)],[0 lambda]])
title('cross-section')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    FILTERING STARTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%  RI Estimation %%%%%%%%%%%%%%%%%%%%%
disp('  ')
disp('starting RI ...     ')
function_DeblurringPoisson_RI ; % RI WORKING PROGRAM
RI_toc=toc;
disp(sprintf(repmat('\b',[1 2+20]))),
disp(['RI completed in ',num2str(RI_toc),' seconds.           '])
[Err_RI,Err_labels]=function_Errors(y/lambda,y_hat_RI/lambda,z/lambda);  %% computes error criteria
%%% PRINT RESULTS TO SCREEN
number_of_digits=5;
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
figure,  imshow(y_hat_RI/lambda), title('LPA-ICI Regularized Inverse (RI) estimate'),
if do_wiener>=1    %%%%%%%%%%%%  RW Estimation %%%%%%%%%%%%%%%%%%%%%
    Wiener_Pilot=abs(fft2(y_hat_RI));   %%% WIENER PILOT ESTIMATE
    function_DeblurringPoisson_RW % RW WORKING PROGRAM
    RW_toc=toc;
    Err_RW=function_Errors(y/lambda,y_hat_RW/lambda,z/lambda);  %% computes error criteria
    %%% PRINT RESULTS TO SCREEN
    disp(sprintf(repmat('\b',[1 numel(tab_content)+size(tab_content,1)+numel(tab_title)+5+20]))),
    disp(['RW completed in ',num2str(RW_toc-RI_toc),' seconds.   Total time: ',num2str(RW_toc),' seconds.'])
    disp(' ');
    disp(tab_title);
    disp([Err_labels,repmat(' ',[size(Err_RI,1) 1]),num2str(Err_RI,number_of_digits),repmat(' ',[size(Err_RI,1) 2]),num2str(Err_RW,number_of_digits)]);
    disp(' ')
    y_hat_RW=y_hat_RW/lambda;
    zRW=zRW/lambda;
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
y_hat_RI=y_hat_RI/lambda;
zRI=zRI/lambda;
z=z/lambda;
y=y/lambda;

%%% SOME FIGURES WITH RANDOMLY SELECTED DETAILS OF IMAGES FROM THE ALGORITHM
if (do_wiener>=1)&(ndirRW==8)&(ndirRI==8)
    for aaa=1:2 % HOW MANY DIFFERENT FIGURES (DIFFERENT DETAILS) TO SHOW
        screensize = get(0,'screensize');       % User's screen size [1 1 width height]
        size_patch=screensize([4,3])/14;
        rand('state',sum(100*clock));
        range1=(1+round((size_z_1-1-size_patch(1))*rand))+[1:size_patch(1)];
        range2=(1+round((size_z_2-1-size_patch(2))*rand))+[1:size_patch(2)];
        drawnow
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
        imshow(h_optRWQ(range1,range2,5),[]),
        if matlab_R<14
            title('Adaptive scales  h^{+}( \cdot ,\pi)    (RWI)')
        else
            title('Adaptive scales \ \  $h^{+}(\ \cdot \ ,\pi)$ \ \  (RWI)','interpreter','latex')
        end
        axes('position',[0.5 0.94 0.0001 0.0001]); axis off, title('Detailed view of a fragment of the image')
    end %%% LOOP ON RANDOM DETAILS
end %%% FIGURES WITH DETAILS
%%% end of code

