% Pointwise Shape-Adaptive DCT Color image denoising (using low-complexity SA-DCT)
%
% Alessandro Foi - Tampere University of Technology - 2005 - 2006   Public release v1.40 (October 2006)   C-code implementation by Kostadin Dabov
% -----------------------------------------------------------------------------------------------------------------------------------------------
%
%
%  IMPORTANT VARIABLES:
%
%  yRGB         :  original image (if available) in RGB space
%  zRGB         :  noisy observation in RGB space   (assumes AWGN model zRGB = yRGB + noise, where noise is iid gaussian)
%  yRGB_hat     :  hard-thresholding estimate in RGB space (obtained from yLumChrom_hat)
%  yRGB_hat_wi  :  Wiener-filter estimate in RGB space (final estimate, obtained from yLumChrom_hat_wi)
%  yLumChrom         :  original image (if available) in luminance-chrominance space
%  zLumChrom         :  noisy observation in luminance-chrominance space
%  yLumChrom_hat     :  hard-thresholding estimate in luminance-chrominance space
%  yLumChrom_hat_wi  :  Wiener-filter estimate in luminance-chrominance space (obtained using yLumChrom_hat as reference estimate in Wiener-filtering)
%  h_opt_Q   :  array with the adaptive-scale indexes for the kernels used in the construction of the adaptive-shape transform support
%
%
% The code implements the algorithm and reproduces the results published in:
% Foi, A., V. Katkovnik, and K. Egiazarian, “Pointwise Shape-Adaptive DCT for High-Quality Denoising and Deblocking of Grayscale and Color Images”,
% (accepted) IEEE Trans. Image Process., 2006.
% Foi, A., V. Katkovnik, and K. Egiazarian, “Pointwise Shape-Adaptive DCT Denoising with Structure Preservation in Luminance-Chrominance Space”,
% Proc. of the 2nd Int. Workshop on Video Process. and Quality Metrics for Consumer Electronics, VPQM2006, Scottsdale, January 2006.
%

clear trans_matrix Trian TrianRot_Int8 TrianW TrianRot_Int8W gh sigma_color yLumChrom_hat yLumChrom_hat_wi error_criteria;
close all; drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Options for this demonstration software    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_from_file=0;     %% loads original noise-free bitmap image from user-selected file (recommended =1)
enable_crop=1;        %% enables cropping of the initial image before beginning of algorithm
sigma_noise=25/255;   %% standard-deviation of the added noise (note: image range is assumed [0 1], hence the division by 255)
estimate_sigma=0;     %% estimate noise std from image?        (in this demo leave =0)
do_wiener=1;          %% enables Wiener-filter stage  (recommended =1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Filtering performance vs. complexity  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
speedup_factor=1;   %% any number between 1 and 5. speedup_factor>1 enables some speed-ups in the algorithm
%                   %% algorithm will run about speedup_factor times faster at the expense of image quality
%                   %% (recommended =1)        (note: number can be fractional, e.g. speedup_factor=4.4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Visualization options  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FiguresLumChromObs=0;    %% enables display of individual luminance and chrominance channel observations  (recommended =0)
FiguresLumChromEst=1;    %% enables display of individual luminance and chrominance channel estimates     (recommended =1)
ReUseFiguresLumChrom=1;  %% replaces figures of estimates during channel-by-channel processing (recommended =1)
figures_y_hats=1;        %% enables display of original, noisy, and SA-DCT estimate image     (recommended =1)
compute_errors=1;        %% enables calculation of error criteria (PSNR, MSE, etc.)  (recommended =1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Alternative filenames which can be used directly if load_from_file==0  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 image_filename='image_House256rgb.png';
% image_filename='image_Peppers512rgb.png';
% image_filename='image_Lena512rgb.png';
% image_filename='image_Baboon512rgb.png';
% image_filename='image_F16_512rgb.png';
% image_filename='image_Lake512rgb.png';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Advanced options to replicate most basic SA-DCT hardware implementations  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blocksize=0;   %% Constrain block-size  (recommended =0)
%              %% enforce  blocksize x blocksize  maximum block-size   (e.g. blocksize=16 or blocksize=8)
%              %% if blocksize<1 the maximum block-size is  2*max(h1)-1 x 2*max(h1)-1
%              %% ( h1 is defined below.  NOTE: h1 is not h1full, but a subset of h1full)
coef_align=1;  %% enable or disable coefficient alignment (recommended =1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  END OF DEMO OPTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Denoising algorithm's parameters  %%%%    SHOULD NOT BE MODIFIED   !!!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1full =  [1 2 3 5 7 9 12];   %% complete set of scales to be used for hard-thresholding      (recommended =[1 2 3 5 7 9 12])
h1Wfull = [1 2 3 5 6 8];      %% complete set of scales to be used for Wiener-filter          (recommended =[1 2 3 5 6 8])
sharparams=[-0.75 -0.70 -0.85 -0.90 -0.97 -1.00 -1.00];  %% order-mixture parameters (define the balance between zero- and first- order polynomial fitting for each scale)  (recommended =[-0.75 -0.70 -0.85 -0.9 -0.97 -1 -1])
DCTthrCOEF=0.77;    %% threshold-parameter  (recommended =0.77)
colormode='opp';            %% luminance-chrominance space to use. 'dct' or 'opp' (opponent) or 'pca'  (recommended ='opp')
chromDCTthrfactor=[1 1 1];  %% factors for threshold-parameter for the three channels in luminance chrominance  (recommended =[1 1 1])
max_overlap=[70 70 70];         %% limits overcompleteness (speed-up)  (recommended =70)  (note: lowering the value leads to a considerable speed-up in Pointwise SA-DCT hard-thresholding)
max_overlapW=[70 70 70];        %% limits overcompleteness (speed-up)  (recommended =70)  (note: lowering the value leads to a considerable speed-up in Pointwise SA-DCT Wiener-filtering)
if speedup_factor~=1   %%%  Speed-ups for "filtering performance vs. complexity" trade-off
    h1full=h1full(round(1:5/round((6-speedup_factor)):6));              %% complete set of scales to be used for hard-thresholding
    h1Wfull=h1Wfull(round(1:4/round((5-(4/5)*speedup_factor)):5));      %% complete set of scales to be used for Wiener-filter
    max_overlap=[round((6-speedup_factor).^2.64), ceil(0.25*(6-speedup_factor).^3.5), ceil(0.25*(6-speedup_factor).^3.5)];         %% limits overcompleteness (speed-up)  (note: lowering the value leads to a considerable speed-up in Pointwise SA-DCT hard-thresholding)
    max_overlapW=[ceil(0.25*(6-speedup_factor).^3.5), ceil(0.081*(6-speedup_factor).^4.2), ceil(0.081*(6-speedup_factor).^4.2)];   %% limits overcompleteness (speed-up)  (note: lowering the value leads to a considerable speed-up in Pointwise SA-DCT Wiener-filtering)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Pointwise Shape-Adaptive DCT Color Denoising   -   Public demo release')
disp('-----------------------------------------------------------------------------------')
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------
% Choose image from file (gets filename)
%---------------------------------------------------------
if load_from_file==1
    [bitmap_name,bitmap_path] = uigetfile('*.*','Select noise-free original image');
    if isequal(bitmap_name,0) || isequal(bitmap_path,0)
        disp(' ');
        disp('User pressed cancel');
        return
    else
        disp(['User selected   ', fullfile(bitmap_path,bitmap_name)]);
        disp(' ');
        image_filename=[bitmap_path,bitmap_name];
    end
else
    if exist(image_filename,'file')
        disp(['Loading image   ', image_filename]);
        disp(' ');
    else
        disp(['  !!!   Specified file "',image_filename ,'" does not exist   !!!'])
        disp(' ')
        disp(' ')
        return
    end
end
%---------------------------------------------------------
% load image
%---------------------------------------------------------
yRGB=im2double(imread(image_filename));

%---------------------------------------------------------
% Checks if color or grayscale image
%---------------------------------------------------------
if size(yRGB,3)==1  %% checks if image is grayscale
    disp(' ')
    disp(['  !!!   This demo is meant for color images.  Please use  demo_SADCT_denoising.m  for grayscale images   !!!'])
    disp(' ')
    disp(' ')
    yRGB=cat(3,yRGB,yRGB,yRGB);   %% builds RGB image out of grayscale one.
end

%---------------------------------------------------------
% Cropping & original image display
%---------------------------------------------------------
if enable_crop
    display('---- PLEASE CROP REGION FROM IMAGE  ----')
    figure
    imshow(yRGB);
    title('---- PLEASE CROP REGION FROM IMAGE  ----')
    hold on
    yRGB=imcrop(yRGB);
    hold off
    disp(sprintf([repmat('\b',[1,41]),' ']));
    close
end
if figures_y_hats
    figure
    imshow(yRGB), title('Original image')
    drawnow
end

%---------------------------------------------------------
% Addition of noise  (generates noisy observation)
%---------------------------------------------------------
init=0;
randn('seed', init);
zRGB = yRGB + sigma_noise*randn(size(yRGB));        %%% z = y + noise

[size_z_1,size_z_2,size_z_3]=size(zRGB);
if min(size_z_1,size_z_2)<2
    disp(' '); disp(['  !!!   Image is a vector,  exiting...  !!!']); disp(' '); disp(' ')
    return
end

%---------------------------------------------------------
% Display noisy observations and calculate errors
%---------------------------------------------------------
if compute_errors
    disp('Initial (noisy) observation criteria values:')
    function_Errors(yRGB,zRGB,2);
    disp(' ');
end
if figures_y_hats
    figure
    imshow(zRGB), title('Noisy observation')
    drawnow
end

if strcmp(colormode,'pca')
    [zLumChrom colormode l2normLumChrom]=function_rgb2LumChrom(zRGB,colormode);
else
    [zLumChrom l2normLumChrom]=function_rgb2LumChrom(zRGB,colormode);
end
yLumChrom=function_rgb2LumChrom(yRGB,colormode);

%%% STANDARD-DEVIATION %%%
if estimate_sigma==1
    for colorchannel=[1 2 3];
        sigma_color(colorchannel)=function_stdEst2D(zLumChrom(:,:,colorchannel))/l2normLumChrom(colorchannel);   %%% estimates standard deviation from image (assumes perfect AWGN model)
    end
    sigma=min(sigma_color)*l2normLumChrom(1);
else
    sigma=sigma_noise*l2normLumChrom(1);
end
%%% sigma is the standard-deviation for the luminance channel
%%% the standard-deviations for the chrominance channels is expressed as sigma*l2normLumChrom(colorchannel)/l2normLumChrom(1)

gammaICI=max(0.8,2.4./(log(1+60*sigma)));    %%% the maximum is used in order to ensure large enough gamma for strong levels of noise (above sigma=75)


%---------------------------------------------------------
% Definition of the set of scales h1
%---------------------------------------------------------
if (sigma>(40/255))||(speedup_factor~=1)     %%% if noise std is large use one extra scale (slows down the algorithm)
    h1=h1full;            %%% set of scales to be used in the hard-thresholding algorithm
    h1W=h1Wfull;          %%% set of scales to be used in the Wiener-filter algorithm
else                        %%% if noise std is not too large, fewer scales can be used (speed-up)
    h1=h1full(find(h1full<10));    %%% set of scales to be used in the hard-thresholding algorithm
    h1W=h1Wfull(find(h1Wfull<7));   %%% set of scales to be used in the Wiener-filter algorithm
end
if blocksize>0   %% limits the maximum scale so not to exceed the chosen maximum block-size
    h1boundU=ceil((blocksize+1)/2);
    h1boundL=floor((blocksize+1)/2);
    h1(find(h1>=h1boundU))=h1boundU;
    h1=h1(1:min(find(h1==max(h1))));
    h1W(find(h1W>=h1boundU))=h1boundU;
    h1W=h1W(1:min(find(h1W==max(h1W))));
end

directional_resolution=8;   %%% DO NOT CHANGE THE NUMBERS OF DIRECTIONS!!  = 8
lenh=numel(h1);
lenhW=numel(h1W);
h_max=h1(lenh);
h_maxW=h1W(lenhW);
SSS=[1:(2*h_max-1)^2];              %%% vector with the possible sizes of the adaptive-shapes
T=DCTthrCOEF*sqrt(2*log(SSS)+1);    %%% Universal-threshold (depends on the size of the adaptive-shape)

%---------------------------------------------------------
% Display information about algorithm's key parameters
%---------------------------------------------------------
disp(['noise sigma in RGB = ',num2str(sigma*255/l2normLumChrom(1)),'   noise sigma in luminance = ',num2str(sigma*255)]);
disp(['GammaICI = ',num2str(gammaICI),'   DCTthrCOEF = ',num2str(DCTthrCOEF),'    h_max = ',num2str(h_max),'   h_maxW = ',num2str(h_maxW)]);
if speedup_factor~=1
    disp(['   !!  Speed-up factor = ',num2str(speedup_factor), '  !!'])
end
if blocksize>0
    disp(['   !!  User-imposed constrained maximum block-size  ',num2str(blocksize),'x',num2str(blocksize), '  !!'])
end
disp(' ');
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LPA-ICI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------
% Kernels construction
%---------------------------------------------------------
% calling kernel creation function
[kernels, kernels_higher_order]=function_CreateLPAKernels([0 0],h1,ones(size(h1)),10,1,1,ones(2,lenh),1);
[kernelsb, kernels_higher_orderb]=function_CreateLPAKernels([1 0],h1,ones(size(h1)),10,1,1,ones(2,lenh),1);
Ker1toc=toc;
disp(['LPA kernels created in ',num2str(Ker1toc),' seconds.   Total time: ',num2str(Ker1toc),' seconds.'])
for s2=1:lenh     % kernel size index
    gha=kernels_higher_order{1,s2,1}(:,:,1);   % gets single kernel from the cell array (ZERO ORDER)
    ghb=kernels_higher_orderb{1,s2,1}(:,:,1);  % gets single kernel from the cell array (FIRST ORDER)
    gh{s2}=(1+sharparams(s2))*ghb-sharparams(s2)*gha; % combines kernels into "order-mixture" kernel
    gh{s2}=single(gh{s2}((end+1)/2,(end+1)/2:end));
end
%---------------------------------------------------------
% Anisotropic LPA-ICI
%---------------------------------------------------------
h_opt_Q=function_AnisLPAICI8(single(zLumChrom(:,:,1)),gh,single(sigma),single(gammaICI));
%%%%%%%%%%%%%   END OF LPA-ICI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LPAtoc=toc;
disp(['LPA-ICI completed in ',num2str(LPAtoc-Ker1toc),' seconds.   Total time: ',num2str(LPAtoc),' seconds.'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SA-DCT HARD-THRESHOLDING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear trans_matrix Trian trans_matrix TrianRot_Int8 TrianRot_Int8W
% BUILDS DCT BASES (ALL SIZES ARE NEEDED)
for h=1:2*max(h_max,h_maxW)-1;
    trans_matrix{h}=dct(eye(h));
end
% BUILDS TRIANGLE MASKS FOR STARSHAPED SET  (kernel lengths as verteces)
for h_opt_1=h1
    for h_opt_2=h1
        Trian{h_opt_1,h_opt_2}=zeros(2*h_max-1);
        for i1=h_max-h_opt_2+1:h_max
            for i2=2*h_max-i1:(h_max-1+h_opt_1-(h_max-i1)*((h_opt_1-h_opt_2)/(h_opt_2-1+eps)))
                Trian{h_opt_1,h_opt_2}(i1,i2)=1;
            end
        end
    end
end
% BUILDS ROTATED TRIANGLE MASKS   (for the eight directions)
for ii=1:8
    for h_opt_1=h1
        for h_opt_2=h1
            if mod(ii,2)==0
                TrianRot_Int8{h_opt_1,h_opt_2,ii}=int8(rot90(-Trian{h_opt_2,h_opt_1}',mod(2+floor((ii-1)/2),8)));
            else
                TrianRot_Int8{h_opt_1,h_opt_2,ii}=int8(rot90(-Trian{h_opt_1,h_opt_2},mod(floor((ii-1)/2),8)));
            end
            if blocksize>0
                TrianRot_Int8{h_opt_1,h_opt_2,ii}([1:h_max-h1boundL, end-h_max+h1boundU+1:end],:)=int8(0);
                TrianRot_Int8{h_opt_1,h_opt_2,ii}(:,[1:h_max-h1boundL, end-h_max+h1boundU+1:end])=int8(0);
            end
        end
    end
end
LPA1toc=toc;
textlabel=('Pointwise SA-DCT hard-thresholding running ...');
disp(textlabel);
if (FiguresLumChromObs||FiguresLumChromEst)&&ReUseFiguresLumChrom
    LumChromFig=figure;
end
for colorchannel=[1 2 3];
    if FiguresLumChromObs
        if ReUseFiguresLumChrom==0;
            figure;
        else
            figure(LumChromFig)
        end
        imshow(zLumChrom(:,:,colorchannel));
        title(['Noisy observations, color channel #' num2str(colorchannel)]);
        drawnow
    end
    T=DCTthrCOEF*sqrt(2*log(SSS)+1)*chromDCTthrfactor(colorchannel);
    yLumChrom_hat(:,:,colorchannel) = function_SADCT_thresholding_fast(size_z_1, size_z_2, h_opt_Q, h_max, TrianRot_Int8, zLumChrom(:,:,colorchannel),0, sigma*l2normLumChrom(colorchannel)/l2normLumChrom(1), trans_matrix, h1, max_overlap(colorchannel),T,coef_align);   %%%% Pointwise SA-DCT hard-thresholding function
    if FiguresLumChromEst
        if ReUseFiguresLumChrom==0;
            figure;
        else
            figure(LumChromFig)
        end
        imshow(yLumChrom_hat(:,:,colorchannel));
        if compute_errors
            Errors = function_Errors(yLumChrom(:,:,colorchannel),yLumChrom_hat(:,:,colorchannel),zLumChrom(:,:,colorchannel));
            title(['Pointwise SA-DCT hard-thresholding est., color chnl.#' num2str(colorchannel),', PSNR=' num2str(Errors(3)) ', ISNR=' num2str(Errors(1))]);
        else
            title(['Pointwise SA-DCT hard-thresholding est., color chnl.#' num2str(colorchannel)]);
        end
        drawnow
    end
end
clear y_hat;
AVGtoc=toc;
disp(sprintf([repmat('\b',[1,numel(textlabel)+1]),'Pointwise SA-DCT hard-thresholding completed in ',num2str(AVGtoc-LPA1toc),' seconds.   Total time: ',num2str(AVGtoc),' seconds.']))
yRGB_hat=function_LumChrom2rgb(yLumChrom_hat,colormode);
if figures_y_hats
    figure
    imshow(yRGB_hat),
    if compute_errors
        Errors = function_Errors(yRGB,yRGB_hat,zRGB);
        title(['Pointwise SA-DCT hard-thresholding color estimate, PSNR=' num2str(Errors(3)) ', ISNR=' num2str(Errors(1))]);
    else
        title(['Pointwise SA-DCT hard-thresholding color estimate']);
    end
end
if do_wiener
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SA-DCT WIENER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     disp(' ');
    %     disp('---------------  Wiener Stage  ----------------------------------------------------');
    % BUILDS TRIANGLE MASKS FOR STARSHAPED SET  (kernel lengths as verteces)
    for h_opt_1=h1W
        for h_opt_2=h1W
            TrianW{h_opt_1,h_opt_2}=zeros(2*h_maxW-1);
            for i1=h_maxW-h_opt_2+1:h_maxW
                for i2=2*h_maxW-i1:(h_maxW-1+h_opt_1-(h_maxW-i1)*((h_opt_1-h_opt_2)/(h_opt_2-1+eps)))
                    TrianW{h_opt_1,h_opt_2}(i1,i2)=1;
                end
            end
        end
    end
    % BUILDS ROTATED TRIANGLE MASKS  (for the eight directions)
    for ii=1:8
        for h_opt_1=h1W
            for h_opt_2=h1W
                if mod(ii,2)==0
                    TrianRot_Int8W{h_opt_1,h_opt_2,ii}=int8(rot90(-TrianW{h_opt_2,h_opt_1}',mod(2+floor((ii-1)/2),8)));
                else
                    TrianRot_Int8W{h_opt_1,h_opt_2,ii}=int8(rot90(-TrianW{h_opt_1,h_opt_2},mod(floor((ii-1)/2),8)));
                end
                if blocksize>0
                    TrianRot_Int8W{h_opt_1,h_opt_2,ii}([1:h_maxW-h1boundL, end-h_maxW+h1boundU+1:end],:)=int8(0);
                    TrianRot_Int8W{h_opt_1,h_opt_2,ii}(:,[1:h_maxW-h1boundL, end-h_maxW+h1boundU+1:end])=int8(0);
                end
            end
        end
    end
    LPA2toc=toc;
    textlabel=('Pointwise SA-DCT Wiener-filter running ...');
    disp(textlabel);
    for colorchannel=[1 2 3];
        yLumChrom_hat_wi(:,:,colorchannel) = function_SADCT_wiener_fast(min(h_opt_Q,lenhW), TrianRot_Int8W, zLumChrom(:,:,colorchannel), sigma*l2normLumChrom(colorchannel)/l2normLumChrom(1), trans_matrix, h1W, yLumChrom_hat(:,:,colorchannel), max_overlapW(colorchannel),coef_align); %%%% Pointwise SA-DCT Wiener-filter function
        if FiguresLumChromEst
            if ReUseFiguresLumChrom==0;
                figure;
            else
                figure(LumChromFig)
            end
            imshow(yLumChrom_hat_wi(:,:,colorchannel));
            if compute_errors
                WErrors = function_Errors(yLumChrom(:,:,colorchannel),yLumChrom_hat_wi(:,:,colorchannel),zLumChrom(:,:,colorchannel));
                title(['Pointwise SA-DCT Wiener-filter est., color chnl.#' num2str(colorchannel),', PSNR=' num2str(WErrors(3)) ', ISNR=' num2str(WErrors(1))]);
            else
                title(['Pointwise SA-DCT Wiener-filter est., color chnl.#' num2str(colorchannel)]);
            end
            drawnow
        end

    end  %% colorchannels
    AVGWtoc=toc;
    disp(sprintf([repmat('\b',[1,numel(textlabel)+1]),'Pointwise SA-DCT Wiener-filter completed in ',num2str(AVGWtoc-LPA2toc),' seconds.   Total time: ',num2str(AVGWtoc),' seconds.']))
    yRGB_hat_wi=function_LumChrom2rgb(yLumChrom_hat_wi,colormode);
end
if (FiguresLumChromObs||FiguresLumChromEst)&&ReUseFiguresLumChrom
    close(LumChromFig)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DENOISING COMPLETED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% what follows is for the display of results %%%
if compute_errors
    for colorchannel=1:3;
        if do_wiener==1
            [error_criteria{colorchannel},Error_labels]=function_Errors(yRGB(:,:,colorchannel),yRGB_hat_wi(:,:,colorchannel),zRGB(:,:,colorchannel));
        else
            [error_criteria{colorchannel},Error_labels]=function_Errors(yRGB(:,:,colorchannel),yRGB_hat(:,:,colorchannel),zRGB(:,:,colorchannel));
        end
        error_criteria{colorchannel}=[repmat(' ',[size(error_criteria{colorchannel},1) 2]),num2str(error_criteria{colorchannel})];
    end
    if do_wiener==1
        error_yRGB_hat_final=function_Errors(yRGB,yRGB_hat_wi,zRGB);
    else
        error_yRGB_hat_final=function_Errors(yRGB,yRGB_hat,zRGB);
    end
    error_noisy=function_Errors(yRGB,zRGB);
    error_criteria{4}=error_yRGB_hat_final;
    error_criteria{4}=[repmat(' ',[size(error_criteria{4},1) 5]),num2str(error_criteria{4})];
    disp(' ');
    disp(' Error criteria:   (R, G, B, and overall)');
    disp([Error_labels,error_criteria{1},error_criteria{2},error_criteria{3},error_criteria{4}])
    disp(' ');
end
if figures_y_hats
    figure
    if do_wiener
        imshow(yRGB_hat_wi),
        if compute_errors
            title(['Pointwise SA-DCT Wiener-filter color estimate, PSNR=' num2str(error_yRGB_hat_final(3)) ', ISNR=' num2str(error_yRGB_hat_final(1))]);
        else
            title(['Pointwise SA-DCT Wiener-filter color estimate']);
        end
    else
        imshow(yRGB_hat),
        if compute_errors
            title(['Pointwise SA-DCT hard-thresholding color estimate, PSNR=' num2str(error_yRGB_hat_final(3)) ', ISNR=' num2str(error_yRGB_hat_final(1))]);
        else
            title(['Pointwise SA-DCT hard-thresholding color estimate']);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  END OF PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%