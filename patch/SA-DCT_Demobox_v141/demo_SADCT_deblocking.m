% Pointwise Shape-Adaptive DCT deblocking of Block-DCT compressed images (e.g. JPEG) (using low-complexity SA-DCT)
%
% Alessandro Foi - Tampere University of Technology - 2005 - 2006   Public release v1.40 (October 2006)   C-code implementation by Kostadin Dabov
% -----------------------------------------------------------------------------------------------------------------------------------------------
%
%  IMPORTANT VARIABLES:
%
%  yRGB         :  original image (if available) in RGB space
%  zRGB         :  B-DCT domain compressed image (observation) in RGB space
%  yRGB_hat     :  Pointwise SA-DCT final estimate in RGB space (final estimate, obtained from yLumChrom_hat_wi or yLumChrom_hat)
%  yLumChrom         :  original image (if available) in luminance-chrominance space (YUV)
%  zLumChrom         :  B-DCT domain compressed image (observation) in luminance-chrominance space (YUV)
%  yLumChrom_hat     :  hard-thresholding estimate in luminance-chrominance space (YUV)
%  yLumChrom_hat_wi  :  Wiener-filter estimate in luminance-chrominance space (YUV) (final estimate, obtained using yLumChrom_hat as reference estimate in Wiener-filtering)
%  h_opt_Q   :  array with the adaptive-scale indexes for the kernels used in the construction of the adaptive-shape transform support
%
%     Note:  if input image is grayscale it is used for the luminance channel Y and no RGB variables are used.
%
%
% The code implements the algorithm and reproduces the results published in:
% Foi, A., V. Katkovnik, and K. Egiazarian, “Pointwise Shape-Adaptive DCT for High-Quality Denoising and Deblocking of
% Grayscale and Color Images”, (accepted) IEEE Trans. Image Process., 2006.
% Foi, A., V. Katkovnik, and K. Egiazarian, “Pointwise Shape-Adaptive DCT for high-quality deblocking of compressed color images”,
% Proc. 14th European Signal Process. Conf., EUSIPCO 2006, Florence, September 2006.
%

clear trans_matrix Trian TrianRot_Int8 TrianW TrianRot_Int8W gh sigma_color yLumChrom_hat yLumChrom_hat_wi error_criteria zLumChrom meanQuant error_criteriaYUVRGB downsampling_rate;
close all; drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Options for this demonstration software    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
create_JPEG=1;           %%   0 -> loads JPEG user-selected file (no reference method)
%                        %%   1 -> create JPEG non-compressed bitmap  (using IJG-JPEG specifications)
%                        %%   2 -> performs B-DCT compression using predefined quantization tables
JPEG_Quality=10;         %%  JPEG-Quality factor, between 1 and 100 (used if create_JPEG==1)
Q_Table_QType='Q1';      %%   'Q1', 'Q2', or 'Q3'    (used if create_JPEG==2)
load_from_file=0;        %% loads original non-compressed bitmap image from user-selected file   (recommended =1)
enable_crop=0;           %% enables cropping of the initial image before beginning of algorithm
do_wiener=1;             %% enables Wiener-filter stage  (recommended =1)  (used for luminance channel only)

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
ReUseFiguresLumChrom=1;  %% reuse figure during channel-by-channel processing (avoids too many open figure windows)   (recommended =1)
figures_y_hats=1;        %% enables display of original, B-DCT compressed, and deblocked SA-DCT estimate image   (recommended =1)
compute_errors=1;        %% enables calculation of error criteria (PSNR, MSE, etc.)   (recommended =1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Alternative filenames which can be used directly if load_from_file==0  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% COLOR IMAGES %%%%%%%%%%%%%%%%
image_filename='image_House256rgb.png';
% image_filename='image_Peppers512rgb.png';
% image_filename='image_Lena512rgb.png';
% image_filename='image_Baboon512rgb.png';
% image_filename='image_F16_512rgb.png';
%%%%%%% GRAYSCALE IMAGES %%%%%%%%%%%%%
% image_filename='image_Lena512.png';
% image_filename='image_Cameraman256.png';
% image_filename='image_Peppers256.png';
% image_filename='image_GreenPeppers512.png';
% image_filename='image_Barbara512v2.png';
% image_filename='image_Peppers512.png';




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
%%  Deblocking/Deringing Custom parameters  %%%%    SHOULD NOT BE MODIFIED    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
structural_oversmoothing_factor=1.0;  %% increase this if JPEG is compressed from an already compressed original (e.g. jpeg of a manipulated jpeg)
textural_oversmoothing_factor=1.0;    %% increase this if image-data is noisy (e.g. JPEG taken with digital camera) 1.5 - 3.0 are reasonable values
%                                     %% visually pleasing results are obtained with some textural oversmoothing, say, 1.2 - 1.5
%                                     %% Restoration of cartoons typically require a rather larger textural oversmoothing ( 2.0 - 4.0 )
Quant_scale_factor_multiplier= 1.0;  %% MUST BE EQUAL TO 1 !!!!   DO NOT CHANGE IT UNLESS YOU KNOW WHAT YOU ARE DOING !!!! :)
%                                    %% it may be changed to a higher value (2 3, 10, 20, 30, ... up to 100 and beyond) in order to handle strong distortions due to over processing of a original image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Denoising algorithm's parameters  %%%%    SHOULD NOT BE MODIFIED   !!!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1full =  [1 2 3 5 7 9 12];   %% complete set of scales to be used for hard-thresholding      (recommended =[1 2 3 5 7 9 12])
h1Wfull = [1 2 3 5 6 8];      %% complete set of scales to be used for Wiener-filter          (recommended =[1 2 3 5 6 8])
sharparams=[-0.75 -0.70 -0.85 -0.90 -0.97 -1.00 -1.00];  %% order-mixture parameters (define the balance between zero- and first- order polynomial fitting for each scale)  (recommended =[-0.75 -0.70 -0.85 -0.9 -0.97 -1 -1])
DCTthrCOEF=0.925;           %% threshold-parameter  (recommended =0.925)
colormode='yuv';            %% luminance-chrominance space to use. 'dct' or 'opp' (opponent) or 'pca'  (recommended ='opp')
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
disp('------------------------------------------------------------------------------------------------------------------')
disp(' Pointwise Shape-Adaptive DCT Deblocking and Deringing for Block-DCT compressed images   -   Public demo release')
disp('------------------------------------------------------------------------------------------------------------------')
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------
% Choose image from file (gets filename)
%---------------------------------------------------------
if load_from_file&&create_JPEG~=0
    [bitmap_name,bitmap_path] = uigetfile('*.*','Select original image');
    if isequal(bitmap_name,0) || isequal(bitmap_path,0)
        disp(' ');
        disp('User pressed cancel');
        return
    else
        disp(['User selected   ', fullfile(bitmap_path,bitmap_name)]);
        disp(' ');
        image_filename=[bitmap_path, bitmap_name];
    end
elseif create_JPEG~=0
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

if create_JPEG==0
    [jpg_name,jpg_path] = uigetfile('*.jp*g','Select JPEG Image');
    if isequal(jpg_name,0) || isequal(jpg_path,0)
        disp(' ');
        disp('User pressed cancel');
        return
    else
        disp(' ');
        disp(['User selected   ', fullfile(jpg_path,jpg_name)]);
        disp(' ');
        image_filename=[jpg_path, jpg_name];
    end
end

%---------------------------------------------------------
% load image
%---------------------------------------------------------
yRGB=im2double(imread(image_filename));    %% note: at this point the loaded image can be either non-compressed or compressed bitmap

%---------------------------------------------------------
% Cropping
%---------------------------------------------------------
if enable_crop
    [size_z_1_full,size_z_2_full,size_z_3_full]=size(yRGB);   %% the sizes before cropping are used later in order to calculate downsampling rate and bpp
    display('---- PLEASE CROP REGION FROM IMAGE  ----')
    figure
    imshow(yRGB(:,:,1));
    title('---- PLEASE CROP REGION FROM IMAGE  ----')
    hold on
    yRGB=imcrop(yRGB);
    hold off
    disp(sprintf([repmat('\b',[1,41]),' ']));
    close
    [size_z_1,size_z_2,size_z_3]=size(yRGB);
else
    [size_z_1,size_z_2,size_z_3]=size(yRGB);
    size_z_1_full=size_z_1;
    size_z_2_full=size_z_2;
end

%---------------------------------------------------------
% Checks if color or grayscale image
%---------------------------------------------------------
if min(size_z_1,size_z_2)<2
    disp(' '); disp(['  !!!   Image is a vector,  exiting...  !!!']); disp(' '); disp(' ')
    return
end
if size_z_3==1
    image_is_grayscale=1;
    colorchannels=[1];
    yLumChrom=yRGB(:,:,1);
    clear yRGB;
else
    if create_JPEG==2
        disp(['  !!!   Predefined Quantization tables ("create_JPEG=2") can be used only with grayscale images  !!!'])
        disp(' ')
        disp(' ')
        return
    end
    image_is_grayscale=0;
    colorchannels=[1 2 3];
end
%---------------------------------------------------------
% Display original image
%---------------------------------------------------------
if figures_y_hats&create_JPEG~=0
    figure,
    if image_is_grayscale
        imshow(yLumChrom(:,:,1))
    else
        imshow(yRGB)
    end
    title(['Original image']);
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if create_JPEG==1    %%% Create JPEG image
    if image_is_grayscale
        imwrite(yLumChrom(:,:,1),'matlab_temporary_jpeg_SA_DCT.jpg','jpg','Quality',JPEG_Quality);
        zLumChrom=im2double(imread('matlab_temporary_jpeg_SA_DCT.jpg'));
        %             keyboard
    else
        imwrite(yRGB,'matlab_temporary_jpeg_SA_DCT.jpg','jpg','Quality',JPEG_Quality);
        zRGB=im2double(imread('matlab_temporary_jpeg_SA_DCT.jpg'));
    end
    image_filename='matlab_temporary_jpeg_SA_DCT.jpg';
elseif create_JPEG==0    %%% JPEG loaded from file
    if image_is_grayscale
        zLumChrom(:,:,1)=yLumChrom; clear yLumChrom  %% RENAME VARIABLE TO CORRECT NOTATION, USING z (and not y) to indicate compressed observation
    else
        zRGB=yRGB; clear yRGB;                       %% RENAME VARIABLE TO CORRECT NOTATION, USING z (and not y) to indicate compressed observation
    end
elseif create_JPEG==2   %%% performs B-DCT compression using predefined quantization tables
    if strcmp(Q_Table_QType,'Q1')         %% Q1 quantization table
        QTable_Y=[50 60 70 70 90 120 255 255;60 60 70 96 130 255 255 255;70 70 80 120 200 255 255 255;70 96 120 145 255 255 255 255;90 130 200 255 255 255 255 255;120 255 255 255 255 255 255 255;255 255 255 255 255 255 255 255;255 255 255 255 255 255 255 255];
    elseif strcmp(Q_Table_QType,'Q2')     %% Q2 quantization table
        QTable_Y=[86 59 54 86 129 216 255 255;64 64 75 102 140 255 255 255;75 70 86 129 216 255 255 255;75 91 118 156 255 255 255 255;97 118 199 255 255 255 255 255;129 189 255 255 255 255 255 255;255 255 255 255 255 255 255 255;255 255 255 255 255 255 255 255];
    elseif strcmp(Q_Table_QType,'Q3')     %% Q3 quantization table
        QTable_Y=[110 130 150 192 255 255 255 255;130 150 192 255 255 255 255 255;150 192 255 255 255 255 255 255;192 255 255 255 255 255 255 255;255 255 255 255 255 255 255 255;255 255 255 255 255 255 255 255;255 255 255 255 255 255 255 255;255 255 255 255 255 255 255 255];
    else
        disp(' '), disp(' '), disp(' !!!!!!!    Q_Table_QType  must be  Q1, Q2, or Q3    (used if create_JPEG==2)   !!!!!!! '),  disp(' '), disp(' ')
        return
    end
    %%% B-DCT compression
    DCT_matrix8=dct(eye(8));
    iDCT_matrix8=DCT_matrix8'; % =inv(DCT_matrix8);
    zLumChrom=zeros(size_z_1,size_z_2);
    for i1=[1:8:size_z_1-8,size_z_1-8+1],
        for i2=[1:8:size_z_2-8,size_z_2-8+1]
            zBLOCK=yLumChrom(i1:i1+7,i2:i2+7,1);
            win=DCT_matrix8*zBLOCK*iDCT_matrix8;           % forward transform
            winh=round(255*win./QTable_Y).*QTable_Y/255;   % quantization
            qBLOCK=iDCT_matrix8*winh*DCT_matrix8;          % inverse transform
            zLumChrom(i1:i1+7,i2:i2+7)=qBLOCK;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  ALGORITHM BEGINS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMAGE ANALYSIS (Downsampling & QTables) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if create_JPEG==1||create_JPEG==0
    [fileinfo_temp]=dir(image_filename);fileinfo_temp=struct2cell(fileinfo_temp);bpp_jpeg=fileinfo_temp{3}*8/(size_z_1_full*size_z_2_full);
    JPEG_header_info = jpeg_read(image_filename);   %% get information about JPEG file
    %% downsampling rate
    for colorchannel=colorchannels
        if ceil([size_z_1_full size_z_2_full]/8)*8==size(JPEG_header_info.coef_arrays{colorchannel})
            downsampling_rate(colorchannel)=1;
        elseif ceil([size_z_1_full size_z_2_full]/16)*8==size(JPEG_header_info.coef_arrays{colorchannel})
            downsampling_rate(colorchannel)=2;
        else
            downsampling_rate(colorchannel)=ceil(max([size_z_1_full size_z_2_full]./size(JPEG_header_info.coef_arrays{colorchannel})));
            disp(' ')
            disp(['  !!! WARNING:  Unusual downsampling rate,  size of coded coefficients array does not match with size of image   !!!'])
            disp(' ')
            disp(' ')
        end
    end
    %% Check if image is grayscale
    if 1==JPEG_header_info.image_components;
        image_is_grayscale=1;
        colorchannels=[1];
    end
    %% Quantization tables
    QTables=JPEG_header_info.quant_tables;
    clear JPEG_header_info;
    for colorchannel=colorchannels
        meanQuant(colorchannel)=mean(mean(QTables{min(numel(QTables),colorchannel)}(1:3,1:3)));
    end

    if ~image_is_grayscale
        if strcmp(colormode,'pca')  %% gets color-transformation matrix for PCA
            [zLumChrom colormode]=function_rgb2LumChrom(zRGB,colormode);
        else
            zLumChrom=function_rgb2LumChrom(zRGB,colormode);
        end
        if create_JPEG~=0
            yLumChrom=function_rgb2LumChrom(yRGB,colormode);
        end
    end
end
if create_JPEG==2
    disp(['Performing ad-hoc B-DCT using quant. table ', Q_Table_QType]);
    disp([' first row:   ',Q_Table_QType,'(1,1:8) = [ ' num2str(QTable_Y(1,1:8)),' ]']);
    disp(' ');
    meanQuant=mean(mean(QTable_Y(1:3,1:3)));
    downsampling_rate=[1];  %% no downsampling
end
if (create_JPEG==1||create_JPEG==2)&&compute_errors
    disp('Initial (JPEG) observation criteria values:')
    if image_is_grayscale
        error_noisy=function_Errors(yLumChrom(:,:,1),zLumChrom(:,:,1),2);
    else
        error_noisy=function_Errors(yRGB,zRGB,2);
    end
    disp(' ');
end
%---------------------------------------------------------
% Display observations (compressed image)
%---------------------------------------------------------
if figures_y_hats
    figure
    if image_is_grayscale
        imshow(zLumChrom(:,:,1))
    else
        imshow(zRGB)
    end
    if create_JPEG==1
        if compute_errors
            title(['Initial data (observations),  bits-per-pixel (bpp) = ',num2str(bpp_jpeg),',  PSNR = ' num2str(error_noisy(2)),' dB']);
        else
            title(['Initial data (observations),  bits-per-pixel (bpp) = ',num2str(bpp_jpeg)]);
        end
    elseif create_JPEG==0
        title(['Initial data (observations),  bits-per-pixel (bpp) = ',num2str(bpp_jpeg)]);
    elseif create_JPEG==2
        if compute_errors
            title(['Initial data (observations),   PSNR = ' num2str(error_noisy(2)),' dB']);
        else
            title(['Initial data (observations)']);
        end
    end
    drawnow
end
%%%% Models B-DCT compression as Additive Gaussian Noise with the following sigmas
sigmaLumChrom=((meanQuant*7.65*Quant_scale_factor_multiplier).^0.65)/1150;
sigmaLumChrom=sigmaLumChrom.*sqrt(downsampling_rate);   %% for chrominances multiplies sigma by sqrt(2) in order to address downsampling
sigma=sigmaLumChrom(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gammaICI=max(0.8,structural_oversmoothing_factor*2.4./(log(1+60*sigma)));    %%% the maximum is used in order to ensure large enough gamma for strong levels of noise (above sigma=75)


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
if create_JPEG==1
    disp(['Quality = ',num2str(JPEG_Quality), '  Quant. = [',num2str(num2str(meanQuant,4)),']  bpp = ',num2str(bpp_jpeg)]);
elseif create_JPEG==0
    disp(['JPEG from user-selected file','  Quant. = [',num2str(num2str(meanQuant,4)),']  bpp = ',num2str(bpp_jpeg)]);
elseif create_JPEG==2
    disp(['Quant.Table: ',Q_Table_QType,'  Quant. = [',num2str(num2str(meanQuant,4))]);
end
disp(['sigma = ',num2str(sigma*255),'  GammaICI = ',num2str(gammaICI),'  DCTthrCOEF = ',num2str(DCTthrCOEF),'    h_max = ',num2str(h_max),'   h_maxW = ',num2str(h_maxW)]);
disp(' ');
if (structural_oversmoothing_factor~=1)||(textural_oversmoothing_factor~=1);
    disp(['    structural_oversmoothing_factor = ',num2str(structural_oversmoothing_factor),'    textural_oversmoothing_factor = ',num2str(textural_oversmoothing_factor)]); disp(' ');
end
if Quant_scale_factor_multiplier~=1
    disp(['    Quant_scale_factor_multiplier = ',num2str(Quant_scale_factor_multiplier),'       !!!!!!!!!!!!!!!!!!!!!!!!!!']);        disp(' ');
end
if image_is_grayscale
    disp(['  Image is grayscale  ( used for Luminance channel only ) ']);        disp(' ');
end
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
textlabel=('Pointwise SA-DCT hard-thresholding running ...');
disp(textlabel);
if (FiguresLumChromObs||FiguresLumChromEst)&&ReUseFiguresLumChrom
    LumChromFig=figure;
end
for colorchannel=colorchannels;
    if FiguresLumChromObs&&~image_is_grayscale
        if ReUseFiguresLumChrom==0;
            figure;
        else
            figure(LumChromFig)
        end
        imshow(zLumChrom(:,:,colorchannel));
        title(['Initial data (observations), color channel #' num2str(colorchannel)]);
        drawnow
    end
    T=DCTthrCOEF*sqrt(2*log(SSS)+1)*chromDCTthrfactor(colorchannel)*textural_oversmoothing_factor;
    yLumChrom_hat(:,:,colorchannel)=function_SADCT_thresholding_fast(size_z_1, size_z_2, h_opt_Q, h_max, TrianRot_Int8,zLumChrom(:,:,colorchannel),0, sigmaLumChrom(colorchannel), trans_matrix, h1, max_overlap(colorchannel),T,coef_align);
    if FiguresLumChromEst
        if ReUseFiguresLumChrom==0;
            figure;
        else
            figure(LumChromFig)
        end
        imshow(yLumChrom_hat(:,:,colorchannel));
        if (create_JPEG==1||create_JPEG==2)&&compute_errors
            Errors = function_Errors(yLumChrom(:,:,colorchannel),yLumChrom_hat(:,:,colorchannel),zLumChrom(:,:,colorchannel));
            title(['Pointwise SA-DCT hard-thresholding est., color chnl.#' num2str(colorchannel), ', PSNR = ', num2str(Errors(3)), ' dB,  ISNR = ' num2str(Errors(1)),' dB']);
        else
            title(['Pointwise SA-DCT hard-thresholding estimate, color channel #' num2str(colorchannel)]);
        end
        drawnow
    end
end
AVGtoc=toc;
disp(sprintf([repmat('\b',[1,numel(textlabel)+1]),'Pointwise SA-DCT hard-thresholding completed in ',num2str(AVGtoc-LPAtoc),' seconds.   Total time: ',num2str(AVGtoc),' seconds.']))
if image_is_grayscale
else
    yRGB_hat=function_LumChrom2rgb(yLumChrom_hat,colormode);
end
if do_wiener
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SA-DCT WIENER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    for colorchannel=[1];
        yLumChrom_hat_wi(:,:,colorchannel) = function_SADCT_wiener_fast(min(h_opt_Q,lenhW), TrianRot_Int8W, zLumChrom(:,:,colorchannel), sigmaLumChrom(colorchannel), trans_matrix, h1W, yLumChrom_hat(:,:,colorchannel), max_overlapW(colorchannel),coef_align); %%%% Pointwise SA-DCT Wiener-filter function
        if FiguresLumChromEst
            if ReUseFiguresLumChrom==0;
                figure;
            else
                figure(LumChromFig)
            end
            imshow(yLumChrom_hat_wi(:,:,colorchannel));
            if (create_JPEG==1||create_JPEG==2)&&compute_errors
                WErrors = function_Errors(yLumChrom(:,:,colorchannel),yLumChrom_hat_wi(:,:,colorchannel),zLumChrom(:,:,colorchannel));
                title(['Pointwise SA-DCT Wiener-filter est., color chnl.#' num2str(colorchannel), ',  PSNR = ' num2str(WErrors(3)) ' dB,  ISNR = ' num2str(WErrors(1)),' dB']);
            else
                title(['Pointwise SA-DCT Wiener-filter estimate, color channel #' num2str(colorchannel)]);
            end
            drawnow
        end
    end  %% colorchannels
    AVGWtoc=toc;
    disp(sprintf([repmat('\b',[1,numel(textlabel)+1]),'Pointwise SA-DCT Wiener-filter completed in ',num2str(AVGWtoc-LPA2toc),' seconds.   Total time: ',num2str(AVGWtoc),' seconds.']))
end
if (FiguresLumChromObs||FiguresLumChromEst)&&ReUseFiguresLumChrom
    close(LumChromFig)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DEBLOCKING COMPLETED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------
% go back to RGB space (if image is color)
%---------------------------------------------------------
if ~image_is_grayscale
    if do_wiener
        if size(yLumChrom_hat_wi,3)==3    %%% Wiener was performed also for chrominances
            yRGB_hat=function_LumChrom2rgb(yLumChrom_hat_wi,colormode);
        else                              %%% Wiener was performed only to lumincance
            yRGB_hat=function_LumChrom2rgb(cat(3,yLumChrom_hat_wi(:,:,1),yLumChrom_hat(:,:,2),yLumChrom_hat(:,:,3)),colormode);
        end
    else
        yRGB_hat=function_LumChrom2rgb(yLumChrom_hat,colormode);
    end
end

%---------------------------------------------------------
% computes & displays error criteria
%---------------------------------------------------------
if (create_JPEG==1||create_JPEG==2)&&compute_errors
    disp(' ');
    for colorchannel=colorchannels;  %% computes error criteria for YUV
        if do_wiener
            if colorchannel>size(yLumChrom_hat_wi,3)
                [error_criteriaYUVRGB{colorchannel},Error_labels]=function_Errors(yLumChrom(:,:,colorchannel),yLumChrom_hat(:,:,colorchannel),zLumChrom(:,:,colorchannel));
            else
                [error_criteriaYUVRGB{colorchannel},Error_labels]=function_Errors(yLumChrom(:,:,colorchannel),yLumChrom_hat_wi(:,:,colorchannel),zLumChrom(:,:,colorchannel));
            end
        else
            [error_criteriaYUVRGB{colorchannel},Error_labels]=function_Errors(yLumChrom(:,:,colorchannel),yLumChrom_hat(:,:,colorchannel),zLumChrom(:,:,colorchannel));
        end
        error_criteriaYUVRGB{colorchannel}=[repmat(' ',[size(Error_labels,1) 2]),num2str(error_criteriaYUVRGB{colorchannel})];
    end
    if ~image_is_grayscale  %% computes error criteria for RGB
        for colorchannel=1:3   % RGB
            [error_criteriaYUVRGB{colorchannel+3},Error_labels]=function_Errors(yRGB(:,:,colorchannel),yRGB_hat(:,:,colorchannel),zRGB(:,:,colorchannel));
            error_criteriaYUVRGB{colorchannel+3}=[repmat(' ',[size(Error_labels,1) 2]),num2str(error_criteriaYUVRGB{colorchannel+3})];
        end
        error_yRGB_hat=function_Errors(yRGB,yRGB_hat,zRGB);
        error_criteriaYUVRGB{7}=[repmat(' ',[size(Error_labels,1) 5]),num2str(error_yRGB_hat)];
        disp(' Error criteria:   (Y, U, V,  R, G, B,  and overall)');
        disp([Error_labels,error_criteriaYUVRGB{1},error_criteriaYUVRGB{2},error_criteriaYUVRGB{3},repmat(' ',[size(Error_labels,1) 2]),error_criteriaYUVRGB{4},error_criteriaYUVRGB{5},error_criteriaYUVRGB{6},error_criteriaYUVRGB{7}])
    else
        disp(' Error criteria  (for luminance only - image is grayscale) :');
        disp([Error_labels,repmat(' ',[size(Error_labels,1) 2]),num2str(error_criteriaYUVRGB{1})]);
    end
    disp(' ');
end

%---------------------------------------------------------
% figure of final estimate
%---------------------------------------------------------
if figures_y_hats
    figure
    if image_is_grayscale
        if do_wiener
            imshow(yLumChrom_hat_wi(:,:,1));
        else
            imshow(yLumChrom_hat(:,:,1));
        end
    else
        imshow(yRGB_hat);
    end
    if (create_JPEG==1||create_JPEG==2)&&compute_errors
        if image_is_grayscale
            title(['Pointwise SA-DCT final estimate,  PSNR = ', error_criteriaYUVRGB{1}(3,:), ' dB,  ISNR = ', error_criteriaYUVRGB{1}(1,:),' dB']);
        else
            title(['Pointwise SA-DCT final estimate,  PSNR = ', num2str(error_yRGB_hat(3)), ' dB,  ISNR = ', num2str(error_yRGB_hat(1)),' dB']);
        end
    else
        title(['Pointwise SA-DCT final estimate.']);
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