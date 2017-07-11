% Pointwise Shape-Adaptive DCT image denoising (using low-complexity SA-DCT)
%
% Alessandro Foi - Tampere University of Technology - 2005 - 2006   Public release v1.40 (October 2006)   C-code implementation by Kostadin Dabov
% -----------------------------------------------------------------------------------------------------------------------------------------------
%
%
%  IMPORTANT VARIABLES:
%
%  y         :  original image (if available)
%  z         :  noisy observation      (assumes AWGN model z = y + noise, where noise is iid gaussian)
%  y_hat     :  hard-thresholding estimate
%  y_hat_wi  :  Wiener-filter estimate (final estimate, obtained using y_hat as reference estimate in Wiener-filtering)
%  h_opt_Q   :  array with the adaptive-scale indexes for the kernels used in the construction of the adaptive-shape transform support
%
%
% The code implements the algorithm and reproduces the results published in:
% Foi, A., V. Katkovnik, and K. Egiazarian, “Pointwise Shape-Adaptive DCT for High-Quality Denoising and Deblocking
% of Grayscale and Color Images”, (accepted) IEEE Trans. Image Process., 2006.
% Foi, A., K. Dabov, V. Katkovnik, and K. Egiazarian, “Shape-Adaptive DCT for Denoising and Image Reconstruction”,
% Proc. SPIE Electronic Imaging 2006, Image Processing: Algorithms and Systems V, 6064A-18, San Jose, January 2006.
%
addpath('C:\Users\Iqbal\Documents\EE565\workspace\LASIP_Image_Restoration_DemoBox_v113');
clear trans_matrix Trian TrianRot_Int8 TrianW TrianRot_Int8W gh;
close all; drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Options for this demonstration software    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_from_file=0;     %% loads original noise-free bitmap image from user-selected file (recommended =1)
enable_crop=0;        %% enables cropping of the initial image before beginning of algorithm
sigma_noise=25/255;   %% standard-deviation of the added noise     (note: image range is assumed [0 1])
estimate_sigma=0;     %% estimate noise std from image?            (in this demo leave =0)
do_wiener=1;          %% enables Wiener-filter stage                          (recommended =1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Filtering performance vs. complexity  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
speedup_factor=1;   %% any number between 1 and 5. speedup_factor>1 enables some speed-ups in the algorithm
%                   %% algorithm will run about speedup_factor times faster at the expense of image quality
%                   %% (recommended =1)        (note: number can be fractional, e.g. speedup_factor=4.4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Visualization options  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figures_y_hats=1;  %% enables display of noisy image and SA-DCT estimates      (recommended =1)
compute_errors=1;  %% enables calculation of error criteria (PSNR, MSE, etc.)  (recommended =1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Alternative filenames which can be used directly if load_from_file==0  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_filename='../intro/cameraman.tif';
% image_filename='image_Cameraman256.png';
% image_filename='image_Peppers256.png';
% image_filename='image_Lena512.png';
% image_filename='image_Boats512.png';
% image_filename='image_House256.png';
% image_filename='image_Barbara512.png';
% image_filename='image_Cheese128.png';



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
max_overlap=70;         %% limits overcompleteness (speed-up)  (recommended =70)  (note: lowering the value leads to a considerable speed-up in Pointwise SA-DCT hard-thresholding)
max_overlapW=70;        %% limits overcompleteness (speed-up)  (recommended =70)  (note: lowering the value leads to a considerable speed-up in Pointwise SA-DCT Wiener-filtering)
if speedup_factor~=1   %%%  Speed-ups for "filtering performance vs. complexity" trade-off
    h1full=h1full(round(1:5/round((6-speedup_factor)):6));              %% complete set of scales to be used for hard-thresholding
    h1Wfull=h1Wfull(round(1:4/round((5-(4/5)*speedup_factor)):5));      %% complete set of scales to be used for Wiener-filter
    max_overlap=round((6-speedup_factor).^2.64);         %% limits overcompleteness (speed-up)  (note: lowering the value leads to a considerable speed-up in Pointwise SA-DCT hard-thresholding)
    max_overlapW=ceil(0.25*(6-speedup_factor).^3.5);     %% limits overcompleteness (speed-up)  (note: lowering the value leads to a considerable speed-up in Pointwise SA-DCT Wiener-filtering)
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
disp(' Pointwise Shape-Adaptive DCT Denoising   -   Public demo release')
disp('-----------------------------------------------------------------------------------')
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------
% Choose image from file (gets filename)
%---------------------------------------------------------
if load_from_file
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
y=im2double(imread(image_filename));

%---------------------------------------------------------
% Checks if color or grayscale image
%---------------------------------------------------------
if size(y,3)>1
    disp(' ')
    disp(['  !!!   This demo is for grayscale images only.  Please use  demo_SADCT_color_denoising.m  for color images   !!!'])
    disp(' ')
    disp(' ')
    return
end

%---------------------------------------------------------
% Cropping & original image display
%---------------------------------------------------------
if enable_crop
    display('---- PLEASE CROP REGION FROM IMAGE  ----')
    figure
    imshow(y);
    title('---- PLEASE CROP REGION FROM IMAGE  ----')
    hold on
    y=imcrop(y);
    hold off
    disp(sprintf([repmat('\b',[1,41]),' ']));
    close
end
if figures_y_hats
    figure
    imshow(y), title('Original image')
    drawnow
end

%---------------------------------------------------------
% Addition of noise (generates noisy observation)
%---------------------------------------------------------
init=0;
randn('seed', init);
z = y + sigma_noise*randn(size(y));     %%% z = y + noise
if estimate_sigma
    sigma=function_stdEst2D(z,2);    %%% estimates standard deviation from image (assumes perfect AWGN model)
else
    sigma=sigma_noise;
end

[size_z_1,size_z_2]=size(z);
if min(size_z_1,size_z_2)<2
    disp(' '); disp(['  !!!   Image is a vector,  exiting...  !!!']); disp(' '); disp(' ')
    return
end

%---------------------------------------------------------
% Display noisy observations and calculate errors
%---------------------------------------------------------
if compute_errors
    disp('Initial (noisy) observation criteria values:')
    function_Errors(y,z,z,2);
    disp(' ');
end
if figures_y_hats
    figure
    imshow(z), title('Noisy observation')
    drawnow
end

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
disp(['noise sigma = ',num2str(sigma*255),'    GammaICI = ',num2str(gammaICI), '    DCTthrCOEF = ', num2str(DCTthrCOEF), '    h_max = ',num2str(h_max),'   h_maxW = ',num2str(h_maxW)]);
if speedup_factor~=1
    disp(['   !!  Speed-up factor = ',num2str(speedup_factor), '  !!'])
end
if blocksize>0
    disp(['   !!  User-imposed constrained maximum block-size  ',num2str(blocksize),'x',num2str(blocksize), '  !!'])
end
disp(' ');
tic  %% algorithm begins

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
clear gh
% calling kernel creation function
[kernels, kernels_higher_order]=function_CreateLPAKernels([0 0],h1,ones(size(h1)),10,1,1,ones(2,lenh),1);
[kernelsb, kernels_higher_orderb]=function_CreateLPAKernels([1 0],h1,ones(size(h1)),10,1,1,ones(2,lenh),1);
for s2=1:lenh     % kernel size index
    gha=kernels_higher_order{1,s2,1}(:,:,1);   % gets single kernel from the cell array (ZERO ORDER)
    ghb=kernels_higher_orderb{1,s2,1}(:,:,1);  % gets single kernel from the cell array (FIRST ORDER)
    gh{s2}=(1+sharparams(s2))*ghb-sharparams(s2)*gha; % combines kernels into "order-mixture" kernel
    gh{s2}=single(gh{s2}((end+1)/2,(end+1)/2:end));
end
Ker1toc=toc;
disp(['LPA kernels created in ',num2str(Ker1toc),' seconds.   Total time: ',num2str(Ker1toc),' seconds.'])
%---------------------------------------------------------
% Anisotropic LPA-ICI
%---------------------------------------------------------
h_opt_Q=function_AnisLPAICI8(single(z),gh,single(sigma),single(gammaICI));   %%% Anisotropic LPA-ICI scales for 8 directions

LPAtoc=toc;
disp(['LPA-ICI completed in ',num2str(LPAtoc-Ker1toc),' seconds.   Total time: ',num2str(LPAtoc),' seconds.'])
%%%%%%%%%%%%%   END OF LPA-ICI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SA-DCT HARD-THRESHOLDING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------
% BUILDS DCT BASES (ALL SIZES ARE NEEDED)
%---------------------------------------------------------
for h=1:2*max(h_max,h_maxW)-1;
    trans_matrix{h}=dct(eye(h));
end
%---------------------------------------------------------
% BUILDS TRIANGLE MASKS FOR STARSHAPED SET  (kernel lengths as vertices)
%---------------------------------------------------------
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
%---------------------------------------------------------
% BUILDS ROTATED TRIANGLE MASKS  (for the eight directions)
%---------------------------------------------------------
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
%---------------------------------------------------------
% SA-DCT filtering
%---------------------------------------------------------
textlabel=('Pointwise SA-DCT hard-thresholding running ...');
disp(textlabel);
y_hat = function_SADCT_thresholding_fast(size_z_1, size_z_2, h_opt_Q, h_max, TrianRot_Int8, z,0, sigma, trans_matrix, h1, max_overlap,T,coef_align);  %%%% Pointwise SA-DCT hard-thresholding function
AVGtoc=toc;
disp(sprintf([repmat('\b',[1,numel(textlabel)+1]),'Pointwise SA-DCT hard-thresholding completed in ',num2str(AVGtoc-LPA1toc),' seconds.   Total time: ',num2str(AVGtoc),' seconds.']))
if compute_errors
    disp(' ');
    disp('Pointwise SA-DCT hard-thresholding results:');
    Errors = function_Errors(y,y_hat,z,2);
    disp(' ');
end
if figures_y_hats
    figure;
    imshow(y_hat);
    if compute_errors
        title(['Pointwise SA-DCT hard-thresholding estimate, PSNR = ' num2str(Errors(3)) ' dB, ISNR = ' num2str(Errors(1)),' dB']);
    else
        title(['Pointwise SA-DCT hard-thresholding estimate']);
    end
    drawnow
end

if do_wiener
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SA-DCT WIENER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------
    % BUILDS TRIANGLE MASKS FOR STARSHAPED SET (kernel lengths as vertices)
    %---------------------------------------------------------
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
    %---------------------------------------------------------
    % BUILDS ROTATED TRIANGLE MASKS  (for the eight directions)
    %---------------------------------------------------------
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
    %---------------------------------------------------------
    % SA-DCT filtering
    %---------------------------------------------------------
    textlabel=('Pointwise SA-DCT Wiener-filter running ...');
    disp(textlabel);
    y_hat_wi = function_SADCT_wiener_fast(min(h_opt_Q,lenhW), TrianRot_Int8W, z, sigma, trans_matrix, h1W, y_hat, max_overlapW,coef_align);   %%%% Pointwise SA-DCT Wiener-filter function
    AVGWtoc=toc;
    disp(sprintf([repmat('\b',[1,numel(textlabel)+1]),'Pointwise SA-DCT Wiener-filter completed in ',num2str(AVGWtoc-LPA2toc),' seconds.   Total time: ',num2str(AVGWtoc),' seconds.']))
    if compute_errors
        disp(' ');
        disp('Pointwise SA-DCT Wiener-filter results:');
        WErrors = function_Errors(y,y_hat_wi,z,2);
    end
    if figures_y_hats
        figure;
        imshow(y_hat_wi);
        if compute_errors
            title(['Pointwise SA-DCT Wiener-filter estimate, PSNR = ' num2str(WErrors(3)) ' dB,  ISNR = ' num2str(WErrors(1)),' dB']);
        else
            title(['Pointwise SA-DCT Wiener-filter estimate']);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DENOISING COMPLETED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  END OF PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%