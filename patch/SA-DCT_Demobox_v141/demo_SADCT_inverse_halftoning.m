% Pointwise Shape-Adaptive DCT inverse-halftoning demo (using low-complexity SA-DCT)
%
% Alessandro Foi & Kostadin Dabov  - Tampere University of Technology - 2005 - 2006   Public release v1.40 (October 2006)
% -----------------------------------------------------------------------------------------------------------------------------------------------
%
% PURPOSE:  Reconstruct a continuous-tone grayscale image from an error-diffusion binary halftone
%
%
% Observation model:
%
% Error-diffusion halftoning is modeled as a convolutional process of the form
%           z  =  p * y  +  q * n
% where  z  is the halftone,  y  is the original grayscale image,  p  and  q  are
% convolutional kernels (which depend on the error-diffusion kernel), and  n  is
% Gaussian white noise.  Here  *  stands for the convolution operation.
% (see Kite et al., “Digital Image Halftoning as 2-D Delta-Sigma Modulation“,
%      Proc. IEEE ICIP, 1997.)
%
%
% Other key variables:
% zRI      : "unfiltered" regularized inverse estimate
% y_hat_RI : filtered RI estimate
% zRW      : "unfiltered" regularized Wiener inverse estimate
% y_hat_RW : filtered RW estimate (final estimate)
%
%
% The code implements the algorithm and reproduces the experiments published in:
% Dabov, K., A. Foi, V. Katkovnik, and K. Egiazarian, "Inverse halftoning by pointwise
% shape-adaptive DCT regularized deconvolution," Proc. Int. TICSP Workshop Spectral Meth.
% Multirate Signal Process., SMMSP 2006 , Florence, Italy, September 2006.
%

clear yh_RI yh_RW h_opt_Q stdh_RI stdh_RW gh hadper;
close all; drawnow;

% ----------------------------------------------------------------------
% INPUT SIGNAL SETTING
%-----------------------------------------------------------------------
load_from_file=0;     %% loads original noise-free bitmap image from user-selected file (recommended =0)
enable_crop=0;        %% enables cropping of the initial image before beginning of algorithm

% Original, continuous-tone image's filename to use if load_from_file==0
image_filename='image_Peppers512.png';
% image_filename='image_Boats512.png';
% image_filename='image_Lena512.png';

% Some options
err_kernel                = 2;         % Select an error diffusion filter:
%    1 -> Floyd-Steinberg
%    2 -> Jarvis et al.
do_wiener       = 1; % do_wiener = 0 does not perform RW stage ( = 1 performs, dafault )
detail_figures  = 0; % shows also figures with details of scales and estimates



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ALGORITHM PARAMETERS (it is recommended not to modify the following parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% RI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Regularization_epsilon_RI = 1.90;      % Regularization parameter of the Fourier domain inversion
h1RI                      = [1 2 3 4]; % The set of scales used for LPA in the RI
GammaParameterRI          = 1.44;      % Gamma parameter of LPA-ICI for RI
DCTthrCOEF                = 1.51;      % Empirically found coefficient that scales the universal threshold

%%% RW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Regularization_epsilon_RW = 0.326;       % Regularization parameter of the Fourier domain inversion
h1RW                      = [1 2 3 5]; % set of scales used for LPA in RWI
GammaParameterRW          = 3.23;      % Gamma parameter of LPA-ICI for RWI
DCTwieCOEF                = 1;         % Regularization parameter of the SA-DCT Wiener filtering

%%% Other auxiliary parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
small_size_RI             = 32;        % Size of the DFT used for the fast computation of the SA-DCT coefficients' variance in RI
small_size_RW             = 32;        % Size of the DFT used for the fast computation of the SA-DCT coefficients' variance in RWI

alphaorderRI              = -0.75;     % LPA order-mixture parameter in RI (-1 zero order, 0 first order)
alphaorderRW              = -0.8;      % LPA order-mixture parameter in RWI (-1 zero order, 0 first order)

sigma                     = 0.198;     % The standard deviation is fixed by the used error diffusion linear model

max_overlap=700;    %% limits overcompleteness (marginal speed-up)  (recommended >=49)  (note: lowering the value may give some speed-up in Pointwise SA-DCT RI hard-thresholding)
max_overlapW=700;   %% limits overcompleteness (marginal speed-up)  (recommended >=81)  (note: lowering the value may give some speed-up in Pointwise SA-DCT RW Wiener-filtering)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  INITIALIZATION AND GENERATION OF THE HALFTONE IMAGE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2RI   = ones(size(h1RI)); % To have line LPA kernels in RI
lenhRI = length(h2RI);     % number of scales in RI
h2RW   = ones(size(h1RW)); % To have line LPA kernels in RW
lenhRW = length(h1RW);     % number of scales in RW

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
if err_kernel==1
    error_kernel=[0 0 7; 3 5 1];  % Floyd-Steinberg
    K = 2.03;   % Gain constant
    method = 'Floyd-Steinberg';
elseif err_kernel==2
    error_kernel=[0 0 0 7 5; 3 5 7 5 3; 1 3 5 3 1];  % Jarvis et al.
    K = 4.45;   % Gain constant
    method = 'Jarvis et al.';
end

error_kernel = error_kernel/sum(error_kernel(:));  % normalizes error-kernel

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
figure
imshow(y), title('Original image')
drawnow


z = function_ErrorDiffusion(y,error_kernel);       % generates halftone
[size_z_1,size_z_2]=size(z);

%---------------------------------------------------------
% Display halftone observation
%---------------------------------------------------------
figure
imshow(z), title('Binary halftone observation')
drawnow


disp(' ')
disp(' ')
disp('-------------------------------------------------------------------------------------')
disp(' Pointwise Shape-Adaptive DCT inverse halftoning with deconvolution - Demo')
disp('-------------------------------------------------------------------------------------')
disp(' ')
disp(['Original image: ', image_filename(1:end-4), ' (', sprintf('%d',size_z_1), 'x', sprintf('%d',size_z_2),  ')'])
disp(['Error diffusion method: ', method])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Kite's error-diffusion convolutional model  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Z = PY + QN
h = zeros(size_z_1,size_z_2);
h(1:size(error_kernel,1),1:size(error_kernel,2)) = error_kernel;  %% zero-pads error-kernel
h = circshift(h,[0 -(size(error_kernel,2)-1)/2]);
H = fft2(h);
P = K./(1+(K-1)*H);
Q = (1-H)./(1+(K-1)*H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version -release; % get matlab release
matlab_R = str2num(ans);

ndirRI = 8; % number of directions
ndirRW = 8; % number of directions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    FILTERING STARTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Regularized Inverse   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  ')
disp('starting RI ...     ')
%%% LPA KERNELS
[kernels0, kernels_higher_order0]=function_CreateLPAKernels([0,0],h1RI,h2RI,10,1,1,ones(2,lenhRI),0);
[kernels1, kernels_higher_order1]=function_CreateLPAKernels([1,0],h1RI,h2RI,10,1,1,ones(2,lenhRI),0);

h_max                  = max(h1RI);
h1                     = h1RI;
lenh                   = lenhRI;
directional_resolution = ndirRI;
h_max = max(h1);

clear Trian TrianRot TrianRotu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILDS TRIANGLE MASKS FOR STARSHAPED SET  (h_opts as verteces)
for h_opt_1=h1RI
    for h_opt_2=h1RI
        Trian{h_opt_1,h_opt_2}=zeros(2*h_max-1);
        for i1=h_max-h_opt_2+1:h_max
            for i2=2*h_max-i1:(h_max-1+h_opt_1-(h_max-i1)*((h_opt_1-h_opt_2)/(h_opt_2-1+eps)))
                Trian{h_opt_1,h_opt_2}(i1,i2)=1;
            end
        end
    end
end
% BUILDS ROTATED TRIANGLE MASKS
for ii=1:8
    for h_opt_1=h1RI
        for h_opt_2=h1RI
            if mod(ii,2)==0
                TrianRot{h_opt_1,h_opt_2,ii}=rot90(Trian{h_opt_2,h_opt_1}',mod(2+floor((ii-1)/2),8));
                TrianRotu{h_opt_1,h_opt_2,ii}=uint8(rot90(Trian{h_opt_2,h_opt_1}',mod(2+floor((ii-1)/2),8)));
            else
                TrianRot{h_opt_1,h_opt_2,ii}=rot90(Trian{h_opt_1,h_opt_2},mod(floor((ii-1)/2),8));
                TrianRotu{h_opt_1,h_opt_2,ii}=uint8(rot90(Trian{h_opt_1,h_opt_2},mod(floor((ii-1)/2),8)));
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Transfer matrix for RI
RI = conj(P)./( (abs(P).^2) +(abs(Q)*Regularization_epsilon_RI).^2); % Transfer Matrix for RI, Standard Tikhonov Regularization
zRI = real(ifft2( fft2(z).* RI ));   % Regularized Inverse Estimate (RI OBSERVATION)
RIQ=RI .* Q;
riq=real(ifft2(RIQ));    % Impulse response of the RI times Q
%%  (approximate) square of RIQ modulus in downsampled Fourier domain
small_size = small_size_RI;
tran_fftmatrix = (ifft(small_size*eye(2*h_max-1),small_size))';
tran_fftmatrix_fast = tran_fftmatrix(:,1:end/2+1);
conj_fftmatrix = conj(tran_fftmatrix');
low_pass_k = ones(round(size_z_1/small_size),round(size_z_2/small_size));  %% low-pass kernel
low_pass_k = conv2(low_pass_k,ones(2));
low_pass_k = low_pass_k/(sum(low_pass_k(:)));
RIQsmallabs2 = conv2(repmat(abs(RIQ).^2,[3 3]),low_pass_k,'same');     %% low-pass filtering
RIQsmallabs2 = RIQsmallabs2(size_z_1+1+round((size_z_1/small_size)*(0:small_size-1)),size_z_2+1+round((size_z_2/small_size)*(0:small_size-1)));    %% decimation


%%%%%%%%%%%%%  ANISOTROPIC LPA-ICI  %%%%%%%%%%%%%%%%%%%
for s1=1:ndirRI % cycle along the directions
    for s=1:lenhRI,
        gh=(1+alphaorderRI)*kernels1{1,s}-alphaorderRI*kernels0{1,s}; % gets single kernel from the cell array
        if ~rem(s1,2) % rotates kernel 45 degrees if even direction
            for iziz=1:size(gh,2)
                gh(:,iziz)=circshift(gh(:,iziz),[(size(gh,2)+1)/2-iziz,0]);
            end
        end
        gh=rot90(gh,floor((s1-1)/2)); % rotates kernel 90 degrees every second direction
        %%%%%%%%%%%%%%%%%% LPA-RI  %%%%%%%%%%
        yh_RI(:,:,s)=conv2(zRI-100000,gh,'same')+100000;   %%% LPA Filtering (in spatial domain)
        stdh_RI(:,:,s)=repmat(sqrt(sum(sum(conv2(riq([size_z_1-h1RI(s)+2:size_z_1 1:size_z_1 1:h1RI(s)-1],[size_z_2-h1RI(s)+2:size_z_2 1:size_z_2 1:h1RI(s)-1]),gh,'valid').^2)))*sigma,[size_z_1,size_z_2]);   % STD of THE ESTIMATE LPA-RI
    end %%%%%%%%%%%%%%% end for H  %%%%%%%%%%%%%%%%%%%%

    %%%%%% ICI %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [YICI_RIT,h_optRI,std_optRI1]=function_ICI(yh_RI,stdh_RI,GammaParameterRI,2*(s1-1)*pi/8);
    h_opt_Q(:,:,s1)=h_optRI; %% STORES RESULTS OF ADAPTIVE SCALES (USED FOR ADAPTIVE NEIGHBORHOOD)
end   %%% END THETA LOOP
clear std_optRI1 YICI_RIT stdh_RI yh_RI;
if detail_figures
    h_optRI=int8(h_optRI);  %% kept only for visualization
else
    clear h_optRI;
end   %%%%%%%%%%%%% END OF ANISOTROPIC LPA-ICI %%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  SA-DCT PART STARTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SSS          = [1:(2*h_max-1)^2];
T            = DCTthrCOEF.*sqrt(2*log(SSS)+1);

% BUILDS DCT BASES (ALL SIZES ARE NEEDED)
for hh=1:2*h_max-1;
    hadper{hh}=dct(eye(hh));
end

% We process neighborhoods with identical shape one after another, simply by following the lexicographic order induced by the directional adaptive scales h_opt_Q which define the shape.
h_opt_Q = double(h_opt_Q);
h_opts_matrix=h_opt_Q(:,:,1);
for iiii=2:8
    h_opts_matrix=h_opts_matrix+h_opt_Q(:,:,iiii)*10^(iiii-1);
end
P2P=[[1:(size_z_1*size_z_2)]',h_opts_matrix(:)];
P2P=P2P(randperm(size_z_1*size_z_2),:);
P2P=sortrows(P2P,2);     %% sort according to lexicographic order

% SA-DCT hard-thresholding filtering
y_hat_RI = function_SADCT_RI_thresholding (h1, uint8(h_opt_Q), hadper, TrianRotu, conj_fftmatrix, tran_fftmatrix_fast, small_size, RIQsmallabs2, zRI, uint32(P2P),  sigma, DCTthrCOEF, T, h_max, max_overlap);
clear P2P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End of Regularized Inverse   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
figure, imshow(y_hat_RI), title('RI SA-DCT estimate'), drawnow

if do_wiener>=1    %%%%%%%%%%%%  RW Estimation %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Regularized Wiener   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    directional_resolution = ndirRW;
    lenh = lenhRW;
    h_max = max(h1RW);
    h1 = h1RW;
    %%% LPA KERNELS
    [kernels0, kernels_higher_order0]=function_CreateLPAKernels([0,0],h1RW,h2RW,10,1,1,ones(2,lenhRW),0);
    [kernels1, kernels_higher_order1]=function_CreateLPAKernels([1,0],h1RW,h2RW,10,1,1,ones(2,lenhRW),0);

    clear Trian TrianRot TrianRotu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BUILDS TRIANGLE MASKS FOR STARSHAPED SET  (h_opts as verteces)
    for h_opt_1=h1RW
        for h_opt_2=h1RW
            Trian{h_opt_1,h_opt_2}=zeros(2*h_max-1);
            for i1=h_max-h_opt_2+1:h_max
                for i2=2*h_max-i1:(h_max-1+h_opt_1-(h_max-i1)*((h_opt_1-h_opt_2)/(h_opt_2-1+eps)))
                    Trian{h_opt_1,h_opt_2}(i1,i2)=1;
                end
            end
        end
    end
    % BUILDS ROTATED TRIANGLE MASKS
    for ii=1:8
        for h_opt_1=h1RW
            for h_opt_2=h1RW
                if mod(ii,2)==0
                    TrianRot{h_opt_1,h_opt_2,ii}=rot90(Trian{h_opt_2,h_opt_1}',mod(2+floor((ii-1)/2),8));
                    TrianRotu{h_opt_1,h_opt_2,ii}=uint8(rot90(Trian{h_opt_2,h_opt_1}',mod(2+floor((ii-1)/2),8)));
                else
                    TrianRot{h_opt_1,h_opt_2,ii}=rot90(Trian{h_opt_1,h_opt_2},mod(floor((ii-1)/2),8));
                    TrianRotu{h_opt_1,h_opt_2,ii}=uint8(rot90(Trian{h_opt_1,h_opt_2},mod(floor((ii-1)/2),8)));
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%% Transfer matrix for RW
    Wiener_Pilot = abs(fft2(y_hat_RI));   %%% WIENER PILOT ESTIMATE
    PSD = (abs(Q)*sigma*sqrt(size_z_1*size_z_2)).^2;  %Power Spectrum of colored noise
    RW = conj(P).*Wiener_Pilot.^2./(Wiener_Pilot.^2.*(abs(P).^2) + PSD*(Regularization_epsilon_RW^2));   % Transfer Matrix for RW    %% Standard Regularization 'a-la-Tikhonov'
    zRW = real(ifft2(fft2(z).* RW ));   % RW OBSERVATION
    RWQ = RW.*Q;
    rwq = real(ifft2(RWQ));   % Impulse response of the RI times Q
    small_size=small_size_RW;
    tran_fftmatrix=(ifft(small_size*eye(2*h_max-1),small_size))';
    tran_fftmatrix_fast=tran_fftmatrix(:,1:end/2+1);
    conj_fftmatrix=conj(tran_fftmatrix');
    low_pass_k=ones(round(size_z_1/small_size),round(size_z_2/small_size));  %% low-pass kernel
    low_pass_k=conv2(low_pass_k,ones(2));
    low_pass_k=low_pass_k/(sum(low_pass_k(:)));
    %%  (approximate) square of RW modulus in downsampled Fourier domain
    RWQsmallabs2 = conv2(repmat(abs(RWQ).^2,[3 3]),low_pass_k,'same');     %% low-pass filtering
    RWQsmallabs2 = RWQsmallabs2(size_z_1+1+round((size_z_1/small_size)*(0:small_size-1)),size_z_2+1+round((size_z_2/small_size)*(0:small_size-1)));    %% decimation


    %%%%%%%%%%%%%  ANISOTROPIC LPA-ICI  %%%%%%%%%%%%%%%%%%%
    clear gh
    for s1=1:ndirRW  % cycle along the directions
        for s=1:lenhRW,
            gh=(1+alphaorderRW)*kernels1{1,s}-alphaorderRW*kernels0{1,s}; % gets single kernel from the cell array
            if ~rem(s1,2) % rotates kernel 45 degrees if even direction
                for iziz=1:size(gh,2)
                    gh(:,iziz)=circshift(gh(:,iziz),[(size(gh,2)+1)/2-iziz,0]);
                end
            end
            gh=rot90(gh,floor((s1-1)/2)); % rotates kernel 90 degrees every second direction
            %%%%%%%%%%%%%%%%%% LPA-RW  %%%%%%%%%%
            yh_RW(:,:,s)=conv2(zRW-100000,gh,'same')+100000;   %%% LPA Filtering (in spatial domain)
            stdh_RW(:,:,s)=repmat(sqrt(sum(sum(conv2(rwq([size_z_1-h1RW(s)+2:size_z_1 1:size_z_1 1:h1RW(s)-1],[size_z_2-h1RW(s)+2:size_z_2 1:size_z_2 1:h1RW(s)-1]),gh,'valid').^2)))*sigma,[size_z_1,size_z_2]);   % STD of THE ESTIMATE LPA-RW
        end %%%%%%%%%%%%%%% end for H  %%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%% ICI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [YICI_RWT,h_optRW,std_optRW1]=function_ICI(yh_RW,stdh_RW,GammaParameterRW,2*(s1-1)*pi/8);
        h_opt_Q(:,:,s1)=h_optRW;   %% STORES RESULTS OF ADAPTIVE SCALES (USED FOR ADAPTIVE NEIGHBORHOOD)
    end  % end theta
    %%%%%%%%%%%%% END OF ANISOTROPIC LPA-ICI %%%%%%%%%%%%%%%%
    clear h_optRW std_optRW1 YICI_RWT stdh_RW yh_RW;
    %%%%%%% SA-TRANSFORM FILTERING STARTS %%%%%%%%%%


    % BUILDS DCT BASES (ALL SIZES ARE NEEDED)
    clear hadper;
    for hh=1:2*h_max-1;
        hadper{hh}=dct(eye(hh));
    end

    % We process neighborhoods with identical shape one after another, simply by following the lexicographic order induced by the directional adaptive scales h_opt_Q which define the shape.
    h_opt_Q = double(h_opt_Q);
    h_opts_matrix=h_opt_Q(:,:,1);
    for iiii=2:8
        h_opts_matrix=h_opts_matrix+h_opt_Q(:,:,iiii)*10^(iiii-1);
    end
    P2Pw = [[1:(size_z_1*size_z_2)]',h_opts_matrix(:)];
    P2Pw = P2Pw(randperm((size_z_1*size_z_2)),:);
    P2Pw = sortrows(P2Pw,2);
    clear h_opts_matrix
    % SA-DCT Wiener filtering
    y_hat_RW = function_SADCT_RW_wiener (h1, uint8(h_opt_Q), hadper, TrianRotu, conj_fftmatrix, tran_fftmatrix_fast, small_size, RWQsmallabs2, zRW, uint32(P2Pw), sigma, T, h_max, y_hat_RI, DCTwieCOEF, max_overlapW);
    clear P2Pw
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End of Regularized Wiener   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RW_toc=toc;
    Err_RW=function_Errors(y,y_hat_RW,z);  %% computes error criteria
    %%% PRINT RESULTS TO SCREEN
    disp(sprintf(repmat('\b',[1 numel(tab_content)+size(tab_content,1)+numel(tab_title)+5+20]))),
    disp(['RW completed in ',num2str(RW_toc-RI_toc),' seconds.   Total time: ',num2str(RW_toc),' seconds.'])
    disp(' ');
    disp(tab_title);
    disp([Err_labels,repmat(' ',[size(Err_RI,1) 1]),num2str(Err_RI,number_of_digits),repmat(' ',[size(Err_RI,1) 2]),num2str(Err_RW,number_of_digits)]);
    disp(' ')
    figure, imshow(y_hat_RW);title('RW SA-DCT estimate')
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
if (do_wiener>=1)&&(detail_figures)&&(~enable_crop)
    for aaa=1:2 % HOW MANY DIFFERENT FIGURES (DIFFERENT DETAILS) TO SHOW
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


%
%       AN EXAMPLE OF THE RESULTS OBTAINED WITH THE CURRENT SOFTWARE
%
%
% -------------------------------------------------------------------------------------
%  Pointwise Shape-Adaptive DCT inverse halftoning with deconvolution - Demo
% -------------------------------------------------------------------------------------
%
% Original image: IMAGE_Peppers512 (512x512)
% Error diffusion method: Jarvis et al.
%
% RI completed in 37.4813 seconds.
% RW completed in 36.8902 seconds.   Total time: 74.3715 seconds.
%
%           RI        RW
% ISNR:  24.68341  24.97568
% SNR:   25.92868  26.22095
% PSNR:  31.68525  31.97752
% MSE:   44.11187  41.24088
% RMSE:  6.641677  6.421906
% MAE:   4.149942  3.997811
% MAX:    120.555  120.1357
