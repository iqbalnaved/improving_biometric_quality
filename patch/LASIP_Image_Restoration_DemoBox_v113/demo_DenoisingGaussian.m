% Anisotropic LPA-ICI Denoising Demo (demo_DenoisingGaussian)
%
% Alessandro Foi - Tampere University of Technology - 2003-2005
% -----------------------------------------------------------------------
%
% Performs the anisotropic LPA-ICI denoising on observations which are
% contaminated by additive Gaussian White noise.
%
%
% Observation model:
% z=y+n
% z : noisy observation
% y : true image (assumed as unknown)
% n : Gaussian white noise
%
% Other key variables:
% y_hat   : anisotropic "fused" estimate
% y_hat_Q : adaptive directional estimate
% h_opt_Q : adaptive directional scales
%
%

clear all
close all

sharparam=-1;              % -1 zero order 0 first order (no sharpening) >0 sharpening
gammaICI=1.05;             % ICI Gamma threshold
directional_resolution=8;  % number of directions
fusing=1;                  % fusing type   (1 classical fusing, 2 piecewise regular)
addnoise=1;                % add noise to observation
sigma_noise=0.1;           % standard deviation of the noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% LPA KERNELS SIZES
%--------------------------------------------------------------------------
h1=[1 2 3 5 7 11];
h2=max(1,ceil(h1*tan(0.5*pi/directional_resolution)));  % row vectors h1 and h2 need to have the same lenght
%h2=ones(size(h1));
lenh=length(h1);

%--------------------------------------------------------------------------
% WINDOWS PARAMETERS
%--------------------------------------------------------------------------
sig_winds=[ones(size(h1)); ones(size(h2))];    % Gaussian parameter
beta=1;                     % Parameter of window 6

window_type=112;  % window=1 for uniform, window=2 for Gaussian
% window=6 for exponentions with beta
% window=8 for Interpolation
% window=11 for sectorial windows
% window=112 for sectorial unifrm windows

TYPE=10;            % TYPE IS A SYMMETRY OF THE WINDOW
% 00 SYMMETRIC
% 10 NONSYMMETRIC ON X1 and SYMMETRIC ON X2
% 11 NONSYMMETRIC ON X1,X2  (Quadrants)
%
% for rotated directional kernels the method that is used for rotation can be specified by adding
% a binary digit in front of these types, as follows:
%
% 10
% 11  ARE "STANDARD" USING NN (Nearest Neighb.) (you can think of these numbers with a 0 in front)
% 00
%
% 110
% 111  ARE EXACT SAMPLING OF THE EXACT ROTATED KERNEL
% 100
%
% 210
% 211  ARE WITH BILINEAR INTERP
% 200
%
% 310
% 311  ARE WITH BICUBIC INTERP (not reccomended)
% 300


%--------------------------------------------------------------------------
% MODELLING
%--------------------------------------------------------------------------

y=im2double(imread('image_Cameraman256.png'));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Anisotropic LPA-ICI Denoising Demo ')
disp('-----------------------------------------------------------------------------------')
disp(' ')


tic
%--------------------------------------------------------------------------
% MODELLING
%--------------------------------------------------------------------------

y=im2double(imread('image_Cameraman256.png'));
tic
init=0;%2055615866;
randn('seed', init);
%---------------------------------------------------------
% Images SIMULATION
%---------------------------------------------------------
[size_z_1,size_z_2]=size(y);

if addnoise==1
    n=repmat(sigma_noise,size_z_1,size_z_2).*randn(size(y));
    z = y + n;
else
    z = y;
end
sigma=function_stdEst2D(z);   %%% estimates noise standard deviation


disp('Initial (noisy) observation criteria values:')
function_Errors(y,z,z,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---------------------------------------------------------
% Kernels construction
%---------------------------------------------------------

% calling kernel creation function
[kernels, kernels_higher_order]=function_CreateLPAKernels([0 0],h1,h2,TYPE,window_type,directional_resolution,sig_winds,beta);
[kernelsb, kernels_higher_orderb]=function_CreateLPAKernels([1 0],h1,h2,TYPE,window_type,directional_resolution,sig_winds,beta);
disp('  ')
disp('kernels created');

sigmaiter=repmat(sigma,size_z_1,size_z_2);
stop_condition=0;
toc


clear yh h_opt_Q y_hat_Q var_opt_Q stdh
YICI_Final1=0; var_inv=0;         YICI_Final2=0;
CWW=0;
CWW2=0;


htbar=timebar('LPA-ICI running','Progress');
figure
figure_number=gcf;
imshow([z,zeros(size(z))])
disp(' ')
disp(' ')
disp('Directional estimates criteria values:')
disp(['  (Computing direction ',num2str(1),' ... )'])
%%%%% loops %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s1=1:directional_resolution     % directional index
    for s2=1:lenh     % kernel size index
        gha=kernels_higher_order{s1,s2,1}(:,:,1);   %gets single kernel from the cell array
        ghb=kernels_higher_orderb{s1,s2,1}(:,:,1);
        gh=(1+sharparam)*ghb-sharparam*gha;
        ghorigin(s1,s2)=gh((end+1)/2,(end+1)/2);
        bound1=min([(find(sum(gh~=0,2)));abs(find(sum(gh~=0,2))-size(gh,1)-1)]); % removes unnecessary zeroes
        bound2=min([(find(sum(gh~=0,1))),abs(find(sum(gh~=0,1))-size(gh,2)-1)]); % removes unnecessary zeroes
        gh=gh(bound1:size(gh,1)-bound1+1,bound2:size(gh,2)-bound2+1);            % removes unnecessary zeroes
        yh(1:size_z_1,1:size_z_2,s2)= conv2(z+10000,gh,'same')-10000; % Estimation
        stdh(:,:,s2)=repmat(sigma*(sum(gh(:).^2))^0.5,size_z_1,size_z_2);  % Std of the estimate
    end %% for s2, window sizes
    [YICI,h_opt,std_opt]=function_ICI(yh,stdh,gammaICI,2*(s1-1)*pi/directional_resolution);    %%%% ICI %%%%%
    aaa=reshape(ghorigin(s1,h_opt),size(h_opt));  %origin weight for optimal kernels
    y_hat_Q(:,:,s1)=YICI;   h_opt_Q(:,:,s1)=h_opt;   var_opt_Q(:,:,s1)=(std_opt.^2+eps);
    YICI_Final1=YICI_Final1+y_hat_Q(:,:,s1)./var_opt_Q(:,:,s1);            %% FUSING %%%%%
    YICI_Final2=YICI_Final2+y_hat_Q(:,:,s1)./var_opt_Q(:,:,s1)-z.*aaa./var_opt_Q(:,:,s1);            %% FUSING 2 %%%%%
    var_inv=var_inv+1./var_opt_Q(:,:,s1);
    CWW=CWW+aaa./var_opt_Q(:,:,s1);
    CWW2=CWW2+(aaa./var_opt_Q(:,:,s1)).^2;

    timebar(htbar,s1/directional_resolution);
    if gcf~=figure_number,        figure(figure_number), end
    imshow([YICI,(h_opt-1)/(lenh-1)]);
    title(['Directional adaptive estimates and corresponding adaptive scales (\theta=', num2str(360/directional_resolution*(s1-1)),'^o)'])
    [Err_y_Q,Err_labels]=function_Errors(y,YICI,z);  %% computes error criteria
    %%%% THIS PART HANDLES THE RESULTS TABLE DISPLAY
    if s1<directional_resolution
        Err_Wait=repmat('  wait... ',[size(Err_y_Q,1) 1]);
        size_Err_Wait=size(Err_Wait,2);
        Err_Wait=[repmat(' ',[1 size_Err_Wait]);Err_Wait];
    else
        Err_Wait=repmat([],[size(Err_y_Q,1)+1 1]);
    end
    if s1==1
        number_of_digits=4;  % number of digits for results in the table
        Err_y_Q_str=num2str(Err_y_Q,number_of_digits);
        tab_content=[['Dir.# ',repmat(' ',[1 size(Err_labels,2)-6]);Err_labels],repmat(' ',[size(Err_y_Q,1)+1 1]),[repmat(' ',[1 floor((size(Err_y_Q_str,2))/2)]),num2str(s1),repmat(' ',[1,ceil((size(Err_y_Q_str,2))/2)-1]);Err_y_Q_str],Err_Wait];
    else
        Err_y_Q_str=num2str(Err_y_Q,number_of_digits);
        disp(sprintf(repmat('\b',[1 numel(tab_content)+size(tab_content,1)+32]))) % clears previous table
        tab_content=[tab_content(:,1:end-size_Err_Wait),repmat(' ',[size(Err_y_Q,1)+1 2]),[repmat(' ',[1 1-numel(num2str(s1))+floor((size(Err_y_Q_str,2))/2)]),num2str(s1),repmat(' ',[1,ceil((size(Err_y_Q_str,2))/2)-1]);Err_y_Q_str],Err_Wait];
    end
    if s1==1
        disp(sprintf(repmat('\b',[1 32]))) % clears previous table
    end
    if s1<directional_resolution
        disp(['  (Computing direction ',num2str(s1+1),' ... )'])
    end
    if s1==directional_resolution
        disp('   ')
    end
    disp(tab_content);
    %%%% END OF RESULTS TABLE DISPLAY
end     %for s1 directions
disp(' ')
disp('Filtering completed.         ')
YICI_Final1=YICI_Final1./var_inv;
YICI_Final2=(YICI_Final2+z.*CWW/directional_resolution)./(var_inv-CWW+CWW./directional_resolution);           %% FUSING 2 %%%%%

if fusing==1, y_hat=YICI_Final1; end
if fusing==2, y_hat=YICI_Final2; end
%%%%%%%%%%%%%   END OF ALGORITHM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure_number)
close(htbar);

disp('  ')
disp('  ')
disp('Final anisotropic "fused" estimate criteria values:')
function_Errors(y,y_hat,z,2);
disp('  ')
toc
figure
subplot(1,2,1), imshow(h_opt_Q(:,:,1),[]), title(['Adaptive scales, \theta=', num2str(360/directional_resolution*(1-1)),'^o'])
subplot(1,2,2), imshow(y_hat_Q(:,:,1)), title(['Adaptive directional estimate, (\theta=', num2str(360/directional_resolution*(1-1)),'^o)'])

figure
subplot(1,2,1), imshow(h_opt_Q(:,:,4),[]), title(['Adaptive scales, \theta=', num2str(360/directional_resolution*(4-1)),'^o'])
subplot(1,2,2), imshow(y_hat_Q(:,:,4)), title(['Adaptive directional estimate, (\theta=', num2str(360/directional_resolution*(4-1)),'^o)'])
figure
subplot(1,2,1)
imshow(z); title('Noisy observation')
subplot(1,2,2)
imshow(y_hat); title('Anisotropic "fused" estimate')

