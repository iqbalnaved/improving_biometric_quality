% Anisotropic Gradient Demo (demo_AnisotropicGradient)
%
% Alessandro Foi - Tampere University of Technology - 2003-2005
% -----------------------------------------------------------------------
%
% Demonstrates the Anisotropic Gradient concept using the Riemann surface
% example.
%
%
% This code implements the algorithm and reproduce the experiments published in:
% 
%

clear all
close all





h_max=6;  %% maximum length of kernel (any value between 3 and 10 will do) [DEFAULT h_max=6]
ndir=4;   %% number of directions (4 or 8 are typically enough)            [DEFAULT ndir=4]
m=1;      %% LPA order  (1 or 2; too large order gives noisier estimates)  [DEFAULT m=1]

gammaICI=2.5;   %% ICI Gamma-threshold                                     [DEFAULT gammaICI=2.5]














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  NO SIGNIFICANT PARAMETERS TO BE TUNED BELOW THIS LINE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%--------------------------------------------------------------------------
% LPA ORDER AND KERNELS SIZES
%--------------------------------------------------------------------------
m=[m,0];        % THE VECTOR ORDER OF LPA;
h1=[2:h_max];   % ALL SIZES BETWEEN 2 AND h_max
h2=ones(size(h1));  %% LINE KERNELS
lenh=length(h1);

%--------------------------------------------------------------------------
% WINDOWS PARAMETERS
%--------------------------------------------------------------------------
sig_winds=[ones(2,lenh)*0.6];    % Gaussian parameter
beta=1;                     % Parameter of window 6
window_type=1;  % window=1 for uniform, window=2 for Gaussian
TYPE=10;        % TYPE IS A SYMMETRY OF THE WINDOW
%               % 00 SYMMETRIC
%               % 10 NONSYMMETRIC ON X1 and SYMMETRIC ON X2
%               % 11 NONSYMMETRIC ON X1,X2  (Quadrants)

disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Anisotropic Gradient Demo ')
disp('-----------------------------------------------------------------------------------')
disp(' ')
disp(' Demonstrates the Anisotropic Gradient concept using the Riemann surface example.')
disp(' ')
disp(' ')
version -release; % get matlab release
matlab_R=str2num(ans);

%--------------------------------------------------------------------------
% CONSTRUCTION OF THE RIEMANN SURFACE EXAMPLE
%--------------------------------------------------------------------------
grid_length=200;  
x1=ones(grid_length,1)*[0:1/(grid_length-1):1]-0.5; x2=rot90(x1);
r=sqrt(x1.^2+x2.^2); 
phi=(angle(x1+i*x2)).*(r>0.1); %% RIEMANN SURFACE
%---------------------------------------------------------
% NOISE
%---------------------------------------------------------    
sigma_noise=0.01;  %%% std of the noise
[size_z_1,size_z_2]=size(phi);
init=0; %2055615866; 
randn('seed', init);
n=sigma_noise*randn(size(phi));
z = phi + n;  %% NOISE IS ADDED
sigma=function_stdEst2D(z);  %% NOISE STD IS ESTIMATED

%%% FIGURE OF OBSERVATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,2,1)
imshow(z,[]),colorbar
title('\phi observation')
subplot(1,2,2)
[meshgrid1,meshgrid2]=meshgrid([1:size_z_2],[size_z_1:-1:1]);plot3(meshgrid1(:),meshgrid2(:),z(:),' .','MarkerSize',0.5)
box on;        axis square tight,   camproj('perspective');  view([-43 26])
title('\phi observation')

%%% IDEAL ANALYTICAL GRADIENT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partx1_analytical=(-x2./(x1.^2+x2.^2))/200;
partx2_analytical=(x1./(x1.^2+x2.^2))/200;
partx1_analytical(find(r<=0.1))=0;
partx2_analytical(find(r<=0.1))=0;
Gradient_Analytical=partx1_analytical+i*partx2_analytical+eps*n;

%%% FIGURE OF IDEAL ANALYTICAL GRADIENT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,2,1)
imshow(partx1_analytical,[-0.05 0.05]); title('\partial\phi/\partialx_1'); colorbar
subplot(2,2,2)
imshow(partx2_analytical,[-0.05 0.05]); title('\partial\phi/\partialx_2'); colorbar
subplot(2,2,3)
imshow(abs(Gradient_Analytical),[0 0.1]); title('|\nabla\phi|'); colorbar
subplot(2,2,4)
imshow(abs(angle(Gradient_Analytical)),[0 pi]); title('\angle(\nabla\phi)'); colorbar
colormap(gray);
axes('position',[0.5 0.94 0.0001 0.0001]); axis off, title('Analytical "ideal" gradient \nabla\phi=(\partial\phi/\partialx_1,\partial\phi/\partialx_2)')

%---------------------------------------------------------
% LPA Kernels construction
%---------------------------------------------------------    
tic
disp('Creating LPA Kernels...    ')
[kernels, kernels_higher_order]=function_CreateLPAKernels(m,h1,h2,TYPE,window_type,ndir,sig_winds,beta);
disp(sprintf(['\b\b Kernels created in ',num2str(toc),' sec.']))
disp('Filtering (adaptive-scale directional derivatives) starts...    ')
tic  
%% Anisotropic Gradient Estimation Starts  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s1=1:ndir     % loop on directions
    for s2=1:lenh     % loop on kernel sizes
        gh=kernels_higher_order{s1,s2,1}(:,:,2);   %gets single kernel from the cell array
        bound1=min([(find(sum(gh~=0,2)));abs(find(sum(gh~=0,2))-size(gh,1)-1)]); % removes unnecessary zeroes
        bound2=min([(find(sum(gh~=0,1))),abs(find(sum(gh~=0,1))-size(gh,2)-1)]); % removes unnecessary zeroes
        gh=gh(bound1:size(gh,1)-bound1+1,bound2:size(gh,2)-bound2+1);            % removes unnecessary zeroes
        yh(:,:,s2)= conv2(z-100000000,gh,'same'); % Estimation
        stdh(:,:,s2)=repmat(sigma*(sum(gh(:).^2))^0.5,size_z_1,size_z_2);  % Std of the estimate
    end %% for s2, window sizes
    [YICI,h_opt,std_opt]=function_ICI(yh,stdh,gammaICI);     %%%% ICI for adaptive scale selection %%%%%
    std_opt(find(h_opt==1))=100000000000;    % Enforcing a derivability condition in the discrete domain
    if (stdh(1,1,2)+1/1000000000)>=stdh(1)       %%% this is needed when (e.g. in diagonal directions), because of rounding, the first and second smallest kernels coincide
        std_opt(find(h_opt==2))=100000000000;    % Enforcing a derivability condition in the discrete domain
    end
    DirDerAdapt(:,:,s1)=YICI;   h_opt_Q(:,:,s1)=h_opt;   var_opt_Q(:,:,s1)=(std_opt.^2+eps);
end %% for s1, directions
disp(sprintf(['\b\b Filtering completed in ',num2str(toc),' sec.'])) 
%%% ENDO OF FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Anisotropic gradient is computed (least-squares)...   ');
tic
THETASTEP=2*pi/ndir;
THETA=[0:THETASTEP:2*pi-THETASTEP];
rcos=reshape(repmat(cos(THETA)',1,size_z_1*size_z_2)',[size_z_1 size_z_2 size(THETA,2)]);
rsin=reshape(repmat(sin(THETA)',1,size_z_1*size_z_2)',[size_z_1 size_z_2 size(THETA,2)]);

%%% Least-Squares
lambda=1./var_opt_Q;
c1=sum(lambda.*rcos.*DirDerAdapt,3);
c2=sum(lambda.*rsin.*DirDerAdapt,3);
c3=sum(lambda.*rcos.*rsin,3);
c4=sum(lambda.*rsin.^2,3);
c5=sum(lambda.*rcos.^2,3);
partx1=(c4.*c1-c3.*c2)./(c4.*c5-c3.*c3);
partx2=(c2-partx1.*c3)./c4;    
Gradient_Anis=partx1+i*partx2;    
disp(sprintf(['\b\b Completed in ',num2str(toc),' sec.']'))
%RSS=sum(lambda.*(DirDerAdapt-rcos.*repmat(partx1,[1,1,size(THETA,2)])-rsin.*repmat(partx2,[1,1,size(THETA,2)])).^2,3);  %% Sum of Squared Residual

%%% FIGURE OF ANISOTROPIC GRADIENT
figure
subplot(2,2,1)
imshow(partx1,[-0.05 0.05]); title('\partial\phi/\partialx_1'); colorbar
subplot(2,2,2)
imshow(partx2,[-0.05 0.05]); title('\partial\phi/\partialx_2'); colorbar
subplot(2,2,3)
imshow(abs(Gradient_Anis),[0 0.1]); title('|\nabla\phi|'); colorbar
subplot(2,2,4)
imshow(abs(angle(Gradient_Anis)),[0 pi]); title('\angle(\nabla\phi)'); colorbar
colormap(gray);
axes('position',[0.5 0.94 0.0001 0.0001]); axis off, title('Anisotropic gradient \nabla\phi=(\partial\phi/\partialx_1,\partial\phi/\partialx_2)')

%%% UP TO FOUR FIGUREs OF STANDARD NON-ADAPTIVE GRADIENT ESTIMATES (FROM DIRECTIONAL DERIVATIVES)
for kernlength=1:ceil(round(h_max/2)/4):round(h_max/2)
    kernder1=([-kernlength:kernlength])';    kernder1=kernder1/([-kernlength:kernlength]*kernder1);
    partx1_NA=conv2(z+100000,-kernder1','same');
    partx2_NA=conv2(z+100000,kernder1,'same');
    clear i;
    Gradient_NA=partx1_NA+i*partx2_NA;
    figure
    subplot(2,2,1)
    imshow(partx1_NA,[-0.05 0.05]); title('\partial\phi/\partialx_1'); colorbar
    subplot(2,2,2)
    imshow(partx2_NA,[-0.05 0.05]); title('\partial\phi/\partialx_2'); colorbar
    subplot(2,2,3)
    imshow(abs(Gradient_NA),[0 0.1]); title('|\nabla\phi|'); colorbar
    subplot(2,2,4)
    imshow(abs(angle(Gradient_NA)),[0 pi]); title('\angle(\nabla\phi)'); colorbar
    colormap(gray);
    axes('position',[0.5 0.94 0.0001 0.0001]); axis off, title(['Gradient estimated using non-adaptive derivative kernels of length ',num2str(kernlength*2+1)])
end

%%% FIGURE OF ICI ADAPTIVE SCALES AND EDGES
figure
subplot(2,2,1)
imshow(h1(h_opt_Q(:,:,1)),[]), title(['Adaptive scales h^+, \theta=',num2str(THETA(1))]),colorbar
subplot(2,2,3)
imshow(h1(h_opt_Q(:,:,2)),[]), title(['Adaptive scales h^+, \theta=',num2str(THETA(2))]),colorbar
subplot(2,2,2)
imshow(sum(abs(DirDerAdapt),3),[0 0.5]), title('\Sigma_i|\partial\phi/\partial^+\theta_i,h^+|    (Edges)'),colorbar
subplot(2,2,4)
imshow(z,[]),title('\phi observation'), colorbar

%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')