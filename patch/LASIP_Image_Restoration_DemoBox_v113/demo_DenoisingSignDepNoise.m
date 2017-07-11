% Recursive Anisotropic LPA-ICI Denoising for Signal-Dependent Noise (demo_DenoisingSignDepNoise)
%
% Alessandro Foi - Tampere University of Technology - 2004-2005
% -----------------------------------------------------------------------
%
% Performs the recursive anisotropic LPA-ICI denoising on observations
% which are contaminated by signal-dependent noise (e.g. Poisson,
% Film-Grain, Speckle).
%
%
% This code implements the algorithm and replicates the results presented in
% Foi, A., R. Bilcu, V. Katkovnik, and K. Egiazarian,
% “Anisotropic local approximations for pointwise adaptive signal-dependent noise removal”,
% Proc. XIII European Signal Process. Conf., EUSIPCO 2005, September 2005.
%
% Simulations follow the experimental settings from
% R.M. Rangayyan, M. Ciuc, and F. Faghih,
% “Adaptive neighborhood filtering of images corrupted by signal-dependent noise,”
% Appl. Opt., vol. 37, no. 20, pp. 4477-4487, 1998.
%

clear all
close all

recursion_type=[1 1 1 2 2];               % Recursion type:   1=updates noise only (adaptive-variance)  2=recursive LPA-ICI filtering with adaptive-variance
gammaICIs=[0.6 0.6 0.6 0.8 0.8];          % ICI Gamma thresholds for the many iterations
sharparams=[-0.7 -0.6 -0.5 -0.2 -0.1];    % -1 zero order 0 first order (no sharpening) >0 sharpening
ndir=8;                                   % number of directions
fusing=1;                                 % fusing type   (1 classical fusing, 2 piecewise constant)
iterrule=1;                               % sigmaiter rule
itercoef=2/3;                             % compensating factor for recursive standard deviation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  TYPE OF SIGNAL-DEPENDENT NOISE TO CONSIDER  %%%%%%%%%%%%%%%%%%%%
%%%  uncomment only one of these  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 poisson=1;
% filmgrain=1;
% speckle=1;

%%% IMAGE
% y=im2double(imread('image_Lena512.png'));
% y=im2double(imread('image_Peppers512.png'));
 y=im2double(imread('image_Aerial256.png'));

% y=imresize(y,[256,256],'bilinear');  %% use this for comparing against some of the results (for 256x256 images) presented in Rangayyan's paper.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% LPA KERNELS SIZES
%--------------------------------------------------------------------------
h1=[1 2 3 5 8 11 15];

%h2=ones(size(h1));
h2=max(1,ceil(h1*tan(3*0.5*pi/ndir)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Recursive Anisotropic LPA-ICI Denoising Demo - Signal-Dependant Noise Case')
disp('-----------------------------------------------------------------------------------')
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  ')
if exist('poisson')
    if poisson==1
        speckle=0; filmgrain=0;
        disp(' Poissonian noise ')
    end
end
if exist('speckle')
    if speckle==1
        filmgrain=0; poisson=0;
        disp(' Speckle noise ')
    end
end
if exist('filmgrain')
    if filmgrain==1
        speckle=0; poisson=0;
        disp(' Film-grain noise ')
    end
end
disp('  ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% WINDOWS PARAMETERS
%--------------------------------------------------------------------------
sig_winds=[ones(size(h1))*0.6 ; ones(size(h2))*0.6];    % Gaussian parameter
beta=1;                     % Parameter of window 6

window_type=11;  % window=1 for uniform, window=2 for Gaussian
% window=11 for sectorial windows
% window=112 for sectorial uniform windows

TYPE=10;            % TYPE IS A SYMMETRY OF THE WINDOW
% 00 SYMMETRIC
% 10 NONSYMMETRIC ON X1 and SYMMETRIC ON X2
% 11 NONSYMMETRIC ON X1,X2  (Quadrants)

lenh=length(h1);
niter=min([length(recursion_type),length(sharparams),length(gammaICIs)]);


% init=0;%2055615866;
% rand('seed', init); randn('seed', init);
[size_z_1,size_z_2]=size(y);

if poisson==1
    chi=255/10;
    z=poissrnd(chi*y)/chi;
end

if filmgrain==1;
    K1=3.3; sigma1=1;    sigma2=0;  p=0.5;  %%Values from the experiments in the paper "adaptive neighborhood filtering of images corrupted by signal dependant noise"
    z=(255*y+(K1*(255*y).^p).*sigma1.*randn(size_z_1,size_z_2)+sigma2.*randn(size_z_1,size_z_2))/255;  %film grain noise model
end


if speckle==1;
    %    meanu=0; sigmau=1;
    z=zeros(size_z_1,size_z_2);
    for ii=1:4
        z=z+( (255*y).*exprnd(1,size(y,1),size(y,2)) )/255;  % SPECKLE NOISE model
    end
    z=(z/4);               % multi-shot speckle
end

disp('Initial (noisy) observation criteria values:')
function_Errors(y,z,z,2);
figure
imshow(z)
title('Noisy observation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%---------------------------------------------------------
% Kernels construction
%---------------------------------------------------------

[kernels, kernels_higher_order]=function_CreateLPAKernels([0 0],h1,h2,TYPE,window_type,ndir,sig_winds,beta);
[kernelsb, kernels_higher_orderb]=function_CreateLPAKernels([1 0],h1,h2,TYPE,window_type,ndir,sig_winds,beta);
disp(' ')
disp(' LPA kernels created');

z_iter=z;
if filmgrain==1;
    sigmaiter=(sqrt((K1^2)*(abs(z))*255*(sigma1^2)+sigma2^2))/255;
end
if speckle==1;
    sigmaiter=((sqrt((255*abs(z)).^2))/255)/2;
end
if poisson==1;
    sigmaiter=sqrt(abs(z).*(z>=0))/sqrt(chi);  %poisson noise
end
toc
disp('       ')
disp(['Running iteration ',num2str(1),' ...']);
figure
figure_number=gcf;
imshow(z), title('Noisy observation')
drawnow
htbar=timebar(['Recursive LPA-ICI running   (#iter=',num2str(niter),')'],'Progress');
for momo=1:niter
    if recursion_type(momo)==2
        z_iter=y_hat;
        sigmaiter=itercoef*sigmaiter1;
    end
    sharparam=sharparams(momo);
    gammaICI=gammaICIs(momo);
    clear yh h_opt_Q YICI_Q var_opt_Q stdh
    YICI_Final1=0; var_inv=0;         YICI_Final2=0;
    CWW=0;
    CWW2=0;
    %%%%% loops %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for s1=1:ndir     % directional index
        for s2=1:lenh     % kernel size index
            gha=kernels_higher_order{s1,s2,1}(:,:,1);   %gets single kernel from the cell array
            ghb=kernels_higher_orderb{s1,s2,1}(:,:,1);
            gh=(1+sharparam)*ghb-sharparam*gha;
            bound1=min([(find(sum(gh~=0,2)));abs(find(sum(gh~=0,2))-size(gh,1)-1)]); % removes unnecessary zeroes
            bound2=min([(find(sum(gh~=0,1))),abs(find(sum(gh~=0,1))-size(gh,2)-1)]); % removes unnecessary zeroes
            gh=gh(bound1:size(gh,1)-bound1+1,bound2:size(gh,2)-bound2+1);            % removes unnecessary zeroes
            ghorigin(s1,s2)=gh((end+1)/2,(end+1)/2);
            ziters{momo}=z_iter;
            if 1;%(momo==1)|(recursion_type(momo)==2)
                yh(:,:,s2)= conv2(z_iter+10000,gh,'same')-10000; % Estimation
                yhdd2{s1,s2,momo}=yh(:,:,s2);
                %                 if momo>1
                %                    max(max(abs(yh(:,:,s2)-yhdd(:,:,s2,s1))))
                %                 end
                %                 yhdd(:,:,s2,s1)=yh(:,:,s2);
            else
                yh(:,:,s2)=yhdd2{s1,s2,momo-1};
                yhdd2{s1,s2,momo}=yh(:,:,s2);
            end
            stdh(:,:,s2)=(conv2(sigmaiter.^2,gh.^2,'same')).^0.5;  % Std of the estimate
        end %% for s2, window sizes
        [YICI,h_opt,std_opt]=function_ICI(yh,stdh,gammaICI,2*(s1-1)*pi/ndir);    %%%% ICI %%%%%
        aaa=reshape(ghorigin(s1,h_opt),size(h_opt));
        YICI_Q(:,:,s1)=YICI;   h_opt_Q(:,:,s1)=h_opt;   var_opt_Q(:,:,s1)=(std_opt.^2+eps);

        YICI_Final1=YICI_Final1+YICI_Q(:,:,s1)./var_opt_Q(:,:,s1);            %% FUSING %%%%%
        YICI_Final2=YICI_Final2+YICI_Q(:,:,s1)./var_opt_Q(:,:,s1)-z_iter.*aaa./var_opt_Q(:,:,s1);            %% FUSING 2 %%%%%
        var_inv=var_inv+1./var_opt_Q(:,:,s1);
        CWW=CWW+aaa./var_opt_Q(:,:,s1);
        CWW2=CWW2+(aaa./var_opt_Q(:,:,s1)).^2;
        timebar(htbar,(s1+(momo-1)*ndir)/niter/ndir);
    end     %for s1 directions
    YICI_Final1=YICI_Final1./var_inv;
    YICI_Final2=(YICI_Final2+z_iter.*CWW/ndir)./(var_inv-CWW+CWW./ndir);           %% FUSING 2 %%%%%
    fusingmask=mean(h_opt_Q,3); fusingmask=fusingmask-min(fusingmask(:)); fusingmask=fusingmask/max(fusingmask(:));
    YICI_Final3=YICI_Final2.*fusingmask+YICI_Final1.*(1-fusingmask);
    if fusing==1
        y_hat=YICI_Final1;
    end

    if fusing==2
        y_hat=YICI_Final2;
    end
    if fusing==3
        y_hat=YICI_Final3;
    end
    %%%%%%%%%%%%%   END OF Filtering  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_hats{momo}=y_hat;
    if gcf~=figure_number, figure(figure_number), end
    imshow(y_hat), title(['LPA-ICI estimate (iter. ',num2str(momo),')']);
    [Err_y_hat,Err_labels]=function_Errors(y,y_hat,z);  %% computes error criteria
    Errors{momo}=Err_y_hat;
    if momo<niter
        Err_Wait=repmat('  wait... ',[size(Err_y_hat,1) 1]);
        size_Err_Wait=size(Err_Wait,2);
        Err_Wait=[repmat(' ',[1 size_Err_Wait]);Err_Wait];
    else
        Err_Wait=repmat([],[size(Err_y_hat,1)+1 1]);
    end
    if momo==1
        number_of_digits=4;  % number of digits for results in the table
        Err_y_hat_str=num2str(Err_y_hat,number_of_digits);
        tab_content=[['Iter.#',repmat(' ',[1 size(Err_labels,2)-6]);Err_labels],repmat(' ',[size(Err_y_hat,1)+1 1]),[repmat(' ',[1 floor((size(Err_y_hat_str,2))/2)]),num2str(momo),repmat(' ',[1,ceil((size(Err_y_hat_str,2))/2)-1]);Err_y_hat_str],Err_Wait];
    else
        Err_y_hat_str=num2str(Err_y_hat,number_of_digits);
        disp(sprintf(repmat('\b',[1 numel(tab_content)+size(tab_content,1)+27]))) % clears previous table
        tab_content=[tab_content(:,1:end-size_Err_Wait),repmat(' ',[size(Err_y_hat,1)+1 2]),[repmat(' ',[1 1-numel(num2str(momo))+floor((size(Err_y_hat_str,2))/2)]),num2str(momo),repmat(' ',[1,ceil((size(Err_y_hat_str,2))/2)-1]);Err_y_hat_str],Err_Wait];
    end
    if momo==1
        disp(sprintf(repmat('\b',[1 25]))) % clears previous table
    end
    if momo<niter
        disp(['Running iteration ',num2str(momo+1),' ...'])
    end
    if momo==niter
        disp('completed         ')
    end
    disp(' ');
    disp(tab_content);

    if filmgrain==1;
        sigmaiter=(sqrt((K1^2)*(abs(y_hat))*255*(sigma1^2)+sigma2^2))/255;
    end
    if speckle==1;
        sigmaiter=((sqrt((255*abs(y_hat)).^2))/255)/2;
    end
    if poisson==1;
        sigmaiter=sqrt(abs(y_hat).*(z>=0))/sqrt(chi);  %poisson noise
    end
    %         % %%% UPDATES NOISE BASED ON THE FUSING
    sigmaiter1=sqrt( (1./(var_inv.^2)).*(  var_inv -  (sigmaiter.^2).*CWW2 + (sigmaiter.*CWW).^2));    % standard fusing
end  %%% end of iteration
close(htbar);
close(figure_number);
disp(' ');
toc
figure
imshow(y_hat), title('Recursive LPA-ICI estimate');
version -release; % get matlab release
matlab_R=str2num(ans);

%%%% FIGURES %%%%%%%%%
figure
for iii=1:niter+2
    subplot(ceil((niter+2)/ceil(sqrt(niter+2))),ceil(sqrt(niter+2)),iii)
    if iii==1
        imshow(z)
        title('initial observation')

    else
        if iii<niter+2
            imshow(y_hats{iii-1});
            title(['ICI Estimate  iter=' , num2str(iii-1,3), '   ISNR=', num2str(Errors{iii-1}(1,:),3),'dB   ', num2str(recursion_type(iii-1),3)]);

        else
            imshow(y);
            title('original');

        end

    end
end


vscale=30;
screensize = get(0,'screensize');       % User's screen size [1 1 width height]
size_patch=min(round(min(screensize([3,3]),screensize([3,3]))/30),[size_z_1,size_z_2]);
for nfigures=1:2
    range1=(1+floor((size_z_1-1-size_patch(1))*rand))+[1:size_patch(1)];range2=(1+floor((size_z_2-1-size_patch(2))*rand))+[1:size_patch(2)];
    figure
    subplot(3,3,1);
    surf(flipud(vscale*y(range1,range2))),  axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
    caxis([0 vscale]);
    camproj('perspective');
    axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,      title('original');

    subplot(3,3,2);
    surf(flipud(vscale*z(range1,range2))),   axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
    camproj('perspective');
    caxis([0 vscale]);
    axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,     title('noisy');

    subplot(3,3,3);
    surf(flipud(vscale*y_hats{1}(range1,range2))),  axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
    camproj('perspective');
    caxis([0 vscale]);
    axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,     title('LPA-ICI estimate (1^s^t iter.)');

    subplot(3,3,4);
    surf(flipud(vscale*y_hats{2}(range1,range2))),   axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
    camproj('perspective');
    caxis([0 vscale]);
    axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,     title('LPA-ICI estimate, 2^n^d iter.');

    subplot(3,3,5);
    surf(flipud(vscale*y_hats{3}(range1,range2))), axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
    camproj('perspective');
    caxis([0 vscale]);
    axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,     title('LPA-ICI estimate, 3^r^d iter.');

    subplot(3,3,6);
    surf(flipud(vscale*y_hats{4}(range1,range2))), axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
    camproj('perspective');
    caxis([0 vscale]);
    axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,     title('LPA-ICI estimate, 4^t^h iter.');
    colormap((gray).^0.3)  %% use gamma correction for better visualizations of surfaces

    subplot(3,3,7);
    surf(flipud(vscale*y_hats{5}(range1,range2))), axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
    camproj('perspective');
    caxis([0 vscale]);
    axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,      title('LPA-ICI estimate, 5^t^h iter.');
    colormap((gray).^0.3)  %% use gamma correction for better visualizations of surfaces
    axes('position',[0.5 0.94 0.0001 0.0001]); axis off, title('Detailed view of a fragment of the image')
end  %% figures