% Recursive Anisotropic LPA-ICI Denoising Demo (demo_RecursiveDenoisingGaussian)
%
% Alessandro Foi - Tampere University of Technology - 2003-2005
% -----------------------------------------------------------------------
%
% Performs the recursive anisotropic LPA-ICI denoising on observations
% which are contaminated by additive Gaussian White noise.
%
%
% Observation model:
% z=y+n
% z : noisy observation
% y : true image (assumed as unknown)
% n : Gaussian white noise
%
%
% Other key variables:
% y_hats  : cell containing all the recursive anistropic estimates
% y_hat   : anisotropic estimate obtained after the last iteration
%
%
%
% This code implements the algorithm and replicates the results presented in
% Foi, A., V. Katkovnik, K. Egiazarian, and J. Astola,
% “A novel anisotropic local polynomial estimator based on directional multiscale optimizations”,
% Proc. of the 6th IMA Int. Conf. Math. in Signal Processing, Cirencester (UK), pp. 79-82, 2004.
%

clear all
close all
%clc

if 0       %% put 0 to use High-Quality parameters (slower);
           %% put 1 to use shorter kernels (faster, but lower quality) [as in the referenced paper].

    %%% FAST, LOW QUALITY SETTINGS  (with redundant iteration numbers to show convergence)
    h1=[1 2 3 5];                   % LPA Kernels' length
    sharparams=-1;                  % -1 zero order 0 first order (no sharpening), >0 sharpening
    gammaICI=1;                     % ICI Gamma threshold
    ndir=8;                         % number of directions
    niter=6;                        % number of iterations (3 is generally enough)
    restrict=0;                     % restricts number of iterations if noise level is small
else
    %%% HIGH QUALITY SETTINGS
    h1=[1 2 3 5 7 10];              % LPA Kernels' length
    sharparams=[-0.8 -0.6 -0.3];    % -1 zero order 0 first order (no sharpening), >0 sharpening
    gammaICI=0.8;                   % ICI Gamma threshold
    ndir=16;                        % number of directions
    niter=3;                        % number of iterations (3 is generally enough)
    restrict=1;                     % restricts number of iterations if noise level is small
end


fusing=1;                       % fusing type   (1 classical fusing, 2 piecewise constant)
iterrule=1;                     % sigmaiter rule (should follow the fusing type)
itercoef=2/3;                   % compensating factor for recursive standard deviation (def. 2/3)

addnoise=1;                     % add noise to observation
sigma_noise=0.1;                % standard deviation of the noise   (e.g. 25/255)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%--------------------------------------------------------------------------
% LPA KERNELS SIZES
%--------------------------------------------------------------------------
h2=ones(size(h1));
lenh=length(h1);

%--------------------------------------------------------------------------
% WINDOWS PARAMETERS
%--------------------------------------------------------------------------
window_type=1;      % window=1 for uniform,
TYPE=10;            % TYPE IS A SYMMETRY OF THE WINDOW
%                   % 00 SYMMETRIC
%                   % 10 NONSYMMETRIC ON X1 and SYMMETRIC ON X2
%                   % 11 NONSYMMETRIC ON X1,X2  (Quadrants)

disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Recursive Anisotropic LPA-ICI Denoising Demo ')
disp('-----------------------------------------------------------------------------------')
disp(' ')


%--------------------------------------------------------------------------
% MODELLING
%--------------------------------------------------------------------------

y=im2double(imread('image_Cameraman256.png'));
% y=im2double(imread('image_Lena512.png'));
% y=im2double(imread('image_Boats512.png'));

tic
init=0; %2055615866;
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


disp('Initial (noisy) observation criteria values:')
function_Errors(y,z,z,2);

figure
imshow(z), title('Noisy observation')

sigma=function_stdEst2D(z);

if restrict==1
niter=max(1,min(round(sigma*43),niter));     % restrict number of iterations when noise is not too large
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear yh h_opt_Q YICI_Q var_opt_Q

%---------------------------------------------------------
% Kernels construction
%---------------------------------------------------------

% calling kernel creation function
[kernels, kernels_higher_order]=function_CreateLPAKernels([0 0],h1,h2,TYPE,window_type,ndir,ones(2,lenh),1);
[kernelsb, kernels_higher_orderb]=function_CreateLPAKernels([1 0],h1,h2,TYPE,window_type,ndir,ones(2,lenh),1);
disp(' ')
disp(' LPA kernels created');

z_iter=z;
sigmaiter=repmat(sigma,size_z_1,size_z_2);
sigmaiters{1}=sigmaiter;
toc
disp('       ')
disp(['Running iteration ',num2str(1),' ...']);
figure
figure_number=gcf;
imshow(z), title('Noisy observation')
drawnow
htbar=timebar(['Recursive LPA-ICI running   (#iter=',num2str(niter),')'],'Progress');
for momo=1:niter
    sigmaiters{momo}=sigmaiter;
    sharparam=sharparams(min(momo,numel(sharparams)));
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
            yh(:,:,s2)= conv2(z_iter+10000,gh,'same')-10000; % Estimation
            if momo==1
                stdh(:,:,s2)=repmat(sigma*(sum(gh(:).^2))^0.5,size_z_1,size_z_2);  % Std of the estimate
            else
                stdh(:,:,s2)=(conv2(sigmaiter.^2,gh.^2,'same')).^0.5;  % Std of the estimate
            end

        end %% for s2, window sizes
        [YICI,h_opt,std_opt]=function_ICI(yh,stdh,gammaICI,2*(s1-1)*pi/ndir);    %%%% ICI %%%%%
        aaa=reshape(ghorigin(s1,h_opt),size(h_opt));  %origin weight for optimal kernels
        YICI_Q(:,:,s1)=YICI;   h_opt_Q(:,:,s1)=h_opt;   var_opt_Q(:,:,s1)=(std_opt.^2+eps);
        YICI_Final1=YICI_Final1+YICI_Q(:,:,s1)./var_opt_Q(:,:,s1);            %% FUSING %%%%%
        YICI_Final2=YICI_Final2+YICI_Q(:,:,s1)./var_opt_Q(:,:,s1)-z_iter.*aaa./var_opt_Q(:,:,s1);            %% FUSING 2 %%%%%
        var_inv=var_inv+1./var_opt_Q(:,:,s1);
        CWW=CWW+aaa./var_opt_Q(:,:,s1);

        CWW2=CWW2+(aaa./var_opt_Q(:,:,s1)).^2;
        if s1==1
            h_opts1{momo}=h_opt;
        end
        if s1==2
            h_opts2{momo}=h_opt;
        end
        timebar(htbar,(s1+(momo-1)*ndir)/niter/ndir);
    end     %for s1 directions
    YICI_Final1=YICI_Final1./var_inv;
    YICI_Final2=(YICI_Final2+z_iter.*CWW/ndir)./(var_inv-CWW+CWW./ndir);           %% FUSING 2 %%%%%

    if fusing==1, y_hat=YICI_Final1; end
    if fusing==2, y_hat=YICI_Final2; end

    %%%%%%%%%%%%%   END OF ALGORITHM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    y_hats{momo}=y_hat;
    [Err_y_hat,Err_labels]=function_Errors(y,y_hat,z);  %% computes error criteria

    %% TABLE DISPLAY FOLLOWS
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

    %% STD UPDATE FOLLOWS
    sigmaiterold=sigmaiter;
    if iterrule==1;
        sigmaiter=sqrt( (1./(var_inv.^2)).*(  var_inv -  (sigmaiter.^2).*CWW2 + (sigmaiter.*CWW).^2));    % standard fusing
    end
    if iterrule==2;
        sigmaiter=abs((1./(var_inv-CWW+CWW./ndir)).^0.5);      % piecewise constant
    end
    sigmaiter=(sigmaiter*itercoef);
    sigmaiters{momo+1}=sigmaiter;

    %% OBSERVATION UPDATE
    z_iter=y_hat;

    if gcf~=figure_number, figure(figure_number), end
    imshow(y_hat), title(['LPA-ICI estimate (iter. ',num2str(momo),')']);
end   %%% end of recursive LPA-ICI
close(htbar);
close(figure_number);

disp(' ');
toc
figure
imshow(y_hat), title('Recursive LPA-ICI estimate');
version -release; % get matlab release
matlab_R=str2num(ans);

vscale=30;
screensize = get(0,'screensize');       % User's screen size [1 1 width height]
size_patch=min(round(min(screensize([3,3]),screensize([3,3]))/30),[size_z_1,size_z_2]);
if niter>=3
    for figure_counter=1:2
        if figure_counter==1
            range1=[40:110];range2=[110:180];
        else
            range1=(1+floor((size_z_1-1-size_patch(1))*rand))+[1:size_patch(1)];range2=(1+floor((size_z_2-1-size_patch(2))*rand))+[1:size_patch(2)];
        end
        if screensize(3)>=1280|figure_counter==2
            figure
            subplot(2,3,1);
        else
            figure
            subplot(1,2,1);
        end
        surf(flipud(vscale*y(range1,range2))),  axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
        caxis([0 vscale]);
        colormap(gray);
        camproj('perspective');
        axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,
        title('original');

        if screensize(3)>=1280|figure_counter==2
            subplot(2,3,2);
        else
            subplot(1,2,2);
        end
        surf(flipud(vscale*z(range1,range2))),   axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
        camproj('perspective');
        caxis([0 vscale]);         axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,
        title('noisy');

        if screensize(3)>=1280|figure_counter==2
            subplot(2,3,3);
        else
            colormap((gray).^0.3)  %% use gamma correction for better visualizations of surfaces
            axes('position',[0.5 0.93 0.0001 0.0001]); axis off, title('Detailed view of a fragment of the image')
            figure
            subplot(1,2,1);
        end
        surf(flipud(vscale*y_hats{1}(range1,range2))),  axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
        camproj('perspective');
        caxis([0 vscale]);         axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,
        title('LPA-ICI estimate (1^s^t iter.)');

        if screensize(3)>=1280|figure_counter==2
            subplot(2,3,4);
        else
            subplot(1,2,2);
        end
        surf(flipud(vscale*y_hats{2}(range1,range2))),   axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
        camproj('perspective');
        caxis([0 vscale]);
        axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,
        title('LPA-ICI estimate, 2^n^d iter.');

        if screensize(3)>=1280|figure_counter==2
            subplot(2,3,5);
        else
            colormap((gray).^0.3)  %% use gamma correction for better visualizations of surfaces
            axes('position',[0.5 0.93 0.0001 0.0001]); axis off, title('Detailed view of a fragment of the image')
            figure
            subplot(1,2,1);
        end
        surf(flipud(vscale*y_hats{3}(range1,range2))), axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
        camproj('perspective');
        caxis([0 vscale]);
        axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,        title('LPA-ICI estimate, 3^r^d iter.');
        if niter>=4
            if screensize(3)>=1280|figure_counter==2
                subplot(2,3,6);
            else
                subplot(1,2,2);
            end
            surf(flipud(vscale*y_hats{4}(range1,range2))), axis equal, set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[]);set(gca,'ztick',[0 0.5*vscale vscale],'zticklabel',[0 0.5 1]); axis([[1 length(range2)],[1 length(range1)],[0 vscale]]); view([-30,35]); box on;
            camproj('perspective');
            caxis([0 vscale]);
            axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.91 0.0001 0.0001]);  axis off,        title('LPA-ICI estimate, 4^t^h iter.');
        end
        colormap((gray).^0.3)  %% use gamma correction for better visualizations of surfaces
        axes('position',[0.5 0.93 0.0001 0.0001]); axis off, title('Detailed view of a fragment of the image')
    end
end

if niter>=3
    figure
    imshow([[z,y_hats{1},y_hats{2},y_hats{3}];([sigmaiters{1},sigmaiters{2},sigmaiters{3},sigmaiters{4}]).^.5/max(sigmaiters{1}(:))^.5]);
    title('Recursive estimates (starting from the initial noisy observation) and their corresponding estimated standard deviations')
    if niter>=4
        figure
        imshow([[h_opts1{1},h_opts1{2},h_opts1{3},h_opts1{4}];[h_opts2{1},h_opts2{2},h_opts2{3},h_opts2{4}]],[])
        title(sprintf('ICI adaptive directional scales (first four iterations, first two directions)\nWhite denotes large scales (strong smoothing), black small scales (no filtering)'))
    else
        imshow([[h_opts1{1},h_opts1{2},h_opts1{3}];[h_opts2{1},h_opts2{2},h_opts2{3}]],[])
        title(sprintf('ICI adaptive directional scales (first three iterations, first two directions)\nWhite denotes large scales (strong smoothing), black small scales (no filtering)'))
    end
end
