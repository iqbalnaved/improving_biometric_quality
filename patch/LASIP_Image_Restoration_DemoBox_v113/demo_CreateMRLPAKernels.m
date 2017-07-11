%  Multiresolution (MR) LPA Kernels demo  (demo_Create_MR_LPAKernels.m)
%
%  Alessandro Foi, Vladimir Katkovnik - Tampere University of Technology - September 2005
% ------------------------------------------------------------------------------------------ 
%
%  Creates and draws multiresolution (MR) LPA two-dimensional kernels.
%  Kernels can be smoothing and differentiating.
%  Window functions allow to obtain symmetric, nonsymmetric and sectorial support kernels.
% 
% ------------------------------------------------------------------------------------------
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%--------------------------------------------------------------------------
% LPA ORDER AND KERNELS SIZES
%--------------------------------------------------------------------------
m=[1,1];            % Order of LPA
h1=[2 4 6 9 13 17]; % sizes of the kernel; for differentiation the smallest scale h(1) should be large enough
%                   % (depending on the order of the derivative)
h2=h1;              %   row vectors h1 and h2 need to have the same length

DerivativeOrder=1;        % Derivative order of the differentiating kernels
DerivativeOrder_imaged=1; % Derivative order of the differentiating kernels shown in images
%                         % 1 for the smoothing kernel
%                         % 2 for the first derivative and so on
%                         % available DerivativeOrder can be seen from the order map given by
%                         % kernels_higher_order{direction,size,2}

%--------------------------------------------------------------------------
% WINDOWS PARAMETERS
%--------------------------------------------------------------------------
sig_winds=[ones(size(h1))*0.3 ; ones(size(h2))*0.3];    % Gaussian parameter
beta=1;                     % Parameter of window 6

window_type=2 ;  % window_type=1 for uniform,
% window_type=2 for Gaussian
% window_type=6 for exponentions with beta
% window_type=8 for Interpolation

TYPE=00;         % TYPE IS A SYMMETRY OF THE WINDOW
% 00 SYMMETRIC
% 10 NONSYMMETRIC ON X1 and SYMMETRIC ON X2
% 11 NONSYMMETRIC ON X1,X2  (Quadrants)

%--------------------------------------------------------------------------
% DIRECTIONAL PARAMETERS
%--------------------------------------------------------------------------
ndir=1;       % number of directions



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lenh=length(h1);
[kernels, kernels_higher_order]=function_CreateLPAKernels(m,h1,h2,TYPE,window_type,ndir,sig_winds,beta);
if DerivativeOrder>1
    'Table of the derivative orders'
    kernels_higher_order{1,1,2}     % gives Table of the derivative orders
    %%% In order to select the derivative to be shown select the number
    %%% of the line of the table
    %%% DerivativeOrder_imaged=1
    'Derivative order selected'
    kernels_higher_order{1,1,2}(DerivativeOrder_imaged,1:2)
end

for s1=1:ndir     % directional index
    for s2=1:lenh     % kernel size index
        gh=kernels_higher_order{s1,lenh-s2+1,1}(:,:,DerivativeOrder_imaged);   %gets single kernel from the cell array
        if s2==1, clear GH
            [maxh11,maxh12]=size(gh); GH(1:maxh11,1:maxh11,1:lenh)=0;
        end
        [N1gh, N2gh]=size(gh);
        LL= [(maxh11+1)/2-(N1gh-1)/2:(maxh11+1)/2+(N1gh-1)/2] ;
        GH(LL,LL,s2)=gh;
    end %% for s2, window sizes
    
    DELTA_GH=GH(:,:,1);
    for jj=2:lenh % Differences of the kernels
        delta_gh=GH(:,:,jj)-GH(:,:,jj-1);

        DELTA_GH(:,:,jj)=delta_gh;
    end
[GhOrth, D]=qr(reshape(DELTA_GH,[size(DELTA_GH,1)*size(DELTA_GH,2),size(DELTA_GH,3)]),0); % Orthonormalization of the differences
    GhOrth=reshape(-GhOrth,[size(DELTA_GH,1),size(DELTA_GH,2),size(DELTA_GH,3)]);
    figure
    max_absGhOrth=max(abs(GhOrth(:)));
    for iii=1:min(lenh,6) 
        subplot(2,3,iii), surf(GhOrth(:,:,iii));   % shading interp
        max_absGhOrth=max(max(abs(GhOrth(:,:,iii))));
          caxis([-max_absGhOrth max_absGhOrth]);
        camproj('perspective'); axis tight, box on
        axis_pos=get(gca,'position'); axes('position',[axis_pos(1)+axis_pos(3)/2 axis_pos(2)+axis_pos(4)*0.95 0.0001 0.0001]);  axis off,        
        title(['scale ', num2str(iii)]),
    end
    figure
for iii=1:min(lenh,6)
    max_absGhOrth=max(max(abs(GhOrth(:,:,iii))));
     subplot(2,3,iii), imshow(GhOrth(:,:,iii),[]);   caxis([-max_absGhOrth max_absGhOrth]); title(['scale ', num2str(iii)]), colorbar('EastOutside'),
end
end

