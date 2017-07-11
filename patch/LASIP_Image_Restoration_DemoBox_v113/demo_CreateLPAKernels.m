% Creates LPA kernels cell array demo  (demo_CreateLPAKernels)
%
% Alessandro Foi - Tampere University of Technology - 2003-2005
% -----------------------------------------------------------------------
%
%  Builds kernels cell arrays kernels{direction,size}
%                  and        kernels_higher_order{direction,size,1:2}
%               kernels_higher_order{direction,size,1}  is the 3D matrix
%                   of all kernels for that particular direction/size
%               kernels_higher_order{direction,size,2}  is the 2D matrix
%                   containing the orders indices for the kernels
%                   contained in kernels_higher_order{direction,size,1}
%
%   ---------------------------------------------------------------------
% 
%   kernels_higher_order{direction,size,1}(:,:,1) is the funcion estimate kernel
%   kernels_higher_order{direction,size,1}(:,:,2) is a first derivative estimate kernel
%
%   kernels_higher_order{direction,size,1}(:,:,n) is a higher order derivative estimate kernel
%   whose orders with respect to x and y are specified in
%   kernels_higher_order{direction,size,2}(n,:)=
%                           =[xorder yorder xorder+yorder]
%   
%

%--------------------------------------------------------------------------
% LPA ORDER AND KERNELS SIZES
%--------------------------------------------------------------------------
    m=[2,1];        % THE VECTOR ORDER OF LPA;
    
    h1=[1 5 15];    %   sizes of the kernel
    h2=[1 3 6];    %   row vectors h1 and h2 need to have the same lenght


%--------------------------------------------------------------------------
% WINDOWS PARAMETERS
%--------------------------------------------------------------------------
    sig_winds=[ones(size(h1))*0.3 ; ones(size(h2))*0.3];    % Gaussian parameter
    beta=1;                     % Parameter of window 6

    window_type=2 ;  % window_type=1 for uniform, window_type=2 for Gaussian
                     % window_type=6 for exponentions with beta
                     % window_type=8 for Interpolation

    TYPE=00;         % TYPE IS A SYMMETRY OF THE WINDOW
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
% DIRECTIONAL PARAMETERS
%--------------------------------------------------------------------------
    directional_resolution=4;       % number of directions

    
                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(' Demo: Creates LPA kernels ')
disp('-----------------------------------------------------------------------------------')
disp('  ')
disp('Creating LPA kernels...  ')
tic
[kernels, kernels_higher_order]=function_CreateLPAKernels(m,h1,h2,TYPE,window_type,directional_resolution,sig_winds,beta);
disp(sprintf('\b\b Kernels created.'))
toc
disp(' ')
disp('Use the command   utility_DrawLPAKernels   to draw figures of the created LPA kernels.');
disp(' ')