% Generates error-diffusion halftoned image.
%
% Alessandro Foi - Tampere University of Technology - 2003-2005
% -------------------------------------------------------------
%
%   z = function_ErrorDiffusion(y,kernel)
%
% INPUTS:
%    y       : original (continuous-tone) image
%    kernel  : error-kernel to use  (OPTIONAL)
%
% OUTPUT:
%    z       : halftone image
%
%
% !! Original image y must have range [0,1]
%
%
% The supplied error-kernel must be causal and its width must be odd.
% The center of the kernel is middle element of the top row.
%
% Two well-known error-diffusion kernels are already implemented:
%
% kernel=1  uses the Floyd-Steinberg kernel [0 0 7;3 5 1]/16
% kernel=2  uses the Jarvis et al. kernel [0 0 0 7 5;3 5 7 5 3;1 3 5 3 1]/48
%
% If kernel is not specified the Floyd-Steinberg kernel is used as default.
%
%
% References:
% R.W. Floyd and L.Steinberg, "An adaptive algoritm for spatial grayscale",
% Proc. SID, vol.17, no.2, pp.75-77, 1976.
% J.Jarvis, C.Judice, and W.Ninke, "A survey of techniques for the
% display of continuous tone pictures on bilevel displays",
% Comp. Graph and Image  Proc., vol.5, pp.13-40, 1976.
%


function  z=function_ErrorDiffusion(y,kernel)

if nargin==1
    kernel=1;  %% defaults to Floyd-Steinberg kernel
end
if numel(kernel)==1&kernel==1
    error_kernel=[0 0 7;3 5 1];  %% Floyd-Steinberg
elseif numel(kernel)==1&kernel==2
    error_kernel=[0 0 0 7 5;3 5 7 5 3;1 3 5 3 1];  %% Jarvis et al.
elseif nargin>1&numel(kernel>1)
    error_kernel=kernel;
end
[ek1,ek2]=size(error_kernel);
ek1=ek1-1; ek2=(ek2-1)/2;
if round(ek2)~=ek2
    disp(' ')
    disp('The width of the error-kernel must be odd.  Center of the kernel is middle element of the top row.');
    disp(' ')
    return
else
    if (sum(abs(error_kernel(1,1:ek2+1)))>0)
        disp(' ')
        disp('Error-kernel is not causal!');
        disp('Error-kernel must be causal and its width must be odd.  Center of the kernel is middle element of the top row.');
        disp(' ')
        return
    end
end

error_kernel=error_kernel/sum(error_kernel(:));  %% normalizes the error kernel
[size_y_1,size_y_2]=size(y);
z=zeros(size_y_1+ek1,size_y_2+ek2+ek2);
z(1:end-ek1,1+ek2:end-ek2)=y;  %% initial condition is original image
y=z;  %% pads original image borders
for pos1 = 1:size_y_1
    for pos2 = 1+ek2:size_y_2+ek2
        z(pos1,pos2)=(y(pos1,pos2)>=0.5);  %% binarization
        e = -z(pos1,pos2) + y(pos1,pos2);  %% error
        y(pos1:pos1+ek1,pos2-ek2:pos2+ek2)=error_kernel*e+y(pos1:pos1+ek1,pos2-ek2:pos2+ek2);  %% ERROR DIFFUSION
    end
end
z=z(1:size_y_1,1+ek2:size_y_2+ek2);  %% removes borders
