% Multiresolution Filtering, based on the MR analysis, thresholding and synthesis. (function_MR_filtering)
% see Katkovnik, “Multiresolution nonparametric regression: a new approach to pointwise spatial adaptation,”
% Digital Signal Processing vol. 15, pp. 73–116, 2005.
%
% [Y_out,Var_out]=function_MR_Filtering(g,thresholding_type, threshold_parameter);
%
% INPUTS:
% g is a 3D matrix of the kernels obtained for all set of scales, size N1*N2*N_h
% thresholding_type indicated a type of the thresholding used in filtering
%
%
% OUTPUTS:
% Y_out is a matrix (N1*N2) of the function estimates
% Var_out is a matrix (N1*N2) of the variances of the estimates
%---------------------------------------------------------------


function [Y_out,Var_out]=function_MR_Filtering(gh,thresholding_type, threshold_parameter);

global z y lenh sigma h1

%-------------------------------------------------------------------------------------
% Initialization
t=ones(1,lenh)*threshold_parameter; % vector of threshold parameters (allows to have different threshold for different subbands)
[xN,yN]=size(y); % image size

%--------------------------------------
% ANALYSIS
%--------------------------------------
% Kernels in gh are ordered starting from the largest h
DELTA_GH=gh(:,:,1);
for jj=2:lenh % Differences of the kernels
    DELTA_GH(:,:,jj)=gh(:,:,jj)-gh(:,:,jj-1);
end

[GhOrth, D]=qr(reshape(DELTA_GH,[size(DELTA_GH,1)*size(DELTA_GH,2),size(DELTA_GH,3)]),0); % Orthonormalization of the differences
GhOrth=reshape(sign(D(1,1))*GhOrth,[size(DELTA_GH,1),size(DELTA_GH,2),size(DELTA_GH,3)]);
q_coefficient=sign(D(1,1))*sum(D,2);     % coefficients used in the synthesis (Proposition 17, page 394)
% for j=1:lenh % Estimates and variances
%     Y_spec(1:xN,1:yN,j)=conv2(z,GhOrth(:,:,j),'same'); % spectral MR estimates
% end % END of j
%---------------------------------------------------
% SYNTHESIS
%-------------------------------------------------------
t(1)=0; % do not threshold approximation coefficients
Y_out=zeros(size(y)); Var_out=Y_out;
for s_spec=1:lenh
    %    [AA,BB]=function_MR_Thresholding(Y_spec(1:xN,1:yN,s_spec),sigma*t(s_spec),thresholding_type); % Universal Stochastic Thresholding
    [AA,BB]=function_MR_Thresholding(conv2(z,GhOrth(:,:,s_spec),'same'),sigma*t(s_spec),thresholding_type); % Universal Stochastic Thresholding
    if s_spec==1
        figure,  figure_number2=gcf;
    else
        if gcf~=figure_number2,  figure(figure_number2),        end
    end
    if s_spec<=5
        subplot(2,3,s_spec)
        imshow(AA,[]), title(['Subband ', num2str(s_spec),' after thresholding']),
    end
    Y_out=Y_out+q_coefficient(s_spec)*AA; % synthesis (estimate calculation)
    Var_out=Var_out+BB*q_coefficient(s_spec)^2;% variances of the estimates
end
subplot(2,3,6)
imshow(Y_out), title(['Estimate reconstructed from ', num2str(lenh),' subbands']),
