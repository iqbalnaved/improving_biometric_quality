% Multiresolution thresholding (function_MR_Thresholding)
% [y, index]=MR_thresholding(x,t,thresholding_type)
%
% INPUTS: 
% x  input matrix; 
% t  threshold parameter,
%         thresholding types: thresholding_type=1  for hard-thresholding;
%                             thresholding_type=2  for soft-thresholding;
%                             thresholding_type=3  for Stein's thresholding;
%                             thresholding_type=31 for smoothed Stein's thresholding.
%
% OUTPUTS:
% y     is a matrix of the thresholded MR spectrums 
% index is a matrix of indices: 1 indicate that the element of the matrix is included in the estimates; 
%                               0 means that the elements is killed
% Last Corrections: V. Katkovnik, October 28, 2002, Gwangju, South Korea


function [y, index]=MR_thresholding(x,t,thresholding_type)
global sigma q SmoothedSteinThreshold


if thresholding_type==1 % hard-thresholding
    
index=(abs(x)>t);

    y=x.*index;
end



if thresholding_type==2  % soft-thresholding
    
    index=(abs(x)>t);

y=(abs(x)-t).* index.*sign(x);
end

if thresholding_type==21 %  Abramovich rule for the threshold
 
   
  xx=(x(:)); nn=[1:length(xx)]; 

  [x_sort,i_sort]=sort(abs(xx));
   
   x_sort1=fliplr(x_sort');% i_sort1=fliplr(i_sort)'; 
    
   prob1=2*(1-normcdf(x_sort1,0,t));
      
 
   x_thresh=(prob1<=nn*q/length(xx));
     if  x_thresh==0
  t_factor=x_sort(1);
     else
 
    tt=max(find(x_thresh==1));
   
    t_factor=x_sort1(tt);
end
index=(abs(x)>=t_factor);
 y=x.*index;
    
 end

if thresholding_type==3  % Stein's thresholding
   
C=1-(t^2)./((abs(x)+0.0001).^2); index=(C>0);
y=x.*C.* index;
end
if thresholding_type==31  % Smoothed Stein Thresholding
% includes linear and median filtering defined by matrix B

B=ones(3,3); % line-wise filetring is used with the line-wise kernels
B=B/sum(B(:));
[X1]=filter2(B,x.^2,'same');
C=1-(t^2)./(X1+0.0001); index=(C>0);
y=x.*C.* index;
end


if thresholding_type==4  % eneralized Stein's thresholding
   alpha=2;

 B=ones(3,3); % line-wise filetring is used with the line-wise kernels
 B=B/sum(B(:));
 [X1]=filter2(B,x.^2,'same');

 C=1-(t^alpha)./(X1+0.0001); index=(C>0);
 y=x.*C.* index;
end

