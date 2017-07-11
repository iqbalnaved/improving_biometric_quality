% pm.m - Anisotropic Diffusion 

function ZN = pm(ZN,iterate);

[m,n] = size(ZN);

lambda = 0.250;
%lambda = .125;
th=.0001;
ZN_save=ZN;
%e_save=100;
for i = 1:iterate, 
   %i 
   %ZN=ZN+lambda*div_grad_vese(ZN, 1, 0.0001); 
   ZN=ZN+lambda*div_grad_cs(ZN, 1, 0.0001); 
%     e(i)=var2(ZN-Z);
%     if i>10&e(i)>max(e(i-3:i-1))
%         break;
%     end
%    else
%    ZN_save=ZN;
%    %e_save=e(i);
% end
end