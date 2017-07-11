function y=unsharp_masking(x,lambda)

% Laplacian operation
h=[0 -1 0;-1 4 -1;0 -1 0]/4;
dx=filter2(h,x);
y=x+lambda*dx;