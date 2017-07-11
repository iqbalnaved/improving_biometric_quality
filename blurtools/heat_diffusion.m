function y=heat_diffusion(x,iter)

% Laplacian operation
h=[0 -1 0;-1 4 -1;0 -1 0]/4;
lambda=0.1;
y=x;
for i=1:iter
dx=imfilter(y,h,'symmetric');
y=y-lambda*dx;
end