function x_hat=dct_denoise(y,th)

X=blkproc(y,[8 8],@dct2);
Y=X.*(abs(X)>th);
x_hat=blkproc(Y,[8 8],@idct2);
