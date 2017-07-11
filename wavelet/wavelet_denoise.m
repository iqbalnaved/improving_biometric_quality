function x_hat=wavelet_denoise(y,wname,th)

[Y,S]=wavedec2(y,4,wname);
%th=100;
Y=Y.*(abs(Y)>th);
x_hat=waverec2(Y,S,wname);
