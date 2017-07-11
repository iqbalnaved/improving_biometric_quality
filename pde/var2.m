function v=var2(x)

[M,N]=size(x);
m=mean2(x);
e=x-m;
v=sum(sum(e.*e))/(M*N);