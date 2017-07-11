function r = shiftImage( f, di, dj )
% Translation of an image
[N,M]=size(f);

if di>0
    iind = [(di+1):N, N*ones(1,di)];
elseif di<0
    iind = [ ones(1,-di), 1:(N+di) ];
else
    iind = 1:N;
end
if dj>0
    jind = [(dj+1):M, M*ones(1,dj)];
elseif dj<0
    jind = [ ones(1,-dj), 1:(M+dj) ];
else
    jind = 1:M;
end
r = f( iind, jind );
