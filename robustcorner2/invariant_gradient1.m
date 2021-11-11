function [dx, dy, aSmooth, u, va, xi] = invariant_gradient1(aSmooth, iEDGE, bin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Changlin Wan
% E-mail: wancl@21cn.com
% Function: invarient_gradient
% Version: 1.0
% Date: 2021/11/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(aSmooth);

if bin > 1
    h = fspecial('average',bin);
    aSmooth = imfilter(aSmooth,h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 不同尺度下的8邻域

nc = [bin bin 0 -bin -bin -bin 0 bin];
nr = [0 -bin -bin -bin 0 bin bin bin];

nn = numel(nc);

a2 = zeros(m*n,numel(nc));

[rx,cx] = ind2sub([m,n],iEDGE);
ix = rx > bin & cx > bin & rx <= m-bin & cx <= n-bin;
rx = rx(ix);
cx = cx(ix);
iEDGE = iEDGE(ix);
for i=1:nn
    ix = sub2ind([m,n],rx+nr(i),cx+nc(i));
    iEDGE = [iEDGE;ix];
end

x = zeros(m,n);
x(iEDGE) = 1;
[rx,cx] = find(x(bin+1:m-bin,bin+1:n-bin));
rx = rx+bin;
cx = cx+bin;
xi = sub2ind([m,n],rx,cx);

nx = numel(xi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 不变梯度需要的计算的像素点位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8邻域灰度值
nix = zeros(nx,nn);
for i=1:nn
    nix(:,i) = sub2ind([m,n],rx+nr(i),cx+nc(i));
end
a2(xi,:) = aSmooth(nix);
a = aSmooth(:);
va = bsxfun(@minus, a2, a);

u = zeros(m*n,nn);
A2 = a2(xi,:);
A = a(xi);
for i=1:nn
    ix = sub2ind([m,n],rx+nr(i),cx+nc(i));   
    u(xi,i) = sum((a2(ix,:)-A2).^2,2)+(a(ix)-A).^2;
end

dx = [];
dy = [];

end

