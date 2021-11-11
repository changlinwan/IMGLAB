function [radius, P, V, T] = turnradiusblock2(aSmooth, tx, ty, u, va, iv, iEDGE2, bin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Changlin Wan
% E-mail: wancl@21cn.com
% Function: turnradiusblocks2
% Version: 1.0
% Date: 2021/11/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(tx);

nc = [bin bin 0 -bin -bin -bin 0 bin];
nr = [0 -bin -bin -bin 0 bin bin bin];

% nc = [bin*ones(1,bin), bin:-1:-bin,-bin*ones(1,2*bin-1),-bin:bin,bin*ones(1,bin-1)];
% nr = [0:-1:-bin, -bin*ones(1,2*bin-1),-bin:bin, bin*ones(1,2*bin-1),bin:-1:1];

nn = numel(nc);

a = aSmooth(:);

[rx,cx] = ind2sub([m,n],iEDGE2);
ix = rx > 2*bin & cx > 2*bin & rx <= m-2*bin & cx <= n-2*bin;
rx = rx(ix);
cx = cx(ix);
xi = iEDGE2(ix);

nx = numel(xi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 最8方向灰度改变值从小到大排序

mx = zeros(m*n,nn);
mx1 = zeros(m*n, nn);
for i=1:nn
    mx1(:,i) = i;
end
[x,y] = sort(u(iEDGE2,:),2);
mx(iEDGE2,:) = x;
mx1(iEDGE2,:) = y;

my = sort(8*va.*va,2);

m8 = mx(xi,end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 计算8邻域中灰度变化最小的两个方向的cosine夹角

u2 = [tx(:),ty(:)];

kk = 2;

cosn = ones(m*n,kk);
miw = zeros(nx,1);
mjw = zeros(nx,kk);

A = u2(xi,:);
U2 = sum(u2.*u2,2) + eps;
A2 = U2(xi);

dv = zeros(nx,kk);
% dt = zeros(nx,kk);


for k=1:kk
    rr = nr(mx1(xi,k))';
    cc = nc(mx1(xi,k))';
    i = rx+rr; j=cx+cc;
    ix = sub2ind([m,n],i,j);

    i = rx+2*rr; j=cx+2*cc;
    ix2 = sub2ind([m,n],i,j);

    cosn(xi,k) = sum(u2(ix,:).*A,2)./sqrt(U2(ix).*A2);
    
    dv(:,k) = (a(xi)-2*a(ix)+a(ix2)).^2;
    
    mjw(:,k) = m8-mx(xi,k);
    
    miw = miw + mjw(:,k);
end

beta = zeros(m*n,kk);
beta(xi,:) = bsxfun(@rdivide, mjw, miw+eps);

dvv = sum(beta(xi,:).*dv,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 计算角点特征值

V = zeros(m,n);
T = zeros(m,n);

m2 = mx(:,1)+mx(:,2);
n2 = my(:,1)+my(:,2);

V(xi) = sqrt(m2(xi).*n2(xi));
V(xi) = sqrt(V(xi).*dvv)+V(xi);

sinn = sqrt(1-cosn.*cosn);
sinn2 = sqrt(0.5*(1-cosn));

T(:) = sum(beta.*(sinn+sinn2),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 计算8邻域中最大相似区域
CHR = local_connection(va(xi,:),xi,[m,n]);

P = (1. + CHR);

radius = (.5+P).*V.*T;
% radius = T;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 角点极值
% 
% rr1 = nr(mx1(xi,1));
% cc1 = nc(mx1(xi,1));
% rr2 = nr(mx1(xi,2));
% cc2 = nc(mx1(xi,2));
% 
% ixx = abs(rr1+rr2) > bin/2 | abs(cc1+cc2) > bin/2;
% 
% LE2 = zeros(m,n);
% 
% LE2(sub2ind([m,n],rx(ixx),cx(ixx))) = 1;
% 
