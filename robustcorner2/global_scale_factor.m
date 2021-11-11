function xx = global_scale_factor(curv, P, T, LE, LE2, CC2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Changlin Wan
% E-mail: wancl@21cn.com
% Function: global_scale_factor
% Version: 1.0
% Date: 2021/11/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LNS = unique(sort(LE(CC2>0)));
LNS(LNS==0) = [];

% LNS = 1:max(LE(:));

numv_n = zeros(max(LE(:)),1);
curv_u = numv_n;
curv_t = numv_n;
numv_b = numv_n;
numv_c = numv_n;

idx = find(LE>0);
for k=1:numel(LNS)
    i = LNS(k);
    ix = LE(idx)==i;
    H = idx(ix);
    idx(ix) = [];
    
    num = numel(H);
    if num == 0
        continue;
    end
    cv = curv(H);
    tv = T(H);

    ix1 = cv > mean(cv);
    ix2 = tv > mean(tv);
    
    curv_u(i) = .5*(mean(cv)+std(cv(~ix1))+mean(cv(ix1)));
    curv_t(i) = .5*(mean(tv)+std(tv(~ix2))+mean(tv(ix2)));
    numv_b(i) = num-8;
    numv_c(i) = sum(cv>30)./numv_b(i);

    UN = unique(sort(LE2(H)));
    numv_n(i) = numel(UN);
end

ixx = numv_n < 5;

numv_nn = 2*(numv_b./(numv_n));
numv_nn(ixx) = (numv_b(ixx));

ma = mean(numv_nn);
numv_nn(numv_nn>2*ma) = 2*ma;
numv_nn = numv_nn/mean(numv_nn);

numv_c(numv_c==0) = max(numv_c)/100;

xx = [curv_u,curv_t,numv_n,numv_b,numv_c,numv_nn];
xx(xx<0) = 0;

end

