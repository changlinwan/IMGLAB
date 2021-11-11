function xx = contour_scale_factor(curv, P, T, LE, LE2)

ln = max(LE(:));

curv_u = zeros(ln,1);
curv_t = zeros(ln,1);
numv_b = zeros(ln,1);
numv_c = zeros(ln,1);
numv_n = ones(ln,1);

idx = find(LE>0);
for i=1:ln
    ix = LE(idx)==i;
    H = idx(ix);
    idx(ix) = [];
    
    num = numel(H);
    if num == 0
        continue;
    end
%     ix = LE2(H)==0;
    
    tv = T(H);
    cv = curv(H);
    
    ix1 = tv > mean(tv);
    ix2 = cv > mean(cv);
    
    curv_t(i) = .5*(mean(tv(ix1))+mean(tv)+std(tv(~ix1)));
    curv_u(i) = .5*(mean(cv(ix2))+mean(cv)+std(cv(~ix2)));
    
    numv_b(i) = num;
    numv_n(i) = sum(curv(H) < .2);
    numv_c(i) = numv_n(i)/num;
%     UN = unique(sort(LE2(H)));
%     numv_m(i) = numel(UN);
end

numv_n = numv_n/mean(numv_n);
numv_c = numv_c/mean(numv_c);

% ixx = numv_m < 5;
% 
% numv_m = 4*numv_b./numv_m;
% numv_m(ixx) = numv_b(ixx);
% 
% numv_m = numv_m/mean(numv_m);

xx = [curv_u,curv_t,numv_b,numv_c,numv_n];

end

