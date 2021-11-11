function [cm,rn,cn,ix] = corner_curvature12(a, mag2, R2, rr, cc, tx, ty, EDGE, curv, V, P, T, bin, th)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Changlin Wan
% E-mail: wancl@21cn.com
% Function: corner_curvature12
% Version: 1.0
% Date: 2021/11/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(curv);

mag = (tx.^2+ty.^2);
nix = a(:) > 0 | mag(:)>0;

minn = round(16*sqrt(sum(nix(:)))/400);



EDGE2 = double(bwareaopen(EDGE, minn, 8));

LE3 = bwlabel(EDGE2);

curv = curv + rand(size(curv))*1e-8;
T = T + rand(size(curv))*1e-8;

h = fspecial('gaussian',[5 5],1.5);
curv2 = EDGE.*(filter2(h,curv)./(filter2(h,curv>0)+eps));
% tx2 = tx.*EDGE;
% tx2 = EDGE.*(filter2(h,tx2)./(filter2(h,EDGE>0)+eps));
% ty2 = ty.*EDGE;
% ty2 = EDGE.*(filter2(h,ty2)./(filter2(h,EDGE>0)+eps));

if bin < 3
    w = round(3.5*bin);
else
    w = round(4.7*bin);
end
% w = 7*bin;
w2 = 2*w+1;
w3 = w+1;

idx = sub2ind([m,n],rr,cc);

ix = rr > w & cc > w & rr <= m-w & cc <= n-w;
rn = rr(ix);
cn = cc(ix);
idx = idx(ix);

X = ones(numel(rn),12);
X(:,1) = 1:numel(rn);
X(:,2) = rn;
X(:,3) = cn;

CC = zeros(m,n);
CC(idx) = 1;
% CC2 = filter2(ones(7),CC);

jtn2 = P(idx);

[tr,tc] = meshgrid(1:w2,1:w2);
tr = tr(:);
tc = tc(:);

tix = 1:w2*w2;

ix1 = tr == 1 | tr == w2 | tc == 1 | tc == w2;
tix1 = tix(ix1);

ix2 = tr == 5 | tr == w2-4 | tc == 5 | tc == w2-4;
tix2 = tix(ix2);

[tr,tc] = meshgrid(1:w3,1:w3);
step = round(w3/2);

tex1 = sub2ind([w2,w2], tr(:), tc(:));
tex2 = sub2ind([w2,w2], tr(:)+step, tc(:));
tex3 = sub2ind([w2,w2], tr(:)+w, tc(:));
tex4 = sub2ind([w2,w2], tr(:), tc(:)+step);
tex5 = sub2ind([w2,w2], tr(:), tc(:)+w);
tex6 = sub2ind([w2,w2], tr(:)+step, tc(:)+step);
tex7 = sub2ind([w2,w2], tr(:)+w, tc(:)+w);
tex8 = sub2ind([w2,w2], tr(:)+step, tc(:)+w);
tex9 = sub2ind([w2,w2], tr(:)+w, tc(:)+step);

tex = [tex1,tex2,tex3,tex4,tex5,tex7,tex8,tex9];

% figure(15),imshow(segmentation_visualize(LE1));
% figure(12),imshow(segmentation_visualize(LE2));
% figure(13),imshow(segmentation_visualize(LE3));

xx = contour_scale_factor(curv,P,T,LE3, []);

curv_u = xx(:,1);
curv_t = xx(:,2);
numv_n = xx(:,5);

eb3 = false(w2);


mn = mag2(idx);
thr = .1;
if bin == 1 || mean(mn) > 4e-3
    tmh = 0;
else
    tmh = (mean(mag2(nix)))*2e-4/(bin*bin*(mean(mn)));
    for i=1:numel(rn)
        ab = a(rn(i)-w:rn(i)+w,cn(i)-w:cn(i)+w);
        mn(i) = texon_measure(ab,tex,tex6);
    end
    mn = mn.*V(idx)/mean(V(idx));
    if sum(mn<tmh)/numel(mn) > .2
        tmh = tmh*1.5;
    end
end

R3 = R2(idx);
for i=1:numel(rn)
    if (abs(rn(i) - 192) < 5 && abs(cn(i) - 14) < 5)
        disp(i);
    end
    if mn(i) < tmh && R3(i) < thr*.1
        X(i,1) = 0;
        continue;
    end
    
    tb = T(rn(i)-w:rn(i)+w,cn(i)-w:cn(i)+w);
    ub = curv(rn(i)-w:rn(i)+w,cn(i)-w:cn(i)+w);
    eb = EDGE2(rn(i)-w:rn(i)+w,cn(i)-w:cn(i)+w);
    cb = CC(rn(i)-w:rn(i)+w,cn(i)-w:cn(i)+w);
    ub2 = curv2(rn(i)-w:rn(i)+w,cn(i)-w:cn(i)+w);
    
    [r2,c2] = find(eb);
    
    [r2c,c2c] = find(cb);
    
    if numel(r2) == 0
        X(i,1) = 0;
        continue;
    end
    
    if sum(eb(tix1)) < 3 || sum(eb(tix2)) < 3
        L = LE3(rn(i)-w:rn(i)+w,cn(i)-w:cn(i)+w);
        ixn = unique(sort(L(L>0)));
        num = numel(ixn);
    else
        [L,num] = bwlabel(eb);
        ixn = 1:num;
    end
    
    for k=1:num
        ix = find(L==ixn(k));
        if numel(ix)>2
            bb = ub2(ix);
            eb(ix) = sum(bb)/w2+bb;
        else
            eb(ix) = 1e-10;
        end
    end
    eb = eb / mean(eb(eb>0));
    
    [dm,ii] = min((((r2-w3).^2+(c2-w3).^2+4)./(eb(eb>0))));
    
    rr2 = rn(i)+r2(ii)-w3;
    cc2 = cn(i)+c2(ii)-w3;
    
    lns = LE3(rr2,cc2);
    
    if numel(lns) ==0 || dm > 12
        X(i,1) = 0;
        continue;
    end
    
    u = curv_u(lns);
    t = curv_t(lns);
    
    c = numv_n(lns);
    
    LN = L(r2(ii),c2(ii));
    
    eb2 = L==LN;
    
    nb = 1;
    
    n2 = numel(r2c);
    
    if n2 > 1
        ddd = sqrt((w3-r2c).^2+(w3-c2c).^2);
        CC(rn(i),cn(i)) = w3/mean(ddd(ddd>0));
        
        mca = 0;
        mc1 = ub2(r2(ii),c2(ii));
        for j=1:n2
            [~,jj] = min(((r2-r2c(j)).^2+(c2-c2c(j)).^2+4)./(eb(eb>0)));
            
            dd = ddd(j);
            
            if LE3(rn(i)+r2(jj)-w3,cn(i)+c2(jj)-w3) == lns
                if dd > 0 && dd < 5 &&  mc1 < ub2(r2(jj),c2(jj)) && mc1 < .5
                    X(i,1) = 0;
                    continue;
                end
                if dd < 8 && ub2(r2(jj),c2(jj)) > mca
                    mca = ub2(r2(jj),c2(jj));
                end
            end
        end
        if mc1 < mca
            if mc1 < .35
                X(i,1) = 0;
                continue;
            else
                nb = mc1/mca;
            end
        end
    end
    
    if R2(idx(i)) > thr
        continue;
    end
    
    cen = round(w3/3);
    eb3(w3-cen:w3+cen,w3-cen:w3+cen) = eb2(w3-cen:w3+cen,w3-cen:w3+cen);
    
    eeb = eb2;
    eeb(w3-cen:w3+cen,w3-cen:w3+cen) = 0;
    
    ph = 1.3;
    jtn = 0;
    
    dn = sum(eb(:)>0) - sum(eb2(:)>0);
    if (jtn2(i) > ph || dn > 3)
        [L,jtn] = bwlabel(eeb);
        if jtn > 1
            ma = [];
            
            k1 = 0;
            for k=1:jtn
                brn = L == k;
                [r3,c3] = find(brn);
                if numel(r3) > 2
                    ix = sub2ind([m,n],r3+rn(i)-w3,c3+cn(i)-w3);
                    
                    dx = mean(tx(ix));
                    dy = mean(ty(ix));
                    k1 = k1+1;
                    ma(k1,:) = [dx,dy];
                else
                    jtn = jtn-1;
                end
            end
            [r3,c3] = find(eb3);
            k2 = numel(r3);
            cosn2 = [];
            if k2 > 0 && jtn > 0
                ix = sub2ind([m,n],r3+rn(i)-w3,c3+cn(i)-w3);
                bx = repmat(mean(tx(ix)),jtn,1);
                by = repmat(mean(ty(ix)),jtn,1);
                mx = ma(:,1);
                my = ma(:,2);
                cosn2 = (bx.*mx+by.*my)./sqrt((bx.*bx+by.*by).*(mx.*mx+my.*my)+eps);
            end
            
            if jtn > 1
                mx = repmat(ma(:,1),1,k1);
                my = repmat(ma(:,2),1,k1);
                
                cosn = (mx.*mx'+my.*my')./sqrt((mx.*mx+my.*my).*(mx'.*mx'+my'.*my')+eps);
                mc = max(double_sin(cosn(:)));
                if mc > .3
                    curvM = mc;
                else
                    curvM = 0;
                end
                
                if numel(cosn2) > 0
                    ac = acos(cosn2(:));
                else
                    ac = acos(cosn(:));
                end
                
%                 ac(ac<1e-3) = [];
                if numel(ac) > 0
                    curvM2 = mean(ac);
                else
                    curvM2 = 0;
                end
            end
        end
    end
    if jtn < 2
        if sum(eb3(:)) == 0
            curvM = max(tb(eb2));
            curvM2 = max(ub(eb2));
        else
            curvM = max(tb(eb3));
            curvM2 = max(ub(eb3));
        end
    end
    
%     if num > 1
%         if sum(eb3(:)) == 0
%             curvM = .5*(curvM+max(tb(eb2)));
%             curvM2 = .4*(curvM2+max(ub(eb2)));
%         else
%             curvM = .5*(curvM+max(tb(eb3)));
%             curvM2 = .5*(curvM2+max(ub(eb3)));
%         end
%     end
    
    curvM = curvM*jtn2(i);
    
    if dn > 3 || jtn < 2
        curvM2 = curvM2*1.3+.1;
        curvM = curvM*1.3+.1;
    end
    
    if sum(eeb(:)) > 0
        curvM = curvM-2*min(tb(eeb));
        curvM2 = curvM2-2*min(ub(eeb));
        u = u+mean(ub(eb2))+mean(ub(eeb))-std(ub(eb2));
        t = t+mean(tb(eb2))+mean(tb(eeb))-std(tb(eb2));
    else
        u = u+mean(ub(eb2))-std(ub(eb2));
        t = t+mean(tb(eb2))-std(tb(eb2));
    end
    R3(i) = V(idx(i))*curvM;
    
    X(i,2:10) = [rr2, cc2, u, t, lns, nb*curvM, nb*curvM2, lns, c];
end

curv_u = X(:,4)*.33333;
curv_t = X(:,5)*.33333;
curv_n = X(:,10);

% mv = X(:,6);
streight = X(:,7);
circle = X(:,8);

% R3 = streight.*mag(idx);

streight(streight<0) = 0;
circle(circle<0) = 0;

Vx = (mag2(idx).*V(idx).*V(idx).*R2(idx)).^(1/3);

beta = (5e-5+.5*mean(Vx))./(Vx);
den = 3;
beta(beta>den) = den;
% beta(beta<1/den) = 1/den;

h = ones(2*minn+1);
h(minn-1:minn+3,minn-1:minn+3) = 0;

CC3 = filter2(h,CC);
G = (CC3+1);
gx2 = (G(idx) / mean(G(:)));

th1 = (.5+th*3);
th2 = sqrt(gx2.*beta./curv_n).*th1;


ix = (X(:,1) > 0 & (R3 > thr | (streight./curv_t > th1 & circle./curv_u > th1 & streight./curv_t > th2  & circle./curv_u > th2)));

idx = idx(ix);

[rn,cn] = ind2sub([m,n],idx);

cm = sqrt(streight.*circle./th2);
cm = cm(ix);

ix = false(m,n);
ix(idx) = 1;

end

function sinn2 = double_sin(cosn)
sinn2 = 2*sqrt(0.5*(1-cosn));%2*sqrt(1-cosn.*cosn)+;
end

function tm = texon_measure(ab,tex,tex5)
ab5 = sort(ab(tex5));
ab8 = sort(ab(tex),1);
dd = bsxfun(@minus, ab8, ab5);
dd = sort(mean(dd.*dd));
tm = max(dd);
% tm = mean(dd(end-3:end));
end
