function [metric, rr, cc] = robustcorner2(a, nonmax, boundary, th)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Changlin Wan
% E-mail: wancl@21cn.com
% Function: robustconrer2
% Version: 1.0
% Date: 2021/11/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isa(a, 'double')
    a = im2double(a);
end
a(isnan(a)) = 0;
if size(a,3) > 1
    a = rgb2gray(a);
end
[m1,n1] = size(a);

[aSmooth, dx, dy] = smoothGradient2(a,1.,0);

[noise_r, EDGE, curv, CURV] = noise_rate(dx,dy);

noise_r2 = noise_r;

noise_bound = .6;
fine_bound = .1;

if noise_r > noise_bound
    a = medfilt2(a,[3 3]);
end
nix = (aSmooth > 0 | dx > 0 | dy > 0);
alpha = sqrt(sum(nix(:)));

if alpha < 200
    bin = 1;
    gb = 1;
elseif alpha < 900
    bin = 2;
    if noise_r < .15
        gb = 1;
    else
        gb = .9+noise_r;
    end
else
    bin = round(alpha/300);
    if noise_r < .15
        gb = 1.5;
    else
        gb = 1.4+noise_r;
    end
end

if gb ~= 1
    [aSmooth,dx,dy] = smoothGradient2(a,gb,0);
    [noise_r2, EDGE, curv, CURV] = noise_rate(dx,dy);
end

a2 = a;
% if noise_r > .15
%     [aSmooth,dx,dy] = smoothGradient2(a,gb*(0.9+min(noise_bound,noise_r)),0);
%     [noise_r2, EDGE, curv, CURV] = noise_rate(dx,dy);
% end

minn = round(16*alpha/400);

if boundary == 1
    
    a1 = zeros(m1+14,n1+14);
    a1(8:m1+7,8:n1+7) = aSmooth;
    aSmooth = a1;

    a1(8:m1+7,8:n1+7) = a;
    a2 = a1;
    
    a1(8:m1+7,8:n1+7) = EDGE;
    EDGE = a1;
    
    a1(8:m1+7,8:n1+7) = curv;
    curv = a1;
    
    a1(8:m1+7,8:n1+7) = CURV;
    CURV = a1;
    
    a1(8:m1+7,8:n1+7) = dx;
    dx = a1;
    a1(8:m1+7,8:n1+7) = dy;
    dy = a1;
end


if noise_r < fine_bound
    EDGE2 = filter2(ones(7), EDGE) > 0;
else
    EDGE2 = filter2(ones(3), (curv>(.06+noise_r2/5)) & EDGE) > 0;
end

iEDGE = find(EDGE2>0);
idx = iEDGE;

if noise_r > fine_bound
    h1 = fspecial('gaussian',[7 7], 2);
else
    h1 = fspecial('gaussian',[5 5], 1);
end

[~, ~, a1, u, va, v] = invariant_gradient1(aSmooth, iEDGE, bin);
[radius, P, V, T] = turnradiusblock2(a1, dx, dy, u, va, v, iEDGE, bin);

thresh = mean(radius(idx)) - .3*std(radius(idx));
thresh = min(thresh,1e-4);

base =  mean(radius(idx))-.1*std(radius(idx));

thresh2 = min(2e-4,max(2e-5, base));
thresh3 = .15;

if th < 0
    th = min(.5,max((80*base/min(noise_bound,(noise_r*noise_r2)^.5)),.01));
end


if boundary == 1
    metric = radius(8:m1+7,8:n1+7);
else
    metric = radius;
end

if nonmax == 1
    
    R = filter2(h1,radius);
    R2 = R+2*abs(R-radius);
    
    w = min(9,max(3,2*round(alpha/300)+1));
%     w = 3;
    
    
    V2 = ordfilt2(V,9,true(3));
    P2 = ordfilt2(P,9,true(3));
    T2 = ordfilt2(T,9,true(3));
    
    
    if minn > 7
        [~,dx2,dy2] = smoothGradient2(a2,5,0);
        mag2 = sqrt(dx2.^2+dy2.^2);
        
   
        h = ones(2*minn+1);
        h(minn-1:minn+3,minn-1:minn+3) = 0;
        
        R3 = R2.*EDGE;
        G = filter2(h,R3);
        gx2 = G / mean(G(:));
        R2 = R2 ./ gx2;
    end
    Y = ordfilt2(R,w*w,true(w));
    
    [rr,cc] = find(T2 > thresh3 & V2 > thresh2 & R > thresh/3 & R2 > thresh & R == Y & CURV < max(.1,noise_r2));
    
    if th > 0 && minn > 7 %&& base > 1e-6
        [cm, rr, cc] = corner_curvature12(aSmooth, mag2, R2, rr, cc, dx, dy, EDGE, curv, V2, P2, T2, bin, th);
    end
    if boundary == 1
        rr = rr - 7;
        cc = cc - 7;
    end
else
    rr = [];
    cc = [];
end

end


function [aSmooth,dx,dy] = smoothGradient2(a, sigma, dt)

GaussianDieOff = .0001;

% Design the filters - a gaussian and its derivative
pw = 1:30; % possible widths
ssq = sigma*sigma;
width = max(find(exp(-(pw.*pw)/(2*sigma*sigma))>GaussianDieOff));
if isempty(width)
    width = 1;  % the user entered a really small sigma
end

t = (-width:width);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter
% Find the directional derivative of 2D Gaussian (along X-axis)
% Since the result is symmetric along X, we can get the derivative along
% Y-axis simply by transposing the result for X direction.
[x,y]=meshgrid(-width:width,-width:width);
dgau2D=-x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);

% H1=fspecial('gaussian', 2*width+1, sigma);
% H2=fspecial('gaussian', 2*width+1, 2);
%
% dgau2D = H1-H2;
%
% % %smooth the image out
aSmooth=imfilter(a,gau,'conv','replicate');   % run the filter accross rows
aSmooth=imfilter(aSmooth,gau','conv','replicate'); % and then accross columns

%apply directional derivatives
if dt == 0
    dx = imfilter(aSmooth, dgau2D, 'conv','replicate');
    dy = imfilter(aSmooth, dgau2D', 'conv','replicate');
elseif dt == 1
    [dx,dy] = derivative5(aSmooth,'x','y');
else
    [dx,dy] = derivative7(aSmooth,'x','y');
end

end
