function [ eout, mag, ori, dx, dy, idxWeak, bone] = robustedge2(a, sigma, thresh)
%ROBUSTEDGE Summary of this function goes here
%   Detailed explanation goes here

% Transform to a double precision intensity image if necessary
if ~isa(a, 'double')
    a = im2double(a);
end

if size(a,3) > 1
    a = rgb2gray(a);
end


[m,n] = size(a);

% Calculate gradients using a derivative of Gaussian filter
[dx, dy] = smoothGradient(a, sigma);

% Calculate Magnitude of Gradient
magGrad = sqrt((dx.*dx) + (dy.*dy));
magGrad = magGrad / max(magGrad(:));

ori = mod(atan2(dy,dx),pi);

[lowThresh, highThresh] = selectThresholds(magGrad, 0.7, 1.0, 0.4);

idxLocalMax = [];
for dir=1:4
idxLocalMax = [idxLocalMax; cannyFindLocalMaxima2(dir, dx,dy, magGrad)];
end

idxWeak = idxLocalMax(magGrad(idxLocalMax) > lowThresh);

M1 = zeros(m,n);
M1(idxLocalMax)= magGrad(idxLocalMax);
magmax = max(M1(:));
if magmax > 0
    M1 = M1 / magmax;
end

CURV = strengthen(M1, idxWeak, dx, dy, 0);

K = .8*M1+.2*CURV;
idxWeak = idxLocalMax(K(idxLocalMax) > lowThresh);

magmax = max(K(:));
if magmax > 0
    K = K / magmax;
end

magGrad(idxLocalMax) = K(idxLocalMax);
idxLocalMax = [];
for dir=1:4
idxLocalMax = [idxLocalMax; cannyFindLocalMaxima2(dir, dx,dy, magGrad)];
end

[lowThresh, highThresh] = selectThresholds(magGrad, 0.7, 1.0, 0.4);

idxWeak = union(idxWeak, idxLocalMax(magGrad(idxLocalMax)>lowThresh));
% K = zeros(m,n);
% K(idxWeak) = magGrad(idxWeak);

[counts,x]=imhist(K(idxWeak), 64);
len = length(idxWeak)-length(idxLocalMax)*.3;
highThresh2 = find(cumsum(counts) > len*thresh(1),...
    1,'first') / 64;
lowThresh2 = find(cumsum(counts) > len*thresh(2),...
    1,'first') / 64;


E = false(m,n);
E(K > lowThresh2) = 1;

idxStrong = idxLocalMax(K(idxLocalMax) > highThresh2);
idx = strongfilt(idxStrong, m, n, 8);
bone = zeros(m,n);
bone(idx) = 1;
bone = bwmorph(bone,'dilate',2);
bone = bwmorph(bone,'erode',2);

% [rstrong,cstrong] = find(K>highThresh2);
rstrong = rem(idx-1, m)+1;
cstrong = floor((idx-1)/m)+1;
H = bwselect(E, cstrong, rstrong, 8);
eout = bwmorph(H, 'thin', 1);  % Thin double (or triple) pixel wide contours

mag = magGrad;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : cannyFindLocalMaxima
%
function idxLocalMax = cannyFindLocalMaxima2(direction,ix,iy,mag)
%
% This sub-function helps with the non-maximum suppression in the Canny
% edge detector.  The input parameters are:
%
%   direction - the index of which direction the gradient is pointing,
%               read from the diagram below. direction is 1, 2, 3, or 4.
%   ix        - input image filtered by derivative of gaussian along x
%   iy        - input image filtered by derivative of gaussian along y
%   mag       - the gradient magnitude image
%
%    there are 4 cases:
%
%                         The X marks the pixel in question, and each
%         3     2         of the quadrants for the gradient vector
%       O----0----0       fall into two cases, divided by the 45
%     4 |         | 1     degree line.  In one case the gradient
%       |         |       vector is more horizontal, and in the other
%       O    X    O       it is more vertical.  There are eight
%       |         |       divisions, but for the non-maximum suppression
%    (1)|         |(4)    we are only worried about 4 of them since we
%       O----O----O       use symmetric points about the center pixel.
%        (2)   (3)


[m,n] = size(mag);

% Find the indices of all points whose gradient (specified by the
% vector (ix,iy)) is going in the direction we're looking at.

switch direction
    case 1
        idx = find((iy<=0 & ix>-iy)  | (iy>=0 & ix<-iy));
    case 2
        idx = find((ix>0 & -iy>=ix)  | (ix<0 & -iy<=ix));
    case 3
        idx = find((ix<=0 & ix>iy) | (ix>=0 & ix<iy));
    case 4
        idx = find((iy<0 & ix<=iy) | (iy>0 & ix>=iy));
end

% Exclude the exterior pixels
if ~isempty(idx)
    v = mod(idx,m);
    extIdx = (v==1 | v==0 | idx<=m | (idx>(n-1)*m));
    idx(extIdx) = [];
end

ixv = ix(idx);
iyv = iy(idx);
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels
switch direction
    case 1
        d = abs(iyv./ixv);
        gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d;
    case 2
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d;
    case 3
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d;
    case 4
        d = abs(iyv./ixv);
        gradmag1 = mag(idx-m).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+m).*(1-d) + mag(idx+m+1).*d;
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : smoothGradient
%
function [GX, GY] = smoothGradient(I, sigma)

% Create an even-length 1-D separable Derivative of Gaussian filter

% Determine filter length
filterExtent = ceil(4*sigma);
x = -filterExtent:filterExtent;

% Create 1-D Gaussian Kernel
c = 1/(sqrt(2*pi)*sigma);
gaussKernel = c * exp(-(x.^2)/(2*sigma^2));

% Normalize to ensure kernel sums to one
gaussKernel = gaussKernel/sum(gaussKernel);

% Create 1-D Derivative of Gaussian Kernel
derivGaussKernel = gradient(gaussKernel);

% Normalize to ensure kernel sums to zero
negVals = derivGaussKernel < 0;
posVals = derivGaussKernel > 0;
derivGaussKernel(posVals) = derivGaussKernel(posVals)/sum(derivGaussKernel(posVals));
derivGaussKernel(negVals) = derivGaussKernel(negVals)/abs(sum(derivGaussKernel(negVals)));

% Compute smoothed numerical gradient of image I along x (horizontal)
% direction. GX corresponds to dG/dx, where G is the Gaussian Smoothed
% version of image I.
GX = imfilter(I, gaussKernel', 'conv', 'replicate');
GX = imfilter(GX, derivGaussKernel, 'conv', 'replicate');

% Compute smoothed numerical gradient of image I along y (vertical)
% direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
% version of image I.
GY = imfilter(I, gaussKernel, 'conv', 'replicate');
GY  = imfilter(GY, derivGaussKernel', 'conv', 'replicate');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : selectThresholds
%
function [lowThresh, highThresh] = selectThresholds(magGrad, PercentOfPixelsNotEdges, ht, lt)

[m,n] = size(magGrad);

counts=imhist(magGrad, 64);
th = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
    1,'first') / 64;
highThresh = ht*th;
lowThresh = lt*th;



function [CURV, idxWeak] = strengthen(M1, idxLocalMax, ax, ay, thresh)
[m,n] = size(M1);
a = zeros(m,n);
a(idxLocalMax) = 1;

w = 5;
for i=1:w
    a(i,:)=1;
    a(m+1-i,:)=1;
    a(:,i)=1;
    a(:,n+1-i)=1;
end
mn = zeros(m,n);

idxWeak = idxLocalMax(M1(idxLocalMax) > thresh);

% Exclude the exterior pixels
if ~isempty(idxWeak)
    v = mod(idxWeak,m);
    extIdx = (v<=w+1 | v>m-w | idxWeak<=w*m | (idxWeak>(n-w)*m));
    idxWeak(extIdx) = [];
end

mx = zeros(1,length(idxWeak));
for x=1:length(idxWeak)
    [rr,cc]=ind2sub([m,n],idxWeak(x));
    mr = rr + ay(rr,cc);
    mc = cc + ax(rr,cc);
    mo2 = ax(rr,cc)^2+ay(rr,cc)^2;
    weak1 = 0;
    weak2 = 0;
    strong = 0;
    
    for i=rr-w:rr+w
        for j=cc-w:cc+w
            if (M1(i,j) > 0)
                xo2 = (i - rr)*(i - rr) + (j - cc)*(j - cc);
                xm2 = (i - mr)*(i - mr) + (j - mc)*(j - mc);
                if(xo2>0)
                    cos = (mo2 + xo2 - xm2) / (2 * sqrt(mo2) * sqrt(xo2));
                else
                    cos = 0;
                end
                dist = sqrt(xo2) * cos;
                if(dist < 1 && dist > -1)
                    strong = strong + a(i,j);
                end
                if(dist > 1)
                    weak1 = weak1 + a(i,j);
                end
                if(dist < -1)
                    weak2 = weak2 + a(i,j);
                end
            end
        end
    end
    if(strong == (2*w+1)^2)
        mx(x) = -strong;
    else
        mx(x) = strong - weak1 - weak2;
        %mn(rr,cc) = strong - max(weak1, weak2);
        %mn(rr,cc) = min(strong, 2*w+1) - weak1 - weak2;
    end
end

% mn(idxWeak) = mx;
mn = mx;
mn = (mn - min(mn(:)));
magmax = max(mn(:));
if(magmax > 0)
    mn = mn / magmax;
end
% mn = 1.1*mn-.1*dense;
CURV = zeros(m,n);
% CURV(idxWeak) = mn(idxWeak);
CURV(idxWeak) = mn;


function idx = strongfilt(idxStrong, m, n, hl)
bw = zeros(m,n);
bw(idxStrong)=1;
CC = bwconncomp(bw,8);
area = cellfun(@numel, CC.PixelIdxList);
idx = CC.PixelIdxList(area > hl);
idx = vertcat(idx{:});

