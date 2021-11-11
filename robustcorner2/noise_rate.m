function [rate, idxWeak, CURV, curv, M1, highThresh] = noise_rate(ax, ay)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Changlin Wan
% E-mail: wancl@21cn.com
% Function: noise_rate
% Version: 1.0
% Date: 2021/11/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(ax);
w = 2;
CURV = zeros(m,n);

mag = sqrt(ax.*ax+ay.*ay)+eps;
ax = ax./mag;
ay = ay./mag;

magmax = max(mag(:));
if magmax>0
    mag = mag / magmax;   % normalize
end
mag(isnan(mag)) = 0;

[lowThresh, highThresh] = selectThresholds(mag, .7, 1.0, .7);
% lowThresh = .2;
% highThresh = .35;

% h = getDoGKernel(7,1.5,.75);
% mag = filter2(h,mag);

% e = cannyFindLocalMaxima(ax,ay,mag,lowThresh);
% sum(e(:)>0)

idxLocalMax = [];
for dir=1:4
    idxLocalMax = [idxLocalMax; cannyFindLocalMaxima2(dir, ax,ay, mag)];
end

% M1 = zeros(m,n);
% idxWeak = find(e);
% M1(idxWeak) = mag(idxWeak);

M1 = zeros(m,n);
M1(idxLocalMax)= mag(idxLocalMax);
idxWeak = idxLocalMax(M1(idxLocalMax)>lowThresh);


% Exclude the exterior pixels
if ~isempty(idxWeak)
    v = mod(idxWeak,m);
    extIdx = (v<=w+1 | v>m-w | idxWeak<=w*m | (idxWeak>(n-w)*m));
    idxWeak(extIdx) = [];
end

theta = zeros(m,n);
theta(idxWeak) = atan(ay(idxWeak)./(ax(idxWeak)+eps));

[rr,cc]=ind2sub([m,n],idxWeak);

ab = zeros(numel(rr),1);
N = zeros(numel(rr),1);
theta1 = theta(idxWeak);

for i=-w:+w
    for j=-w:+w
        rrr = rr+i;
        ccc = cc+j;
        ix = sub2ind([m,n],rrr,ccc);
        
        angle = abs(theta1-theta(ix));
        angle(M1(ix) == 0) = 0;
        angle(angle > pi/2) = pi - angle(angle>pi/2);
        N = N + (M1(ix) > 0);
        ab = ab+angle;
    end
end
ix = N > 1;
CURV(idxWeak(ix)) = ab(ix)./(N(ix)-1);
CURV(isnan(CURV)) = 0;


rate = mean(CURV(idxWeak));

h = ones(15);
h(5:11,5:11) = 0;
curv = filter2(h/81, CURV);

% D = curv / max(curv(:));
% M1 = 1.2*M1 -  .2*D;
% [lowThresh, highThresh] = selectThresholds(mag(m/4:3*m/4,n/4:3*n/4), .99, .35, .2);

M2 = M1 > highThresh;

idxWeak = bwmorph(M2,'thin',1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : selectThresholds
%
function [lowThresh, highThresh] = selectThresholds(magGrad, PercentOfPixelsNotEdges, ht, lt)

[m,n] = size(magGrad);

counts=imhist(magGrad, 64);

th = find(cumsum(counts) > (PercentOfPixelsNotEdges)*(m*n),...
    1,'first') / 64;
highThresh = ht*th;
lowThresh = lt*th;

end

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
    otherwise
        idx = [];
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
    otherwise
        gradmag1 = [];
        gradmag2 = [];
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2);
end
