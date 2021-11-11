%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show corners into the output images or into the edge-image
function img1=mark(img,c,r,w)
[M,N,C]=size(img);
img1=img;
for i=1:numel(r)
    x = c(i);
    y = r(i);
%     if x <= w/2 || x > M-w/2 || y <= w/2 || y > N-w/2
%         continue;
%     end
if isa(img,'logical')
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<1);
    img1(max(1,x-floor(w/2)+1):min(M,x+floor(w/2)-1),max(1,y-floor(w/2)+1):min(N,y+floor(w/2)-1),:)=...
        img(max(1,x-floor(w/2)+1):min(M,x+floor(w/2)-1),max(1,y-floor(w/2)+1):min(N,y+floor(w/2)-1),:);
else
    H = double(max(img(:)));
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<.5*H)*H;
    img1(max(1,x-floor(w/2)+1):min(M,x+floor(w/2)-1),max(1,y-floor(w/2)+1):min(N,y+floor(w/2)-1),:)=...
        img(max(1,x-floor(w/2)+1):min(M,x+floor(w/2)-1),max(1,y-floor(w/2)+1):min(N,y+floor(w/2)-1),:);
end
end

