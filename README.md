# IMGLAB
Some Matlab Codes of Image Processing (Demo Site http://www.92io.cn)
1. RobustCorner2 -- Feature Fusion based Corner Detection

I = imread('dataset\lab.bmp');
tic;
[~,r2,c2] = robustcorner2(I,1,0,-1);
toc;
figure(1),imshow(mark(I,r2,c2,4));hold on;plot(c2,r2,'ro');hold off;

2. RobustEdge2 -- Geometry Enhanced Edge Detection

I = imread('dataset\lab.bmp');
tic;
E = robustedge2(I,1,[.7 .4]);
toc;
figure(2),imshow(E);


