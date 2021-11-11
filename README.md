# IMGLAB
Some Matlab Code of Image Processing
1. RobustCorner2 -- Corner Detection

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


