function [chr3] = local_connection(ua, iv, sz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Changlin Wan
% E-mail: wancl@21cn.com
% Function: local_connection
% Version: 1.0
% Date: 2021/11/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx,nn] = size(ua);

ar = zeros(nx,nn);
br = zeros(nx,nn);

u1 = ua;
u2 = ua;
th = 0;

for i=1:nn/2
    if i==1
        u3 = ua(:,[i+1:nn,1:i]);
        u4 = ua(:,[nn-i+1:nn,1:nn-i]);
        ixx1 = ua.*u4 >= th;
        ixx2 = ua.*u3 >= th;
%         ixx3 = ar == i-1;
%         ixx4 = br == i-1;
        ar(ixx2) = i;
        br(ixx1) = i;

        du = u1./(u1-u4);
        iix = ~ixx1;
        br(iix) = abs(du(iix));

        du = u2./(u2-u3);
        iix = ~ixx2;
        ar(iix) = abs(du(iix));

%         parts = sum(ar < 1,2)/2;
%         u1(:) = u4(:);
%         u2(:) = u3(:);
    else
        k =  i-1;
        c1 = [k+1:nn,1:k];
        c2 = [nn-k+1:nn,1:nn-k];
%         [r,c] = find(ar == k);
%         ixx = sub2ind([nx,nn], r, c);
%         ixx3 = sub2ind([nx,nn], r, c1(c));
        ar1 = ar(:,c1);
        ixx3 = ar == k;
        ar(ixx3) = k+ar1(ixx3);

%         [r,c] = find(br == k);
%         ixx = sub2ind([nx,nn], r, c);
%         ixx4 = sub2ind([nx,nn], r, c2(c));
        br1 = br(:,c2);
        ixx4 = br == k;
        br(ixx4) = k+br1(ixx4);
%         br(ixx4) = k+br(ixx4(:,[k:nn,1:k-1]));
    end
end

chr = ar+br;
chr2 = max(chr,[],2);
chr = chr2;
chr2(chr2 > nn) = nn;
% chr2 = min(chr,[],2);
chr2 = min(chr2,nn-chr2);
chr2 = min(chr2,nn/2-chr2);
chr2(chr>nn-1) = 2;
chr3 = zeros(sz);
chr3(iv) = chr2;
% 
% chr1 = zeros([numel(chr3),8]);
% chn1 = chr1;
% 
% for i=1:8
%     tt = chm(:,i);
%     chr1(iv,i) = tt(:);
%     tt = chi(:,i);
%     chn1(iv,i) = tt(:);
% end

end

