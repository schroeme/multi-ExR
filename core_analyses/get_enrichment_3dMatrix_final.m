function [Raa,Rab,Rba,Rbb] = get_enrichment_3dMatrix_final(A1, B1, pixel, rmax, xyzshift, distance, step, flag)

% function [xyzshift, d, 0, g] = get_enrichment_3dMatrix(A, B, pixel, rmax, xyzshift,distance,flag)
% calculates enrichment between two 3d Matrix
%       For two clusters supposed to be overlaped, set xyzshift=[0,0,0];
%       For two clusters not overlaped, set xyzshift=[] if the best
%       overlaping vector is unknow. Use distance=[min, max] to limit the
%       |xyzshift|. The function will find the vector xyzshift to get the
%       best overlap between min(A,Amax/mv) and min(B,Bmax/mv). 
%
% INPUTS
% A = matrix of channel A 
% B = matrix of channel  B
% pixel = voxel size of the matrix
% rmax = maximum r value to correlate in units of pixels. default is 100;
% xyzshift = xyz vector to overlap the two clusters.
% distance = [min, max], distance variation allowed for xyzshift. for RIM
%     and PSD95, we use [40,160]. This also incorporate any potential error
%     for dual channel registration.
% flag = display flag.  insert 1 to display errorbar(r, g, dg) after
%   computation. 
%
% OUTPUTS
% xyzshift = [dx, dy, dz]
% d = |xyzshift|
% r = radius values
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 3D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Updated 01.26.10 by Sarah Veatch.
% Modified 02.18.14 by Aihui Tang.
% Last updated 10.21.17 by Aihui Tang.

if nargin~=7, flag = 0; end  % flag for display
n = 7;
k7 = zeros(n,n,n);
for i = 1:n
    for j = 1:n
        for k = 1:n
            if sqrt((i-(n+1)/2)^2+(j-(n+1)/2)^2+(k-(n+1)/2)^2)<=n/2
                k7(i,j,k) = 1;
            end
        end
    end
end
k7 = k7/sum(k7(:));

bka = median(A1(:));
bkb = median(B1(:));
k = [0,1,0;1,1,1;0,1,0]/5;
A1 = max(0,convn(A1-bka,k,'same')); 
B1 = max(0,convn(B1-bkb,k,'same'));

A2 = max(0,A1-max(A1(:))/2); B2 = max(0,B1-max(B1(:))/2); 
A2 = convn(A2,k7,'same'); B2 = convn(B2,k7,'same');
maskA = double(A2>0); maskB = double(B2>0); 
A2 = A1.*maskA; B2 = B1.*maskB;
A3 = min(A2,max(A2(:))/4); B3 = min(B2,max(B2(:))/4); 

Na = sum(sum(sum(A2.*maskA)));  % number of particles within channel 1
Nb = sum(sum(sum(B2.*maskB)));  % number of particles within channel 2
Va = sum(sum(sum(maskA)));      % volume of maskA
Vb = sum(sum(sum(maskB)));      % volume of maskB
%NP = real(fftshift(ifftn(fftn(maskA,size(A2)+rmax).*conj(fftn(maskB,size(B2)+rmax)))))*1.06;
NP = real(fftshift(ifftn(fftn(maskA,size(A2)+rmax).*conj(fftn(maskB,size(B2)+rmax)))));
G1 = Va*Vb/Na/Nb*real(fftshift(ifftn(fftn(A2.*maskA,size(A2)+rmax).*conj(fftn(B2.*maskB,size(B2)+rmax)))))./NP; % 2D G(r) with proper normalization
L0=size(G1);L=floor(L0/2+1);
NPmask = max(0,max(max(max(NP)))/4-NP);
G1(find(NPmask))=0;             %set non-relavent points to 0

md = round(160/pixel);   %detect peak of correlation only in area with d<160 nm to center
if isempty(xyzshift) 
    G0 = real(fftshift(ifftn(fftn(A3.*maskA,size(A3)+rmax).*conj(fftn(B3.*maskB,size(B3)+rmax))))); 
    G0(find(NPmask))=0;
    G = G0*0; 
    pG = find(G==0);
    [x,y,z]=ind2sub(size(G), pG);
    d = pdist2([x, y, z],L);
    temp = find(d > distance(1)/pixel & d < distance(2)/pixel);
    G(pG(temp)) = 1;
    G = G.*G0; 
    Gm = max(max(max(G)));
    [mx, my, mz] = ind2sub(size(G), find(G==Gm));         %find the peak of correlation
    xyzshift = [mx(1)-L(1),my(1)-L(2),mz(1)-L(3)]; 
end

[mxa, mya, mza] = ind2sub(size(A2), find(A2==max(A2(:))));    
mxa = median(mxa); mya = median(mya); mza = median(mza); %peak position
[mxb, myb, mzb] = ind2sub(size(B2), find(B2==max(B2(:))));    
mxb = median(mxb); myb = median(myb); mzb = median(mzb); %peak position

%A2=A1;B2=B1;

raa = zeros(rmax/step+1,1);raa = raa./raa;
rab = raa; rba = raa; rbb = raa; raa1 = raa; rab1 = raa; rba1 = raa; rbb1 = raa; 

avgA2 = Na/Va;
avgB2 = Nb/Vb;
avgMaskA = maskA.*avgA2;
avgMaskB = maskB.*avgB2;

[xa,ya,za] = ind2sub(size(A2), find(A2>0));
[xa1,ya1,za1] = ind2sub(size(avgMaskA), find(avgMaskA>0));
d = pdist2([xa,ya,za],[mxa, mya, mza]).*pixel;
d1 = pdist2([xa1,ya1,za1],[mxa,mya,mza]).*pixel;
radius = 0:step:rmax;
[bincounts, ind] = histc(d, radius);
[bincounts1, ind1] = histc(d1, radius);
for jj = 1:max(ind)
    temp = find(ind == jj);F = [];
    if ~isempty(temp) 
        coor = [xa(temp),ya(temp),za(temp)];
        for jjj = 1:size(coor,1) F(jjj) = A2(coor(jjj,1),coor(jjj,2),coor(jjj,3)); end
        raa(jj) = mean(F); % A1 density along distance to A1 peak
    end
end
for jj = 1:max(ind1)
    temp1 = find(ind1 == jj); F1 = [];
    if ~isempty(temp1)
        coor1 = [xa1(temp1),ya1(temp1),za1(temp1)];
        for iii = 1:size(coor1,1) F1(iii) = avgMaskA(coor1(iii,1),coor1(iii,2),coor1(iii,3)); end
        raa1(jj) = mean(F1);
    end
end
raan = raa./raa1;

[xb,yb,zb] = ind2sub(size(B2), find(B2>0));
[xb1,yb1,zb1] = ind2sub(size(avgMaskB), find(avgMaskB>0));
d = pdist2([xb,yb,zb],[mxa, mya, mza]-xyzshift).*pixel;
d1 = pdist2([xb1,yb1,zb1],[mxa, mya, mza]-xyzshift).*pixel;
radius = 0:step:rmax;
[bincounts, ind] = histc(d, radius);
[bincounts1, ind1] = histc(d1, radius);
for jj = 1:max(ind)
    temp = find(ind == jj); F = [];
    if ~isempty(temp) 
        coor = [xb(temp),yb(temp),zb(temp)];
        for jjj = 1:size(coor,1) F(jjj) = B2(coor(jjj,1),coor(jjj,2),coor(jjj,3)); end
        rba(jj) = mean(F); % B1 density along distance to A1 peak
    end
end
for jj = 1:max(ind1)
    temp1 = find(ind1 == jj); F1 = [];
    if ~isempty(temp1)
        coor1 = [xb1(temp1),yb1(temp1),zb1(temp1)];
        for iii = 1:size(coor1,1) F1(iii) = avgMaskB(coor1(iii,1),coor1(iii,2),coor1(iii,3)); end
        rba1(jj) = mean(F1);
    end
end
rban = rba./rba1;

[xa,ya,za] = ind2sub(size(A2), find(A2>0));
[xa1,ya1,za1] = ind2sub(size(avgMaskA), find(avgMaskA>0));
d = pdist2([xa,ya,za],[mxb, myb, mzb]+xyzshift).*pixel;
d1 = pdist2([xa1,ya1,za1],[mxb, myb, mzb]+xyzshift).*pixel;
radius = 0:step:rmax;
[bincounts, ind] = histc(d, radius);
[bincounts1, ind1] = histc(d1, radius);
for jj = 1:max(ind)
    temp = find(ind == jj);F = [];
    if ~isempty(temp) 
        coor = [xa(temp),ya(temp),za(temp)];
        for jjj = 1:size(coor,1) F(jjj) = A2(coor(jjj,1),coor(jjj,2),coor(jjj,3)); end
        rab(jj) = mean(F); % A1 density along distance to B1 peak
    end
end
for jj = 1:max(ind1)
    temp1 = find(ind1 == jj); F1 = [];
    if ~isempty(temp1)
        coor1 = [xa1(temp1),ya1(temp1),za1(temp1)];
        for iii = 1:size(coor1,1) F1(iii) = avgMaskA(coor1(iii,1),coor1(iii,2),coor1(iii,3)); end
        rab1(jj) = mean(F1);
    end
end
rabn = rab./rab1;

[xb,yb,zb] = ind2sub(size(B2), find(B2>0));
[xb1,yb1,zb1] = ind2sub(size(avgMaskB), find(avgMaskB>0));
d = pdist2([xb,yb,zb],[mxb, myb, mzb]).*pixel;
d1 = pdist2([xb1,yb1,zb1],[mxb, myb, mzb]).*pixel;
radius = 0:step:rmax;
[bincounts, ind] = histc(d, radius);
[bincounts1, ind1] = histc(d1, radius);
for jj = 1:max(ind)
    temp = find(ind == jj);F = [];
    if ~isempty(temp) 
        coor = [xb(temp),yb(temp),zb(temp)];
        for jjj = 1:size(coor,1) F(jjj) = B2(coor(jjj,1),coor(jjj,2),coor(jjj,3)); end
        rbb(jj) = mean(F); % B1 density along distance to B1 peak
    end
end
for jj = 1:max(ind1)
    temp1 = find(ind1 == jj); F1 = [];
    if ~isempty(temp1)
        coor1 = [xb1(temp1),yb1(temp1),zb1(temp1)];
        for iii = 1:size(coor1,1) F1(iii) = avgMaskB(coor1(iii,1),coor1(iii,2),coor1(iii,3)); end
        rbb1(jj) = mean(F1);
    end
end
rbbn = rbb./rbb1;

if flag,
    figure('Color', 'white'); hold on
    plot(radius+2.5, raa,'b--')
    plot(radius+2.5, rba,'b')
    plot(radius+2.5, rab,'r')
    plot(radius+2.5, rbb,'r:')
    axis tight
end

%Rab = rabn; Rab = rab; Rab1 = rab1;%Rbb = rbbn;
Raa = raan; Rab = rabn; Rba = rban; Rbb = rbbn;
end
