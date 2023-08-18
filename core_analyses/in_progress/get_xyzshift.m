function [xyzshift,I_shifted] = get_xyzshift(A1, B1, toshift, pixel, rmax, xyzshift, distance)

% Gets and applies xyzshift for two images that we desire to register

% Based on "get_enrichment_3dMatrix_final.m" (author info below)
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

I_shifted = imtranslate(toshift,xyzshift);

end
